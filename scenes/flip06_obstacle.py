#
# This FLIP example combines narrow band flip, 2nd order wall boundary conditions, and 
# adaptive time stepping.
# 
from manta import *

dim    = 3
res    = 64
#res    = 124
gs     = vec3(res,res,res)
if (dim==2):
	gs.z=1
s = Solver(name='main', gridSize = gs, dim=dim)

narrowBand    = 3
minParticles  = pow(2,dim)
saveParts     = False
frames        = 200

# Adaptive time stepping
s.frameLength = 0.8                 # length of one frame (in "world time")
s.cfl         = 3.0                 # maximal velocity per cell and timestep, 3 is fairly strict
s.timestep    = s.frameLength 
s.timestepMin = s.frameLength / 4.  # time step range
s.timestepMax = s.frameLength * 4.

# prepare grids and particles
flags     = s.create(FlagGrid)
phi       = s.create(LevelsetGrid)
phiParts  = s.create(LevelsetGrid)
phiObs    = s.create(LevelsetGrid)

vel       = s.create(MACGrid)
velOld    = s.create(MACGrid)
velParts  = s.create(MACGrid)
#mapWeights= s.create(MACGrid)

pressure  = s.create(RealGrid)
fractions = s.create(MACGrid)
tmpVec3   = s.create(VecGrid)

pp       = s.create(BasicParticleSystem) 
pVel     = pp.create(PdataVec3) 
mesh     = s.create(Mesh)

# acceleration data for particle nbs
pindex = s.create(ParticleIndexSystem) 
gpi    = s.create(IntGrid)

# scene setup
bWidth=1
flags.initDomain(boundaryWidth=bWidth, phiWalls=phiObs )
fluidVel = 0
fluidSetVel = 0
phi.setConst(999.)

# standing dam
fluidbox1 = Box( parent=s, p0=gs*vec3(0,0,0), p1=gs*vec3(1.0,0.3,1)) 
phi.join( fluidbox1.computeLevelset() )
fluidbox2 = Box( parent=s, p0=gs*vec3(0.1,0,0), p1=gs*vec3(0.2,0.75,1)) 
phi.join( fluidbox2.computeLevelset() )

if 1:
	sphere = Sphere( parent=s , center=gs*vec3(0.66,0.3,0.5), radius=res*0.2)
	phiObs.join( sphere.computeLevelset() )
	#obsbox = Box( parent=s, p0=gs*vec3(0.4,0.2,0), p1=gs*vec3(0.7,0.4,1))
	#obsbox = Box( parent=s, p0=gs*vec3(0.3,0.2,0), p1=gs*vec3(0.7,0.6,1))
	#phiObs.join( obsbox.computeLevelset() )


flags.updateFromLevelset(phi)
phi.subtract( phiObs );
sampleLevelsetWithParticles( phi=phi, flags=flags, parts=pp, discretization=2, randomness=0.05 )

if fluidVel!=0:
	# set initial velocity
	fluidVel.applyToGrid( grid=vel , value=fluidSetVel )
	mapGridToPartsVec3(source=vel, parts=pp, target=pVel )

# also sets boundary flags for phiObs
updateFractions( flags=flags, phiObs=phiObs, fractions=fractions, boundaryWidth=bWidth )
setObstacleFlags(flags=flags, phiObs=phiObs, fractions=fractions)

lastFrame = -1
if 1 and (GUI):
	gui = Gui()
	gui.show()
	#gui.pause()

# save reference any grid, to automatically determine grid size
if saveParts:
	pressure.save( 'ref_flipParts_0000.uni' );

#main loop
while s.frame < frames:
	maxVel = vel.getMax()
	s.adaptTimestep( maxVel )
	mantaMsg('\nFrame %i, time-step size %f' % (s.frame, s.timestep))
	
	# FLIP 
	pp.advectInGrid(flags=flags, vel=vel, integrationMode=IntRK4, deleteInObstacle=False, stopInObstacle=False )
	pushOutofObs( parts=pp, flags=flags, phiObs=phiObs )

	advectSemiLagrange(flags=flags, vel=vel, grid=phi, order=1) # first order is usually enough
	advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2)

	# create level set of particles
	gridParticleIndex( parts=pp , flags=flags, indexSys=pindex, index=gpi )
	unionParticleLevelset( pp, pindex, flags, gpi, phiParts )

	# combine level set of particles with grid level set
	phi.addConst(1.); # shrink slightly
	phi.join( phiParts );
	extrapolateLsSimple(phi=phi, distance=narrowBand+2, inside=True ) 
	extrapolateLsSimple(phi=phi, distance=3 )
	phi.setBoundNeumann(0) # make sure no particles are placed at outer boundary, warning - larger values can delete thin sheets at outer walls...
	flags.updateFromLevelset(phi)

	# combine particles velocities with advected grid velocities
	mapPartsToMAC(vel=velParts, flags=flags, velOld=velOld, parts=pp, partVel=pVel, weight=tmpVec3)
	extrapolateMACFromWeight( vel=velParts , distance=2, weight=tmpVec3 )
	combineGridVel(vel=velParts, weight=tmpVec3 , combineVel=vel, phi=phi, narrowBand=(narrowBand-1), thresh=0)
	velOld.copyFrom(vel)

	# forces & pressure solve
	addGravity(flags=flags, vel=vel, gravity=(0,-0.001,0))

	extrapolateMACSimple( flags=flags, vel=vel , distance=2, intoObs=True )
	setWallBcs(flags=flags, vel=vel, fractions=fractions, phiObs=phiObs)	

	solvePressure(flags=flags, vel=vel, pressure=pressure, phi=phi, fractions=fractions )

	extrapolateMACSimple( flags=flags, vel=vel , distance=4, intoObs=True )
	setWallBcs(flags=flags, vel=vel, fractions=fractions, phiObs=phiObs)

	if (dim==3):
		# mis-use phiParts as temp grid to close the mesh
		phiParts.copyFrom(phi)
		phiParts.setBound(0.5,0)
		phiParts.createMesh(mesh)

	# set source grids for resampling, used in adjustNumber!
	pVel.setSource( vel, isMAC=True )
	adjustNumber( parts=pp, vel=vel, flags=flags, minParticles=1*minParticles, maxParticles=2*minParticles, phi=phi, exclude=phiObs, narrowBand=narrowBand ) 
	flipVelocityUpdate(vel=vel, velOld=velOld, flags=flags, parts=pp, partVel=pVel, flipRatio=0.97 )
	
	s.step()

	if (lastFrame!=s.frame):
		# generate data for flip03_gen.py surface generation scene
		if saveParts:
			pp.save( 'flipParts_%04d.uni' % s.frame ); 
		if 0 and (GUI):
			gui.screenshot( 'flip06_%04d.png' % s.frame );

	#s.printMemInfo()
	lastFrame = s.frame;


