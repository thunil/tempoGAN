#
# Flip scene with adaptive time stepping (otherwise similar to flip02)
# 
from manta import *

# solver params
dim = 3
res = 80
gs = vec3(res,res,res)
if (dim==2):
	gs.z=1
s = Solver(name='main', gridSize = gs, dim=dim)

# how many frames to calculate 
frames    = 200

# adaptive time stepping
s.frameLength = 0.6   # length of one frame (in "world time")
s.timestepMin = 0.1   # time step range
s.timestepMax = 2.0
s.cfl         = 1.5   # maximal velocity per cell
s.timestep    = (s.timestepMax+s.timestepMin)*0.5

minParticles = pow(2,dim)
timings = Timings()

# size of particles 
radiusFactor = 1.0

# prepare grids and particles
flags    = s.create(FlagGrid)
phi      = s.create(LevelsetGrid)

vel      = s.create(MACGrid)
velOld   = s.create(MACGrid)
pressure = s.create(RealGrid)
tmpVec3  = s.create(VecGrid)
tstGrid  = s.create(RealGrid)
phiObs   = s.create(LevelsetGrid)

pp       = s.create(BasicParticleSystem) 
pVel     = pp.create(PdataVec3) 
# test real value, not necessary for simulation
pTest    = pp.create(PdataReal) 
mesh     = s.create(Mesh)

# acceleration data for particle nbs
pindex = s.create(ParticleIndexSystem) 
gpi    = s.create(IntGrid)

# scene setup, 0=breaking dam, 1=drop into pool
# geometry in world units (to be converted to grid space upon init)
setup = 0
flags.initDomain(boundaryWidth=0)
fluidVel = 0
fluidSetVel = 0

if setup==0:
	# breaking dam
	fluidbox = Box( parent=s, p0=gs*vec3(0,0,0), p1=gs*vec3(0.4,0.6,1)) # breaking dam
	phi = fluidbox.computeLevelset()
	#fluidbox2 = Box( parent=s, p0=gs*vec3(0.2,0.7,0.3), p1=gs*vec3(0.5,0.8,0.6)) # breaking dam
	#phi.join( fluidbox2.computeLevelset() )
elif setup==1:
	# falling drop
	fluidBasin = Box( parent=s, p0=gs*vec3(0,0,0), p1=gs*vec3(1.0,0.2,1.0)) # basin
	dropCenter = vec3(0.5,0.5,0.5)
	dropRadius = 0.15
	fluidSetVel= vec3(0,-0.03,0)
	fluidDrop  = Sphere( parent=s , center=gs*dropCenter, radius=res*dropRadius)
	fluidVel   = Sphere( parent=s , center=gs*dropCenter, radius=res*(dropRadius+0.05) )
	phi = fluidBasin.computeLevelset()
	phi.join( fluidDrop.computeLevelset() )
	flags.updateFromLevelset(phi)

flags.updateFromLevelset(phi)

if dim==3:
	# obstacle init needs to go after updateFromLs
	obsBox = Box( parent=s, p0=gs*vec3(0.7,0.0,0.5), p1=gs*vec3(0.8,1.0,0.8)) 
	obsBox.applyToGrid(grid=flags, value=(FlagObstacle) )
	#obsBox.applyToGrid(grid=flags, value=(FlagObstacle|FlagStick) )

sampleLevelsetWithParticles( phi=phi, flags=flags, parts=pp, discretization=2, randomness=0.05 )
mapGridToPartsVec3(source=vel, parts=pp, target=pVel )

if fluidVel!=0:
	# set initial velocity
	fluidVel.applyToGrid( grid=vel , value=gs*fluidSetVel )
	mapGridToPartsVec3(source=vel, parts=pp, target=pVel )

# testing the real channel while resampling - original particles
# will have a value of 0.1, new particle will get a value from the tstGrid
testInitGridWithPos(tstGrid)
pTest.setConst( 0.1 )

lastFrame = -1
if 1 and (GUI):
	gui = Gui()
	gui.show( dim==2 )
	#gui.pause()
	  
	# show all particles shaded by velocity
	gui.nextPdata()
	gui.nextPartDisplay()
	gui.nextPartDisplay()



#main loop
while s.frame < frames:
	maxVel = vel.getMax()
	s.adaptTimestep( maxVel )
	mantaMsg('\nFrame %i, time-step size %f' % (s.frame, s.timestep))

	
	pp.advectInGrid(flags=flags, vel=vel, integrationMode=IntRK4, deleteInObstacle=False )

	# make sure we have velocities throught liquid region
	mapPartsToMAC(vel=vel, flags=flags, velOld=velOld, parts=pp, partVel=pVel, weight=tmpVec3 ) 
	extrapolateMACFromWeight( vel=vel , distance=2, weight=tmpVec3 )  # note, tmpVec3 could be free'd now...
	markFluidCells( parts=pp, flags=flags )

	# create approximate surface level set, resample particles
	gridParticleIndex( parts=pp , flags=flags, indexSys=pindex, index=gpi )
	unionParticleLevelset( pp, pindex, flags, gpi, phi , radiusFactor ) 
	extrapolateLsSimple(phi=phi, distance=4, inside=True); 
	# note - outside levelset doesnt matter...

	# forces & pressure solve
	addGravity(flags=flags, vel=vel, gravity=(0,-0.003,0))
	setWallBcs(flags=flags, vel=vel)	
	solvePressure(flags=flags, vel=vel, pressure=pressure, phi=phi)
	setWallBcs(flags=flags, vel=vel)

	# set source grids for resampling, used in adjustNumber!
	pVel.setSource( vel, isMAC=True )
	pTest.setSource( tstGrid );
	adjustNumber( parts=pp, vel=vel, flags=flags, minParticles=1*minParticles, maxParticles=2*minParticles, phi=phi, radiusFactor=radiusFactor ) 

	# make sure we have proper velocities
	extrapolateMACSimple( flags=flags, vel=vel, distance=(int(maxVel*1.5)+2) ) 
	
	flipVelocityUpdate(vel=vel, velOld=velOld, flags=flags, parts=pp, partVel=pVel, flipRatio=0.97 )

	if 0 and (dim==3):
		phi.createMesh(mesh)
	
	#timings.display()
	#s.printMemInfo()
	s.step()

	# optionally save particle data , or screenshot
	if 0 and (lastFrame!=s.frame):
		pp.save( 'flipParts_%04d.uni' % s.frame ); 
	if 0 and (GUI) and (lastFrame!=s.frame):
		gui.screenshot( 'flip04_%04d.png' % s.frame );

	lastFrame = s.frame;


