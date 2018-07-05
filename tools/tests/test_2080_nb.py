#
# Narrow band test
# 
import sys
from manta import *
from helperInclude import *

# Choose dimension (2 or 3) and resolution
dim = 3
res = 44

narrowBandWidth  = 3
combineBandWidth = narrowBandWidth - 1

# Solver params
gs = vec3(res,res,res)
if dim==2: gs.z = 1;

s = Solver(name='main', gridSize = gs, dim=dim)
s.timestep = 0.9

gravity = (0,-0.003,0)

minParticles = pow(2,dim)

# Prepare grids and particles
flags    = s.create(FlagGrid)

phiParts = s.create(LevelsetGrid) 
phi      = s.create(LevelsetGrid) 
pressure = s.create(RealGrid)

vel      = s.create(MACGrid)
velOld   = s.create(MACGrid)
velParts = s.create(MACGrid)
mapWeights  = s.create(MACGrid)

pp        = s.create(BasicParticleSystem) 
pVel      = pp.create(PdataVec3)
mesh      = s.create(Mesh)

# Acceleration data for particle nbs
pindex = s.create(ParticleIndexSystem)
gpi    = s.create(IntGrid)

# Geometry in world units (to be converted to grid space upon init)
flags.initDomain(boundaryWidth=0)
phi.initFromFlags(flags)

# A simple breaking dam
fluidBasin = s.create(Box, p0=gs*vec3(0,0,0), p1=gs*vec3(1.0,0.15,1.0))
phi.join( fluidBasin.computeLevelset() )

fluidDam = s.create(Box, p0=gs*vec3(0,0.15,0), p1=gs*vec3(0.4,0.5,0.8))
phi.join( fluidDam.computeLevelset() )
	
flags.updateFromLevelset(phi)

sampleLevelsetWithParticles( phi=phi, flags=flags, parts=pp, discretization=2, randomness=0.4 )
mapGridToPartsVec3(source=vel, parts=pp, target=pVel )

if 0 and (GUI):
	gui = Gui()
	gui.show( dim==3 )
	#gui.pause()
		
# Main loop
while s.frame < 10:

	# Advect particles and grid phi
	pp.advectInGrid(flags=flags, vel=vel, integrationMode=IntRK4, deleteInObstacle=False ) 
	advectSemiLagrange(flags=flags, vel=vel, grid=phi, order=1)
	flags.updateFromLevelset(phi)
	
	advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2, clampMode=1)

	# Create level set of particles
	gridParticleIndex( parts=pp , flags=flags, indexSys=pindex, index=gpi )
	unionParticleLevelset( pp, pindex, flags, gpi, phiParts , radiusFactor=1 )

	# Combine level set of particles with grid level set
	phi.addConst(1.); # shrink slightly
	phi.join( phiParts );
	extrapolateLsSimple(phi=phi, distance=narrowBandWidth+2, inside=True )

	extrapolateLsSimple(phi=phi, distance=3 )
	flags.updateFromLevelset(phi)

	# Combine particles velocities with advected grid velocities
	mapPartsToMAC(vel=velParts, flags=flags, velOld=velOld, parts=pp, partVel=pVel, weight=mapWeights)
	extrapolateMACFromWeight( vel=velParts , distance=2, weight=mapWeights )
	combineGridVel(vel=velParts, weight=mapWeights , combineVel=vel, phi=phi, narrowBand=combineBandWidth, thresh=0)
	velOld.copyFrom(vel)
		
	# Forces & pressure solve
	addGravity(flags=flags, vel=vel, gravity=gravity)
	setWallBcs(flags=flags, vel=vel)
	solvePressure(flags=flags, vel=vel, pressure=pressure, phi=phi)
	setWallBcs(flags=flags, vel=vel)

	# Extrapolate velocities
	extrapolateMACSimple( flags=flags, vel=vel, distance=5 )
	
	# Update particle velocities
	flipVelocityUpdate(vel=vel, velOld=velOld, flags=flags, parts=pp, partVel=pVel, flipRatio=0.95 )
		
	# Resample particles
	pVel.setSource( vel, isMAC=True ) # Set source grids for resampling, used in adjustNumber!
	adjustNumber( parts=pp, vel=vel, flags=flags, minParticles=1*minParticles, maxParticles=2*minParticles, phi=phi, narrowBand=narrowBandWidth ) 

	#if dim==3: phi.createMesh(mesh) 
	s.step()
		
doTestGrid( sys.argv[0],"phi"     , s, phi     , threshold=1e-07 , thresholdStrict=1e-10  )
doTestGrid( sys.argv[0],"vel"     , s, vel     , threshold=1e-07 , thresholdStrict=1e-10  )

doTestGrid( sys.argv[0],"phiParts", s, phiParts, threshold=1e-07 , thresholdStrict=1e-10  )
doTestGrid( sys.argv[0],"velParts", s, velParts, threshold=1e-07 , thresholdStrict=1e-10  )
	
