#
# Surface turbulence test
# 
import sys
from manta import *
from helperInclude import *

# solver params
dim = 3
res = 18 
gs = vec3(res,res,res)
if (dim==2):
	gs.z=1
s = Solver(name='main', gridSize = gs, dim=dim)
s.timestep = 0.8
minParticles = pow(2,dim)

# size of particles 
radiusFactor = 1.0

# prepare grids and particles
flags    = s.create(FlagGrid)
phi      = s.create(LevelsetGrid)

vel      = s.create(MACGrid)
velOld   = s.create(MACGrid)
pressure = s.create(RealGrid)
tmpVec3  = s.create(VecGrid)
tmpReal  = s.create(RealGrid)

surfacePointsDisplaced = s.create(BasicParticleSystem)
spdDummy = surfacePointsDisplaced.create(PdataVec3) # dummy for display

pp       = s.create(BasicParticleSystem)
pVel     = pp.create(PdataVec3)
pPrevPos  = pp.create(PdataVec3)

surfacePoints = s.create(BasicParticleSystem)
surfaceNormal = surfacePoints.create(PdataVec3)
surfaceWaveH = surfacePoints.create(PdataReal)
surfaceWaveDtH = surfacePoints.create(PdataReal)
surfaceWaveSource = surfacePoints.create(PdataReal)
surfaceWaveSeedAmplitude = surfacePoints.create(PdataReal)
surfaceWaveSeed = surfacePoints.create(PdataReal)

# acceleration data for particle nbs
pindex = s.create(ParticleIndexSystem) 
gpi    = s.create(IntGrid)

flags.initDomain(boundaryWidth=1)

# falling drop
fluidBasin = s.create(Box, p0=gs*vec3(0,0,0), p1=gs*vec3(1.0,0.2,1.0)) # basin
dropCenter = vec3(0.5,0.4,0.5)
dropRadius = 0.1
fluidDrop  = s.create(Sphere, center=gs*dropCenter, radius=res*dropRadius)
fluidVel   = s.create(Sphere, center=gs*dropCenter, radius=res*(dropRadius+0.05) )
phi = fluidBasin.computeLevelset()
phi.join( fluidDrop.computeLevelset() )

flags.updateFromLevelset(phi)
sampleLevelsetWithParticles( phi=phi, flags=flags, parts=pp, discretization=2, randomness=0.35 )

# testing only
spdDummy2= surfacePointsDisplaced.create(PdataReal) # dummy for test 
dummyFlags = s.create(FlagGrid) # also for test
dummyFlags.initDomain(boundaryWidth=1)

if 0 and (GUI):
	gui = Gui()
	gui.show()

	gui.nextPdata()
	gui.nextPartDisplay()
	gui.nextPartDisplay()
	#gui.pause()
   
   
#main loop
while s.frame<30:
	# FLIP 
	pp.advectInGrid(flags=flags, vel=vel, integrationMode=IntRK4, deleteInObstacle=False )

	# make sure we have velocities throught liquid region
	mapPartsToMAC(vel=vel, flags=flags, velOld=velOld, parts=pp, partVel=pVel, weight=tmpVec3 ) 
	extrapolateMACFromWeight( vel=vel , distance=2, weight=tmpVec3 )  # note, tmpVec3 could be free'd now...
	markFluidCells( parts=pp, flags=flags )

	# create approximate surface level set, resample particles
	gridParticleIndex( parts=pp , flags=flags, indexSys=pindex, index=gpi )
	unionParticleLevelset( pp, pindex, flags, gpi, phi , radiusFactor=1. ) 
	resetOutflow(flags=flags,parts=pp,index=gpi,indexSys=pindex) 
	# extend levelset somewhat, needed by particle resampling in adjustNumber
	extrapolateLsSimple(phi=phi, distance=4, inside=True); 

	# forces & pressure solve
	addGravity(flags=flags, vel=vel, gravity=(0,-0.001,0))
	setWallBcs(flags=flags, vel=vel)	
	solvePressure(flags=flags, vel=vel, pressure=pressure, phi=phi)
	setWallBcs(flags=flags, vel=vel)

	# set source grids for resampling, used in adjustNumber!
	pVel.setSource( vel, isMAC=True )
	adjustNumber( parts=pp, vel=vel, flags=flags, minParticles=1*minParticles, maxParticles=2*minParticles, phi=phi, radiusFactor=1. ) 

	# make sure we have proper velocities
	extrapolateMACSimple( flags=flags, vel=vel )
	
	flipVelocityUpdate(vel=vel, velOld=velOld, flags=flags, parts=pp, partVel=pVel, flipRatio=0.97 )
		
	# approx 1,000 coarse particles and 2,000 surface points (for 'breaking dam' setup)
	particleSurfaceTurbulence( flags=flags, coarseParts=pp, coarsePartsPrevPos=pPrevPos, surfPoints=surfacePoints, surfaceNormals=surfaceNormal, surfaceWaveH=surfaceWaveH, surfaceWaveDtH=surfaceWaveDtH, surfacePointsDisplaced=surfacePointsDisplaced, surfaceWaveSource=surfaceWaveSource, surfaceWaveSeed=surfaceWaveSeed, surfaceWaveSeedAmplitude=surfaceWaveSeedAmplitude, res=res,
		nbSurfaceMaintenanceIterations = 4,
		surfaceDensity = 15,
		dt = 0.005,
		waveSpeed = res,
		waveDamping = 0.1,
		waveSeedFrequency = 4.0,
		waveMaxAmplitude = 0.5,
		waveMaxFrequency = 128.0,
		waveSeedingCurvatureThresholdRegionCenter = 0.025,
		waveSeedingCurvatureThresholdRegionRadius = 0.01,
		waveSeedStepSizeRatioOfMax = 0.05 )

	# older debug checks
	#debugCheckParts(pp, flags)
	#debugCheckParts(surfacePoints, flags)
	#debugCheckParts(surfacePointsDisplaced, flags)
	
	# white color 
	spdDummy.setConst(vec3(1,1,1))

	# use temp grid to test displaced part pos
	spdDummy2.setConst(0.1)
	# map parts needs obstacle flags at the boundary in dummyFlags
	dummyFlags.setConst( FlagFluid )
	mapPartsToGrid(target=tmpReal, flags=dummyFlags, parts=surfacePointsDisplaced, source=spdDummy2 )
	# for sanity check, this should fail:
	#mapPartsToGrid(target=tmpReal, flags=dummyFlags, parts=surfacePoints, source=spdDummy2 )

	s.step()
	
		
doTestGrid( sys.argv[0],"phi" , s, phi     , threshold=1e-07 , thresholdStrict=1e-10  )
doTestGrid( sys.argv[0],"vel" , s, vel     , threshold=1e-07 , thresholdStrict=1e-10  )
# most important:
doTestGrid( sys.argv[0],"surf", s, tmpReal , threshold=1e-07 , thresholdStrict=1e-10  )

