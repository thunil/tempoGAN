#
# Example scene with surface wave turbulence 
# (corresponding paper: Surface Turbulence for Particle-Based Liquid Simulations, Mercier et al., SIGGRAPH Asia 2015)
# 
from manta import *

# solver params
dim = 3
res = 32
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
#tstGrid  = s.create(RealGrid)

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

# scene setup, 0=breaking dam, 1=drop into pool
setup    = 0
bWidth   = 1
fluidVel = 0
fluidSetVel = 0
flags.initDomain(boundaryWidth=bWidth)

if setup==0:
	# breaking dam
	fluidbox = s.create(Box, p0=gs*vec3(0,0,0), p1=gs*vec3(0.4,0.4,1)) # breaking dam
	phi = fluidbox.computeLevelset()
elif setup==1:
	# falling drop
	fluidBasin = s.create(Box, p0=gs*vec3(0,0,0), p1=gs*vec3(1.0,0.1,1.0)) # basin
	dropCenter = vec3(0.5,0.3,0.5)
	dropRadius = 0.1
	fluidDrop  = s.create(Sphere, center=gs*dropCenter, radius=res*dropRadius)
	fluidVel   = s.create(Sphere, center=gs*dropCenter, radius=res*(dropRadius+0.05) )
	fluidSetVel= vec3(0,-1,0)
	phi = fluidBasin.computeLevelset()
	phi.join( fluidDrop.computeLevelset() )

flags.updateFromLevelset(phi)
sampleLevelsetWithParticles( phi=phi, flags=flags, parts=pp, discretization=2, randomness=0.35 )

if fluidVel!=0:
	# set initial velocity
	fluidVel.applyToGrid( grid=vel , value=fluidSetVel )
	mapGridToPartsVec3(source=vel, parts=pp, target=pVel )


if 1 and (GUI):
	gui = Gui()
	gui.show()

	gui.nextPdata()
	gui.nextPartDisplay()
	gui.nextPartDisplay()
	#gui.pause()
   
   
#main loop
for t in range(500):
	
	# FLIP 
	pp.advectInGrid(flags=flags, vel=vel, integrationMode=IntRK4, deleteInObstacle=False )

	# make sure we have velocities throught liquid region
	mapPartsToMAC(vel=vel, flags=flags, velOld=velOld, parts=pp, partVel=pVel, weight=tmpVec3 ) 
	extrapolateMACFromWeight( vel=vel , distance=2, weight=tmpVec3 )  # note, tmpVec3 could be free'd now...
	markFluidCells( parts=pp, flags=flags )

	# create approximate surface level set, resample particles
	gridParticleIndex( parts=pp , flags=flags, indexSys=pindex, index=gpi )
	unionParticleLevelset( pp, pindex, flags, gpi, phi , radiusFactor ) 
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
	adjustNumber( parts=pp, vel=vel, flags=flags, minParticles=1*minParticles, maxParticles=2*minParticles, phi=phi, radiusFactor=radiusFactor ) 

	# make sure we have proper velocities
	extrapolateMACSimple( flags=flags, vel=vel )
	
	flipVelocityUpdate(vel=vel, velOld=velOld, flags=flags, parts=pp, partVel=pVel, flipRatio=0.97 )
			
	particleSurfaceTurbulence( flags=flags, coarseParts=pp, coarsePartsPrevPos=pPrevPos, surfPoints=surfacePoints, surfaceNormals=surfaceNormal, surfaceWaveH=surfaceWaveH, surfaceWaveDtH=surfaceWaveDtH,surfacePointsDisplaced=surfacePointsDisplaced, surfaceWaveSource=surfaceWaveSource, surfaceWaveSeed=surfaceWaveSeed, surfaceWaveSeedAmplitude=surfaceWaveSeedAmplitude, res=res,
		nbSurfaceMaintenanceIterations = 6,
		surfaceDensity = 12,
		outerRadius = 1.0*radiusFactor,
		dt = 0.005,
		waveSpeed = 32, # res,
		waveDamping = 0.05,
		waveSeedFrequency = 4.0,# res/8.0,
		waveMaxAmplitude = 0.5, # res/64.0
		waveMaxSeedingAmplitude = 0.5, # as a multiple of max amplitude
		waveMaxFrequency = 128.0,
		waveSeedingCurvatureThresholdRegionCenter = 0.025, # any curvature higher than this value will seed waves
		waveSeedingCurvatureThresholdRegionRadius = 0.01,
		waveSeedStepSizeRatioOfMax = 0.05 # higher values will result in faster and more violent wave seeding
		)
	
	
	
	# white color for displaced points
	spdDummy.setConst(vec3(1,1,1))
	
	# output
	# surfacePointsDisplaced.save(  'output/surfaceTurbulence/res16/displacedPoints_%04d.uni' % t  );
	# gui.screenshot( 'output/surfaceTurbulence/img_%04d.png' % t );

	s.step()
		


