#
# Test particle save/load
# Unfortunately not directly, but through a helper grid
# 
import sys
from manta import *
from helperInclude import *

if getVisualSetting():
	# skip when in visual mode...
	exit(0);
	
# solver params
res = 50
gs = vec3(res,res,res)
s = Solver(name='main', gridSize = gs, dim=3)
s.timestep = 0.58

# prepare grids and particles
flags    = FlagGrid(parent = s)
vel      = MACGrid(parent = s)
velOld   = MACGrid(parent = s)
pressure = RealGrid(parent = s)
tmpVec3  = VecGrid(parent = s)
density  = RealGrid(parent = s)

pp       = BasicParticleSystem(parent = s) 
# add velocity data to particles
pVel     = pp.create(PdataVec3) 
pDens    = pp.create(PdataReal) 

# scene setup
flags.initDomain(boundaryWidth=0)

# noise field
noise = NoiseField(parent = s, loadFromFile=True )
noise.posScale = vec3(100) # high frequency
noise.clamp = True
noise.clampNeg = 0
noise.clampPos = 1.2
noise.valScale = 0.9
noise.valOffset = 0.15
noise.timeAnim = 0.1

	
if 0 and (GUI):
	gui = Gui()
	gui.show()
	gui.pause()

# generate particle data
genRefFiles = getGenRefFileSetting()
if (genRefFiles==1):
	# enable one of the following
	fluidbox1 = Box(parent=s, p0=gs*vec3(0.2,0.2,0.2), p1=gs*vec3(0.8,0.4,0.8))
	fluidbox2 = Box(parent=s, p0=gs*vec3(0.2,0.6,0.2), p1=gs*vec3(0.8,0.8,0.8))
	phiInit = fluidbox1.computeLevelset()
	phiInit.join( fluidbox2.computeLevelset() )
	flags.updateFromLevelset(phiInit)
	# phiInit is not needed from now on!
	del phiInit

	sampleFlagsWithParticles( flags=flags, parts=pp, discretization=3, randomness=0.2 )
	pDens.setConst( 1.3 ) # for buoyancy

	flags.fillGrid()
	mapPartsToGrid(target=density, flags=flags, parts=pp, source=pDens )

	# one solve step, generate static velocity field
	addBuoyancy(density=density, vel=vel, gravity=vec3(0,-5e-1,0), flags=flags)
	setWallBcs(flags=flags, vel=vel)    
	solvePressure(flags=flags, vel=vel, pressure=pressure)
	setWallBcs(flags=flags, vel=vel)

	# now set pattern to particle density
	setNoisePdata( pp, pDens, noise )
		
	#main loop
	for t in range(5):
		# simply move around particles a bit 
		pp.advectInGrid(flags=flags, vel=vel, integrationMode=IntRK4, deleteInObstacle=False ) 
		s.step()

	density.setConst(-1.)
	mapPartsToGrid(target=density, flags=flags, parts=pp, source=pDens ) 

	pp.save(    referenceFilename(sys.argv[0],"parts" ) );
	pDens.save( referenceFilename(sys.argv[0],"pDens" ) );
else:
	flags.fillGrid()

	# only test restoring a particle system here - no sim!
	pp.load(    referenceFilename(sys.argv[0],"parts" ) );
	pDens.load( referenceFilename(sys.argv[0],"pDens" ) );

	mapPartsToGrid(target=density, flags=flags, parts=pp, source=pDens ) 

	s.step()

# check resulting values, note the strict/double prec threshold is ridiculously un-strict - problem is the float rounding in the uni/raw files, which makes this test pretty meaningless for doubles...
doTestGrid( sys.argv[0],"dens" , s, density  , threshold=1e-05 , thresholdStrict=1e-02 )


