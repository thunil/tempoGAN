#
# 2d flip test
# 
import sys
from manta import *
from helperInclude import *

# solver params
dim = 2
particleNumber = 2
res = 64
gs = vec3(res,res,res)
if (dim==2):
	gs.z=1
	particleNumber = 3 # use more particles in 2d
s = Solver(name='main', gridSize = gs, dim=dim)
s.timestep = 0.7

# prepare grids and particles
flags    = s.create(FlagGrid)
vel      = s.create(MACGrid)
velOld   = s.create(MACGrid)
pressure = s.create(RealGrid)
tmpVec3  = s.create(VecGrid)
dens     = s.create(RealGrid)

pp       = s.create(BasicParticleSystem) 
# add velocity data to particles
pVel     = pp.create(PdataVec3) 
pDens    = pp.create(PdataReal) 

# scene setup
flags.initDomain(boundaryWidth=0)
fluidbox = s.create(Box, p0=gs*vec3(0.1,0,0), p1=gs*vec3(0.4,0.6,1)) # breaking dam
phiInit = fluidbox.computeLevelset()
flags.updateFromLevelset(phiInit)
# phiInit is not needed from now on!

# note, there's no resamplig here, so we need _LOTS_ of particles...
sampleFlagsWithParticles( flags=flags, parts=pp, discretization=particleNumber, randomness=0.2 )
pDens.setConst( 0.5 )

if 0 and (GUI):
    gui = Gui()
    gui.show()
    #gui.pause()
    
#main loop
for t in range(60):
    
    # FLIP 
    pp.advectInGrid(flags=flags, vel=vel, integrationMode=IntRK4, deleteInObstacle=False ) 
    mapPartsToMAC(vel=vel, flags=flags, velOld=velOld, parts=pp, partVel=pVel, weight=tmpVec3 ) 
    extrapolateMACFromWeight( vel=vel , distance=2, weight=tmpVec3 ) 
    markFluidCells( parts=pp, flags=flags )

    # density marker for testing
    mapPartsToGrid(target=dens, flags=flags, parts=pp, source=pDens )

    addGravity(flags=flags, vel=vel, gravity=(0,-0.003,0))

    # pressure solve
    setWallBcs(flags=flags, vel=vel)    
    solvePressure(flags=flags, vel=vel, pressure=pressure)
    setWallBcs(flags=flags, vel=vel)

    # we dont have any levelset, ie no extrapolation, so make sure the velocities are valid
    extrapolateMACSimple( flags=flags, vel=vel )
    
    # FLIP velocity update
    flipVelocityUpdate(vel=vel, velOld=velOld, flags=flags, parts=pp, partVel=pVel, flipRatio=0.97 )
    
    #gui.screenshot( 'flipt_%04d.png' % t );
    s.step()


doTestGrid( sys.argv[0],"dens" , s, dens  , threshold=0.0001 , thresholdStrict=1e-10 )
doTestGrid( sys.argv[0],"vel"  , s, vel   , threshold=0.001  , thresholdStrict=1e-10 )

