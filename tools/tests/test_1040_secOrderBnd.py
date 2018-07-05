#
# Second order boundaries, inner sphere test
#
import sys
from manta import *
from helperInclude import *

new_BC     = True
dim        = 2
res        = 16
gs         = vec3(res,res,res)
if (dim==2): gs.z = 1
s          = FluidSolver(name='main', gridSize = gs, dim=dim)
s.timestep = 1
timings    = Timings()

flags     = s.create(FlagGrid)
vel       = s.create(MACGrid)
pressure  = s.create(RealGrid)
fractions = s.create(MACGrid)
density   = s.create(RealGrid)

flags.initDomain()

center = gs*vec3(0.5,0.5,0.5)
radius = res*0.4
sphere = s.create(Sphere, center=center, radius=radius)
phiObs = sphere.computeLevelset()
phiObs.multConst(-1)

initVortexVelocity(phiObs=phiObs, vel=vel, center=center, radius=radius)

updateFractions( flags=flags, phiObs=phiObs, fractions=fractions)
setObstacleFlags(flags=flags, phiObs=phiObs, fractions=fractions)
flags.fillGrid()

if 0 and (GUI):
	gui = Gui()
	gui.show()
	gui.pause()

#main loop
for t in range(10):

	advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2, orderSpace=1, clampMode=1)  
	advectSemiLagrange(flags=flags, vel=vel, grid=vel    , order=2, strength=1.0, clampMode=1)

	if(new_BC):
		setWallBcs(flags=flags, vel=vel, fractions=fractions, phiObs=phiObs)
		extrapolateMACSimple( flags=flags, vel=vel, distance=1 );

		solvePressure( flags=flags, vel=vel, pressure=pressure, fractions=fractions)

		setWallBcs(flags=flags, vel=vel, fractions=fractions, phiObs=phiObs)
		extrapolateMACSimple( flags=flags, vel=vel, distance=1 );

	else:
		setWallBcs(flags=flags, vel=vel)
		solvePressure( flags=flags, vel=vel, pressure=pressure )
		setWallBcs(flags=flags, vel=vel)
	
	#timings.display()
	s.step()


# check final state
doTestGrid( sys.argv[0],"frac" , s, fractions , threshold=0.0001 , thresholdStrict=1e-10 )
doTestGrid( sys.argv[0],"vel"  , s, vel       , threshold=0.0001 , thresholdStrict=1e-10 )

