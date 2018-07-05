#
# Simple example scene (hello world)
# Simulation of a buoyant smoke density plume (with noise texture as smoke source)
#

#import pdb; pdb.set_trace()

from manta import *

# solver params
res = 64
gs  = vec3(res, int(1.5*res), res)
s   = FluidSolver(name='main', gridSize = gs)

# prepare grids
flags    = s.create(FlagGrid)
vel      = s.create(MACGrid)
density  = s.create(RealGrid)
pressure = s.create(RealGrid)

# noise field, tweak a bit for smoke source
noise = s.create(NoiseField, loadFromFile=True)
noise.posScale = vec3(45)
noise.clamp = True
noise.clampNeg = 0
noise.clampPos = 1
noise.valOffset = 0.75
noise.timeAnim = 0.2

source = s.create(Cylinder, center=gs*vec3(0.5,0.1,0.5), radius=res*0.14, z=gs*vec3(0, 0.02, 0))

flags.initDomain()
flags.fillGrid()

if (GUI):
	gui = Gui()
	gui.show()
	
#main loop
for t in range(250):
	mantaMsg('\nFrame %i' % (s.frame))
	if t<100:
		densityInflow(flags=flags, density=density, noise=noise, shape=source, scale=1, sigma=0.5)
		
	# optionally, enforce inflow velocity
	#source.applyToGrid(grid=vel, value=vec3(0.1,0,0))

	advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2)    
	advectSemiLagrange(flags=flags, vel=vel, grid=vel    , order=2, strength=1.0)
	
	setWallBcs(flags=flags, vel=vel)    
	addBuoyancy(density=density, vel=vel, gravity=vec3(0,-6e-4,0), flags=flags)
	
	solvePressure( flags=flags, vel=vel, pressure=pressure )
	s.step()

