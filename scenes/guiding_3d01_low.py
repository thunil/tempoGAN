#
# Simulation of a 3D buoyant smoke density plume at low resolution,
# to be used as guiding velocities for the high resolution version
#

from manta import *

# solver params
res0 = 40
scale = 1.0
res = res0*scale
gs = vec3(res,2.0*res,res)
s = Solver(name='main', gridSize = gs, dim=3)
s.timestep = 0.65*scale
numFrames = 200
timings = Timings()

output_uni = 'plume3DLowRes_%04d.uni'
output_ppm = 'plume3DLowRes_%04d.ppm'

# prepare grids
flags = s.create(FlagGrid)
vel = s.create(MACGrid)
velT = s.create(MACGrid)
density = s.create(RealGrid)
pressure = s.create(RealGrid)

# noise field
noise = s.create(NoiseField, loadFromFile=True)
noise.posScale = vec3(0)
noise.clamp = True
noise.clampNeg = 0
noise.clampPos = 1
noise.valScale = 1
noise.valOffset = 0.75
noise.timeAnim = 0.2

bWidth=0
flags.initDomain(boundaryWidth=bWidth) 
flags.fillGrid()
setOpenBound(flags, bWidth, 'yY', FlagOutflow|FlagEmpty) 

if (GUI):
	gui = Gui()
	gui.show()

source = s.create(Cylinder, center=gs*vec3(0.5,0.05,0.5), radius=res*0.1, z=gs*vec3(0, 0.02, 0))
	
#main loop
for t in range(int(numFrames*scale)):
	densityInflow(flags=flags, density=density, noise=noise, shape=source, scale=1, sigma=0.5)
		
	advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2)    
	advectSemiLagrange(flags=flags, vel=vel, grid=vel,     order=2, openBounds=True, boundaryWidth=bWidth+1)
	resetOutflow(flags=flags,real=density) 
	
	setWallBcs(flags=flags, vel=vel)
	addBuoyancy(density=density, vel=vel, gravity=vec3(0,-1e-3*scale,0), flags=flags)
	
	solvePressure(flags=flags, vel=vel, pressure=pressure)

	setWallBcs(flags=flags, vel=vel)
	#projectPpmFull( density, output_ppm % (t+1) , 0, 2.0 );
	vel.save( output_uni % (t) )
	
	s.step()

