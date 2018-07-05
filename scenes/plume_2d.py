#
# Simple example scene for a 2D simulation
# Simulation of a buoyant smoke density plume with open boundaries at top & bottom
#
from manta import *

# solver params
res = 64
gs = vec3(res,res,1)
s = Solver(name='main', gridSize = gs, dim=2)
s.timestep = 1.0
timings = Timings()

# prepare grids
flags = s.create(FlagGrid)
vel = s.create(MACGrid)
density = s.create(RealGrid)
pressure = s.create(RealGrid)

bWidth=1
flags.initDomain(boundaryWidth=bWidth) 
flags.fillGrid()

setOpenBound(flags, bWidth,'yY',FlagOutflow|FlagEmpty) 

if (GUI):
	gui = Gui()
	gui.show( True ) 
	#gui.pause()

source = s.create(Cylinder, center=gs*vec3(0.5,0.1,0.5), radius=res*0.14, z=gs*vec3(0, 0.02, 0))
	
#main loop
for t in range(400):
	mantaMsg('\nFrame %i' % (s.frame))

	if t<300:
		source.applyToGrid(grid=density, value=1)
		
	advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2) 
	advectSemiLagrange(flags=flags, vel=vel, grid=vel,     order=2, openBounds=True, boundaryWidth=bWidth)
	resetOutflow(flags=flags,real=density) 

	setWallBcs(flags=flags, vel=vel)    
	addBuoyancy(density=density, vel=vel, gravity=vec3(0,-4e-3,0), flags=flags)

	solvePressure(flags=flags, vel=vel, pressure=pressure)
	
	#timings.display()    
	s.step()

