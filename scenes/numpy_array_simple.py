#
# Simple example scene for a simulation with some numpy array usage
#
from manta import *
import numpy as np

res = 64
gs = vec3(res,res,1)
s = Solver(name='main', gridSize = gs, dim=2)

# prepare grids
flags    = s.create(FlagGrid)
vel      = s.create(MACGrid)
density  = s.create(RealGrid)
pressure = s.create(RealGrid)
tmp      = s.create(RealGrid)

bWidth=1
flags.initDomain(boundaryWidth=bWidth) 
flags.fillGrid()
setOpenBound(flags, bWidth,'yY',FlagOutflow|FlagEmpty) 
source = s.create(Cylinder, center=gs*vec3(0.5,0.1,0.5), radius=res*0.14, z=gs*vec3(0, 0.02, 0))

npArray = np.ones( [res,res], dtype=np.float32 )

if (GUI):
	gui = Gui()
	gui.show( True ) 
	
#main loop
for t in range(400):
	mantaMsg('\nFrame %i' % (s.frame))
	source.applyToGrid(grid=density, value=1)
	advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2) 
	advectSemiLagrange(flags=flags, vel=vel, grid=vel,     order=2, openBounds=True, boundaryWidth=bWidth)
	resetOutflow(flags=flags,real=density) 

	setWallBcs(flags=flags, vel=vel)    
	addBuoyancy(density=density, vel=vel, gravity=vec3(0,-4e-3,0), flags=flags) 
	solvePressure(flags=flags, vel=vel, pressure=pressure)

	# small example function in test.cpp
	numpyTest( density, npArray, 0.01 ) # just adds constant value everywhere

	# grid conversion from numpyconvert.cpp plugins
	copyArrayToGridReal( target=tmp, source=npArray )
	
	s.step()

