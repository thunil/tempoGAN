#
# Obstacle test case
# 

import sys
from manta import *
from helperInclude import *

gs = vec3(31,47,33)
s = Solver(name='main', gridSize = gs)
#s.timestep = 1.3
s.timestep = 0.7

flags = s.create(FlagGrid)
vel = s.create(MACGrid)
density = s.create(RealGrid)
pressure = s.create(RealGrid)

flags.initDomain()
flags.fillGrid()

source = s.create(Box, p0=gs*vec3(0.3,0.1,0.3), p1=gs*vec3(0.7,0.2,0.7) )

obstacle1 = s.create(Box, p0=gs*vec3(0.5,0.5,0.5), p1=gs*vec3(0.8,0.6,0.8) )
obstacle2 = s.create(Box, p0=gs*vec3(0.0,0.8,0.0), p1=gs*vec3(0.4,0.9,0.4) )

obstacle1.applyToGrid( grid=flags, value=FlagObstacle)
obstacle2.applyToGrid( grid=flags, value=FlagObstacle)

if 0 and (GUI):
    gui = Gui()
    gui.show(); gui.pause()
    
#main loop
for t in range(10):
    source.applyToGrid( grid=density, value=3.72)
        
    advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2, clampMode=1)    
    advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2, clampMode=1)
    
    setWallBcs(flags=flags, vel=vel)    
    addBuoyancy(density=density, vel=vel, gravity=vec3(0,-5e-2,0), flags=flags)
    
    solvePressure(flags=flags, vel=vel, pressure=pressure)
    setWallBcs(flags=flags, vel=vel)
    s.step()


# check final state
doTestGrid( sys.argv[0],"dens" , s, density , threshold=0.002 , thresholdStrict=1e-10, debugShowDifference=False )
doTestGrid( sys.argv[0],"vel"  , s, vel     , threshold=0.003 , thresholdStrict=1e-10, debugShowDifference=False )

for t in range(99): # pause for debugging
	s.step()

