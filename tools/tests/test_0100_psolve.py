#
# 3d pressure solve
# 

import sys
from manta import *
from helperInclude import *

# solver params
res = 52
gs = vec3(res,res,res)
s = Solver(name='main', gridSize = gs, dim=3)
s.timestep = 1.0

# prepare grids
flags    = s.create(FlagGrid)
vel      = s.create(MACGrid)
pressure = s.create(RealGrid)
dummy    = s.create(RealGrid)

flags.initDomain()
flags.fillGrid()

if 0 and (GUI):
	gui = Gui()
	gui.show()
	gui.pause()

velSource = s.create(Box, p0=gs*vec3(0.3,0.4,0.3), p1=gs*vec3(0.7,0.8,0.7) )

# ============================

# a simple solve
vel.setConst( vec3(0,0,0) )
velSource.applyToGrid(grid=vel, value=vec3(0.15, 0.3, 0.21) )    
solvePressure(flags=flags, vel=vel, pressure=pressure, cgMaxIterFac=99, cgAccuracy=1e-04, zeroPressureFixing=False)
s.step()

# check - note, unfortunately low threshold here necessary for float<>double comparisons...
doTestGrid( sys.argv[0], "pressure0" , s, pressure , threshold=1e-04, thresholdStrict=1e-10)
doTestGrid( sys.argv[0], "vel0"      , s, vel      , threshold=1e-04, thresholdStrict=1e-10)

# ============================
    
# second solve , with BCs
vel.setConst( vec3(0,0,0) )
velSource.applyToGrid(grid=vel, value=vec3(1.5, 3, 2.1) )

setWallBcs(flags=flags, vel=vel) 
solvePressure(flags=flags, vel=vel, pressure=pressure, cgMaxIterFac=99, cgAccuracy=1e-04, zeroPressureFixing=False)
s.step()

# check final state
doTestGrid( sys.argv[0], "pressure" , s, pressure , threshold=1e-04, thresholdStrict=1e-10)
doTestGrid( sys.argv[0], "vel"      , s, vel      , threshold=1e-04, thresholdStrict=1e-10)

