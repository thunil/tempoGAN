#
# Simple example for free-surface simulation, falling drop
# that shouldnt hit the floor and give straight down velocities

import sys
from manta import *
from helperInclude import *

# solver params
dim    = 3
res    = 45
frames = 18

if getVisualSetting():
	# in visual mode
	res    = 90  * getVisualSetting()
	frames = 100 * getVisualSetting()

gs = vec3(res,res,res)
if (dim==2):
	gs.z=1
s = Solver(name='main', gridSize = gs, dim=dim)
s.timestep = 0.6
accuracy = 5e-5

if getVisualSetting():
	s.timestep = 0.3 / getVisualSetting()

# prepare grids and particles
flags    = s.create(FlagGrid)
vel      = s.create(MACGrid)
pressure = s.create(RealGrid)
# temp grid for visual output
tmp      = s.create(RealGrid) 

# scene setup
flags.initDomain(boundaryWidth=0)
liqDrop = s.create(Box, p0=gs*vec3(0.4,0.75,0.4), p1=gs*vec3(0.6,0.95,0.6))
phi = liqDrop.computeLevelset()
flags.updateFromLevelset(phi)

if 0 and (GUI):
	gui = Gui()
	gui.show()
	gui.pause()
   
#main loop
for t in range(frames):
	
	# update and advect levelset
	phi.reinitMarching(flags=flags, velTransport=vel) #, ignoreWalls=False)
	advectSemiLagrange(flags=flags, vel=vel, grid=phi, order=2, clampMode=1)
	flags.updateFromLevelset(phi)
	
	# velocity self-advection
	advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2, clampMode=1)
	addGravity(flags=flags, vel=vel, gravity=vec3(0,-0.0125,0))
	
	# pressure solve
	setWallBcs(flags=flags, vel=vel)
	solvePressure(flags=flags, vel=vel, pressure=pressure, cgMaxIterFac=0.5, cgAccuracy=accuracy, phi=phi)
	setWallBcs(flags=flags, vel=vel)
	
	s.step()

	if 1 and getVisualSetting() and (t%getVisualSetting()==0):
		tmp.copyFrom( phi )
		tmp.multConst( -1 ); tmp.clamp(0,1.0)
		projectPpmFull( tmp, '%s_%04d.ppm' % (sys.argv[0],t/getVisualSetting()) , 1, 1.0 );


# check final state , eg 32bit and 64bit windows version can have slightly different results...
doTestGrid( sys.argv[0], "phi" , s, phi  , 1e-05  , thresholdStrict=1e-10, debugShowDifference=False )
doTestGrid( sys.argv[0], "vel" , s, vel  , 1e-05  , thresholdStrict=1e-10, debugShowDifference=False )

for t in range(99): # pause for debugging
	s.step()

