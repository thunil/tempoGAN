#
# Simple test for free-surface simulation with ghost fluid boundaries
# 

import sys
from manta import *
from helperInclude import *

# solver params
dim    = 3
res    = 52
frames = 50

if getVisualSetting():
	# in visual mode
	res    = 76  * getVisualSetting()
	frames = 100 * getVisualSetting()

gs = vec3(res,res,res)
if (dim==2):
	gs.z=1
s = Solver(name='main', gridSize = gs, dim=dim)
s.timestep = 0.25
accuracy = 5e-5
if getFloatSetting()==2:
	accuracy = 1e-10

if getVisualSetting():
	s.timestep = 0.28 / getVisualSetting()
	
# prepare grids and particles
flags     = s.create(FlagGrid)
vel       = s.create(MACGrid)
pressure  = s.create(RealGrid)
tmp       = s.create(RealGrid)

# scene setup
flags.initDomain(boundaryWidth=0)
basin = s.create(Box, p0=gs*vec3(0,0,0), p1=gs*vec3(1,0.2,1))
drop  = s.create(Sphere, center=gs*vec3(0.5,0.5,0.5), radius=res*0.15)
phi = basin.computeLevelset()
phi.join(drop.computeLevelset())
flags.updateFromLevelset(phi)

if 0 and (GUI):
	gui = Gui()
	gui.show()
	#gui.pause()
	

#main loop
for t in range(frames):
	
	# update and advect levelset
	phi.reinitMarching(flags=flags, velTransport=vel) #, ignoreWalls=False)
	advectSemiLagrange(flags=flags, vel=vel, grid=phi, order=2, clampMode=1)
	flags.updateFromLevelset(phi)
	
	# velocity self-advection
	advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2, clampMode=1)
	addGravity(flags=flags, vel=vel, gravity=vec3(0,-0.025,0))
	
	# pressure solve
	setWallBcs(flags=flags, vel=vel)
	solvePressure(flags=flags, vel=vel, pressure=pressure, cgMaxIterFac=0.5, cgAccuracy=accuracy, \
				  phi=phi ) # leave gfClamp at default
	setWallBcs(flags=flags, vel=vel)
	
	s.step()
	#gui.pause()
	#gui.screenshot( 'screenOn_%04d.png' % t );

	if 1 and getVisualSetting() and (t%getVisualSetting()==0):
		tmp.copyFrom( phi )
		tmp.multConst( -1 ); tmp.clamp(0,1.0)
		projectPpmFull( tmp, '%s_%04d.ppm' % (sys.argv[0],t/getVisualSetting()) , 1, 1.0 );


# check final state
doTestGrid( sys.argv[0],"phi"  , s, phi  , threshold=1e-07 , thresholdStrict=1e-10   )
doTestGrid( sys.argv[0],"vel"  , s, vel  , threshold=1e-07 , thresholdStrict=1e-10   )

