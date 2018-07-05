#
# Simple 2D buoyant plume
#
import sys
from manta import *
from helperInclude import *

# solver params
res      = 64
frames   = 30
accuracy = 1e-05

if getFloatSetting()==2:
	accuracy = 1e-10

if getVisualSetting():
	# in visual mode
	res    = 128 * getVisualSetting()
	frames = 200 * getVisualSetting()
	#frames = 5 # debug!

gs = vec3(res,res,1)
s = Solver(name='main', gridSize = gs, dim=2)
s.timestep = 1.0

if getVisualSetting():
	s.timestep = 0.5 / getVisualSetting()

# prepare grids
flags = s.create(FlagGrid)
vel = s.create(MACGrid)
density = s.create(RealGrid)
pressure = s.create(RealGrid)

flags.initDomain()
flags.fillGrid()

if 0 and (GUI):
	gui = Gui()
	gui.show()

source = s.create(Cylinder, center=gs*vec3(0.5,0.1,0.5), radius=res*0.14, z=gs*vec3(0, 0.02, 0))

#main loop
for t in range(frames): 
	source.applyToGrid(grid=density, value=1)

	advectSemiLagrange(flags=flags, vel=vel, grid=density, order=1) # note, first order here...
	advectSemiLagrange(flags=flags, vel=vel, grid=vel,     order=1)
	
	setWallBcs(flags=flags, vel=vel)	
	addBuoyancy(density=density, vel=vel, gravity=vec3(0,-9e-3,0), flags=flags)
	
	solvePressure(flags=flags, vel=vel, pressure=pressure, cgAccuracy=accuracy, cgMaxIterFac=5. )
	
	s.step()

	if 1 and getVisualSetting() and (t%(2*getVisualSetting())==0):
		projectPpmFull( density, '%s_%04d.ppm' % (sys.argv[0],t/(2*getVisualSetting()) ) , 0, 1.0 );


# check final state
doTestGrid( sys.argv[0],"dens" , s, density , threshold=0.0005 , thresholdStrict=1e-08 )
doTestGrid( sys.argv[0],"vel"  , s, vel     , threshold=0.0005 , thresholdStrict=1e-08 )

