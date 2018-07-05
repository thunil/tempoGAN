#
# Simple buoyant smoke density plume
#

import sys
from manta import *
from helperInclude import *

# solver params
res    = 60
frames = 15

if getVisualSetting():
	# in visual mode
	res    = 80 * getVisualSetting()
	frames = 75 * getVisualSetting()
	#frames = 3 # debug!

gs = vec3(res,1.25*res,res)
s = Solver(name='main', gridSize = gs)
s.timestep = 0.5

if getVisualSetting():
	s.timestep = 0.5 / getVisualSetting()

# prepare grids
flags = s.create(FlagGrid)
vel = s.create(MACGrid)
density = s.create(RealGrid)
pressure = s.create(RealGrid)

# noise field
noise = s.create(NoiseField, loadFromFile=True )
noise.posScale = vec3(45)
noise.clamp = True
noise.clampNeg = 0
noise.clampPos = 1
noise.valScale = 1
noise.valOffset = 0.75
noise.timeAnim = 0.2

flags.initDomain()
flags.fillGrid()

bWidth = 1
openB = 'xXyYzZ'
setOpenBound(flags,bWidth,openB,FlagOutflow|FlagEmpty) 

if 0 and (GUI):
	gui = Gui()
	gui.show()

source = s.create(Cylinder, center=gs*vec3(0.5,0.1,0.5), radius=res*0.14, z=gs*vec3(0, 0.02, 0))
    
#main loop
for t in range(frames): 
	densityInflow(flags=flags, density=density, noise=noise, shape=source, scale=1, sigma=0.5)
	
	#source.applyToGrid(grid=vel, value=velInflow)
	advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2, clampMode=1)	
	resetOutflow(flags=flags,real=density) 
	# note - this scene uses bWidth+1 , this is unnecessary, but doesnt make a big difference for this test
	advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2, openBounds=True, boundaryWidth=bWidth+1, clampMode=1)
	
	setWallBcs(flags=flags, vel=vel)	
	addBuoyancy(density=density, vel=vel, gravity=vec3(0,-5e-2,0), flags=flags)
	
	solvePressure(flags=flags, vel=vel, pressure=pressure)
	setWallBcs(flags=flags, vel=vel)
	#density.save('den%04d.uni' % t)
	
	s.step()

	if 1 and getVisualSetting() and (t%getVisualSetting()==0):
		projectPpmFull( density, '%s_%04d.ppm' % (sys.argv[0],t/getVisualSetting()) , 0, 4.0 );

# check final state (note - threshold are pretty large)
doTestGrid( sys.argv[0],"dens" , s, density , threshold=0.001 , thresholdStrict=1e-10 )
doTestGrid( sys.argv[0],"vel"  , s, vel     , threshold=0.005 , thresholdStrict=1e-10 )

