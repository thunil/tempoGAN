#
# Lid driven cavity test
#
import sys
from manta import *
from helperInclude import *

visc       = 0.0001  
lidVel     = 1.00

# solver params
res  = 50  
gDim = 2
gs = vec3(res,res,res)
if gDim==2: gs.z = 1;
s = Solver(name='main', gridSize = gs, dim=gDim)

s.frameLength = 0.1 
s.timestepMin = s.frameLength * 0.01
s.timestepMax = s.frameLength * 1.0
s.cfl         = 1.0
s.timestep    = s.frameLength

density  = s.create(RealGrid)
# prepare grids
diff1 = s.create(RealGrid)
flags = s.create(FlagGrid)

flags.initDomain(boundaryWidth=1) 
flags.fillGrid()

vel      = s.create(MACGrid)
pressure = s.create(RealGrid)

timings = Timings()

if 0 and (GUI):
	gui = Gui()
	gui.show( True ) 
	gui.pause()


lid = s.create(Box, p0=gs*vec3(0.0,1.0-(1./float(gs.x)*3.1),0.0), p1=gs*vec3(1.0,1.0,1.0)  )
source = s.create(Cylinder, center=gs*vec3(0.5,0.5,0.5), radius=res*0.10, z=gs*vec3(0, 0.10, 0))
	
#main loop
lastFrame = -1
outcnt    = 0
for t in range(50):
	maxvel = vel.getMax()
	s.adaptTimestep(maxvel)
	#mantaMsg('\nFrame %i, max vel %f, t %f, dt %f ' % (s.frame, maxvel, s.timeTotal, s.timestep ))

	lid.applyToGrid(grid=vel, value=Vec3( lidVel*float(gs.x),0,0) )

	if (lastFrame!=s.frame) and (s.frame%25==0):
		source.applyToGrid(grid=density, value=1)

	advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2, clampMode=2) 
	advectSemiLagrange(flags=flags, vel=vel, grid=vel,     order=2, clampMode=2)
	resetOutflow(flags=flags,real=density) 

	setWallBcs(flags=flags, vel=vel)    
	density.setBound(0.0, 1)

	if visc>0.: 
		alphaV = visc * s.timestep * float(res*res) 
		cgSolveDiffusion( flags, vel, alphaV )

	setWallBcs(flags=flags, vel=vel)    
	solvePressure(flags=flags, vel=vel, pressure=pressure)

	if 0 and (GUI) and (lastFrame!=s.frame) and (s.frame%1==0):
		gui.screenshot( 'ldc06re%05d_%04d.jpg' % (int(Re), s.frame) );
	
	lastFrame = s.frame
	s.step()

# check final state
doTestGrid( sys.argv[0],"vel"  , s, vel     , threshold=0.0005 , thresholdStrict=1e-08 )

