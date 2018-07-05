#
# Simulation of a buoyant smoke with adaptive time-stepping
#

from manta import *

# solver params
dim = 3
res = 64
gs = vec3(res,1.5*res,res)
if (dim==2):
	gs.z=1
s = FluidSolver(name='main', gridSize = gs, dim=dim)

# how many frames to calculate 
frames    = 100

# set time step range
s.frameLength = 1.2   # length of one frame (in "world time")
s.timestepMin = 0.2   # time step range
s.timestepMax = 2.0
s.cfl         = 3.0   # maximal velocity per cell
s.timestep    = (s.timestepMax+s.timestepMin)*0.5

# prepare grids
flags = s.create(FlagGrid)
vel = s.create(MACGrid)
density = s.create(RealGrid)
pressure = s.create(RealGrid)

# noise field
noise = s.create(NoiseField, loadFromFile=True)
noise.posScale = vec3(45)
noise.clamp = True
noise.clampNeg = 0
noise.clampPos = 1
noise.valScale = 1
noise.valOffset = 0.75
noise.timeAnim = 0.2

flags.initDomain()
flags.fillGrid()
timings = Timings()

if (GUI):
	gui = Gui()
	gui.show( dim==2 )

source = s.create(Cylinder, center=gs*vec3(0.5,0.1,0.5), radius=res*0.14, z=gs*vec3(0, 0.02, 0))


#main loop
lastFrame = -1
while s.frame < frames:
	
	maxvel = vel.getMax()
	s.adaptTimestep(maxvel)
	mantaMsg('\nFrame %i, time-step size %f' % (s.frame, s.timestep))

	
	if s.timeTotal<50.:
		densityInflow(flags=flags, density=density, noise=noise, shape=source, scale=1, sigma=0.5)

	advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2)    
	advectSemiLagrange(flags=flags, vel=vel, grid=vel    , order=2)
	
	setWallBcs(flags=flags, vel=vel)    
	addBuoyancy(density=density, vel=vel, gravity=vec3(0,-6e-3,0), flags=flags)
	
	solvePressure( flags=flags, vel=vel, pressure=pressure )
	setWallBcs(flags=flags, vel=vel)

	if 0 and (GUI) and (lastFrame!=s.frame) and (s.frame%1==0):
		gui.screenshot( 'plumead_%04d.jpg' % s.frame );

	#timings.display()
	lastFrame = s.frame 
	s.step()

