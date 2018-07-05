#
# Lid driven cavity example with viscosity
#
from manta import *
#setDebugLevel(2)

# the normalized unit cube in manta has which world space size?
worldScale = 1.0  

# viscosity, in [m^2/s] , rescale to unit cube
# uncomment one of these to select LDC with specific Reynolds nr
# (higher ones will need larger resolution!)
#visc       = 0.0002  / (worldScale*worldScale)  # Re 5k
visc       = 0.0001  / (worldScale*worldScale)  # Re 10k
#visc       = 0.00005 / (worldScale*worldScale)  # Re 20k 
#visc       = 0.00001 / (worldScale*worldScale)  # Re 100k 
#visc       = 0. # off, rely on numerical viscosity, no proper LDC!
# move whole top side in one time unit, ie 1 m/s if domain is 1m 
lidVel     = 1.00
# reynolds nr , characteristic length = domain size 

# to reduce the start up time, enable this to start the sim with a larger CFL number
doQuickstart = True

# compute Reynolds nr
Re = 0.
if visc>0.:
	Re = lidVel * worldScale / visc

# solver params
res  = 100  
gDim = 2
gs = vec3(res,res,res)
if gDim==2: gs.z = 1;
s = Solver(name='main', gridSize = gs, dim=gDim)

s.frameLength = 0.1 
s.timestepMin = s.frameLength * 0.01
s.timestepMax = s.frameLength * 1.0
s.cfl         = 1.0
s.timestep    = s.frameLength

if doQuickstart:
	# cfl reduction, start high, reduce over time
	s.cfl         = 10.0 
	mantaMsg("Note - quickstart active, starting with high CFL number, reduced later on",0)

density  = s.create(RealGrid)
# prepare grids
diff1 = s.create(RealGrid)
flags = s.create(FlagGrid)

flags.initDomain(boundaryWidth=1) 
flags.fillGrid()

vel      = s.create(MACGrid)
pressure = s.create(RealGrid)

timings = Timings()

if 1 and (GUI):
	gui = Gui()
	gui.show( True ) 
	gui.pause()


# geometry
lid = s.create(Box, p0=gs*vec3(0.0,1.0-(1./float(gs.x)*3.1),0.0), p1=gs*vec3(1.0,1.0,1.0)  )
# densities only for visual debugging... no influence!
source = s.create(Cylinder, center=gs*vec3(0.5,0.5,0.5), radius=res*0.10, z=gs*vec3(0, 0.10, 0))


	
#main loop
lastFrame = -1
outcnt    = 0
for t in range(98765):
	maxvel = vel.getMax()
	s.adaptTimestep(maxvel)
	mantaMsg('\nFrame %i, max vel %f, t %f, dt %f ' % (s.frame, maxvel, s.timeTotal, s.timestep ))

	if doQuickstart:
		# slowly reduce cfl number
		if s.cfl>5.0 and s.timeTotal>20.:
			s.cfl = 5.0
		if s.cfl>1.0 and s.timeTotal>30.:
			s.cfl = 1.0

	lid.applyToGrid(grid=vel, value=Vec3( lidVel*float(gs.x),0,0) )
	#if t%50==1:
	if (lastFrame!=s.frame) and (s.frame%25==0):
		source.applyToGrid(grid=density, value=1)

	advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2, clampMode=2) 
	advectSemiLagrange(flags=flags, vel=vel, grid=vel,     order=2, clampMode=2)
	resetOutflow(flags=flags,real=density) 

	# vel diffusion / viscosity!
	if visc>0.:
		# diffusion param for solve = const * dt / dx^2
		alphaV = visc * s.timestep * float(res*res)

		mantaMsg("Viscosity: %f , alpha=%f , Re=%f " %(visc, alphaV, Re), 0 )

		setWallBcs(flags=flags, vel=vel)    
		cgSolveDiffusion( flags, vel, alphaV )

	# solve pressure
	setWallBcs(flags=flags, vel=vel)    
	solvePressure(flags=flags, vel=vel, pressure=pressure)

	if 0 and (GUI) and (lastFrame!=s.frame) and (s.frame%1==0):
		gui.screenshot( 'ldc_re%05d_%04d.jpg' % (int(Re), s.frame) );
	
	lastFrame = s.frame
	s.step()

