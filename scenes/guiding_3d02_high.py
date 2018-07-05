#
# Simulation of a 3D buoyant smoke density plume 
# guided by a low res sim 
#
# parameterss:
#   W = guiding weight (constant value of wScalar)
#   beta = blur radius
#
from manta import *

# params from original simulation
# res1 = before interp, res2 = after interp
timestep = 0.65
res1 = 40
numFrames = 200
factor = 2
res2 = int(res1*factor)

# solver params
gs2 = vec3(res2,int(2.0*res2),res2)
s2 = Solver(name='main', gridSize = gs2, dim=3)
s2.timestep = timestep
timings = Timings()

input_uni  = 'plume3DLowRes_%04d.uni'
output_ppm = 'plume3DHighRes_%04d.ppm'
output_uni = 'plume3DHighRes_%04d.uni'

# PD params
beta = 5
wScalar = 2

tau = 0.58/wScalar
sigma = 2.44/tau
theta = 0.3

# prepare grids
flags = s2.create(FlagGrid)
vel = s2.create(MACGrid)
velT = s2.create(MACGrid)
density = s2.create(RealGrid)
pressure = s2.create(RealGrid)
W = s2.create(RealGrid)

gsLoad = vec3(res1,int(2.0*res1),res1)
sLoader = Solver(name='main', gridSize = gsLoad, dim=3) 
velIn = sLoader.create(MACGrid)

# noise field
noise = s2.create(NoiseField, loadFromFile=True)
noise.posScale = vec3(0)
noise.clamp = True
noise.clampNeg = 0
noise.clampPos = 1
noise.valScale = 1
noise.valOffset = 0.75
noise.timeAnim = 0.2

bWidth=0
flags.initDomain(boundaryWidth=bWidth) 
flags.fillGrid()
setOpenBound(flags,bWidth,'yY',FlagOutflow|FlagEmpty) 

if (GUI):
	gui = Gui()
	gui.show()
	gui.nextVec3Display() # hide velocity display
	gui.nextVec3Display()
	gui.nextVec3Display()
	#gui.pause()

source = s2.create(Cylinder, center=gs2*vec3(0.5,0.05,0.5), radius=res2*0.1, z=gs2*vec3(0, 0.02, 0))
W.multConst(0)
W.addConst(wScalar)
	
#main loop
for t in range(numFrames):
	densityInflow(flags=flags, density=density, noise=noise, shape=source, scale=1, sigma=0.5)
		
	advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2)    
	advectSemiLagrange(flags=flags, vel=vel, grid=vel,     order=2, openBounds=True, boundaryWidth=bWidth+1)
	resetOutflow(flags=flags,real=density) 
	
	setWallBcs(flags=flags, vel=vel)
	addBuoyancy(density=density, vel=vel, gravity=vec3(0,-1e-3*factor,0), flags=flags)
	
	velIn.load( input_uni % (t) )
	interpolateMACGrid( source=velIn, target=velT )
	velT.multConst(vec3(factor))

	PD_fluid_guiding(vel=vel, velT=velT, flags=flags, weight=W, blurRadius=beta, pressure=pressure, \
		tau = tau, sigma = sigma, theta = theta, preconditioner = PcMGStatic, zeroPressureFixing=True )
	
	setWallBcs(flags=flags, vel=vel)
	if 0:
		projectPpmFull( density, output_ppm % (t) , 0, 2.0 );
	density.save(output_uni % (t))
	
	s2.step()

