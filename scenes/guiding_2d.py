#
# Simple example scene for a 2D circularly guided simulation
# Main params:
# 	W = guiding weight, different strengths for top / bottom halves
# 	beta = blur radius
#	scale = up-res factor
#
from manta import *
import math

# solver params
res0 = 64
scale = 2
res = res0*scale
gs = vec3(res,res,1)
s = Solver(name='main', gridSize = gs, dim=2)
s.timestep = 2.0/scale
timings = Timings()

# IOP (=1), ADMM (=2) or PD (=3)

# params
valAtMin = 1
valAtMax = 5
beta = 2
output_ppm = 'guiding2D_%04d.ppm'
output_png = 'guiding2D_%04d.png'
tau = 1.0
sigma = 0.99/tau
theta = 1.0

# prepare grids
flags = s.create(FlagGrid)
vel = s.create(MACGrid)
velT = s.create(MACGrid)
density = s.create(RealGrid)
pressure = s.create(RealGrid)
W = s.create(RealGrid)

bWidth=1
flags.initDomain(boundaryWidth=bWidth) 
flags.fillGrid()

if (GUI):
	gui = Gui()
	gui.show()
	gui.nextVec3Display()
	gui.nextVec3Display()
	gui.nextVec3Display()
	#gui.pause()

source = s.create(Cylinder, center=gs*vec3(0.5,0.2,0.5), radius=gs.y*0.14, z=gs*vec3(0, 0.02*1.5, 0))
getSpiralVelocity2D(flags=flags, vel=velT, strength=0.5*scale)

setGradientYWeight(W=W, minY=0,     maxY=res/2, valAtMin=valAtMin, valAtMax=valAtMin)
setGradientYWeight(W=W, minY=res/2, maxY=res,   valAtMin=valAtMax, valAtMax=valAtMax)

#main loop
for t in range(100*scale):
	resetOutflow(flags=flags,real=density)

	if t<100*scale:
		source.applyToGrid(grid=density, value=1)
	
	advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2)    
	advectSemiLagrange(flags=flags, vel=vel, grid=vel,     order=2)
	
	setWallBcs(flags=flags, vel=vel)
	addBuoyancy(density=density, vel=vel, gravity=vec3(0,0.25*scale*-4e-3,0), flags=flags)

	PD_fluid_guiding(vel=vel, velT=velT, flags=flags, weight=W, blurRadius=beta, pressure=pressure, \
		tau = tau, sigma = sigma, theta = theta, preconditioner = PcMGStatic, zeroPressureFixing=True ) 

	setWallBcs(flags=flags, vel=vel)
	if 0 and (t%scale==0):
		projectPpmFull( density, output_ppm % (t/scale) , 0, 1.0 );
		#gui.screenshot(output_png % (t/scale) )
	
	#timings.display()
	s.step()

