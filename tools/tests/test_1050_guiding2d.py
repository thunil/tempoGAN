#
# Test in line with 2d guiding example
#

import math
import sys
from manta import *
from helperInclude import *

# solver params
res0 = 30
scale = 2
res = res0*scale
gs = vec3(res,res,1)
s = Solver(name='main', gridSize = gs, dim=2)
s.timestep = 2.0/scale

# IOP (=1), ADMM (=2) or PD (=3)

# params
valAtMin = 1
valAtMax = 5
beta = 2
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

if 0 and (GUI):
	gui = Gui()
	gui.show()

source = s.create(Cylinder, center=gs*vec3(0.5,0.3,0.5), radius=gs.y*0.14, z=gs*vec3(0, 0.04*1.5, 0))
getSpiralVelocity2D(flags=flags, vel=velT, strength=1.5*scale)
setGradientYWeight(W=W, minY=0,     maxY=res/2, valAtMin=valAtMin, valAtMax=valAtMin)
setGradientYWeight(W=W, minY=res/2, maxY=res,   valAtMin=valAtMax, valAtMax=valAtMax)

#main loop
for t in range(5):
	resetOutflow(flags=flags,real=density)

	source.applyToGrid(grid=density, value=1)
	
	advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2, clampMode=1)    
	advectSemiLagrange(flags=flags, vel=vel, grid=vel,     order=2, clampMode=1)
	
	setWallBcs(flags=flags, vel=vel)
	addBuoyancy(density=density, vel=vel, gravity=vec3(0,0.25*scale*-1e-2,0), flags=flags)

	PD_fluid_guiding(vel=vel, velT=velT, flags=flags, weight=W, blurRadius=beta, pressure=pressure, \
		tau = tau, sigma = sigma, theta = theta, preconditioner = 1 )

	setWallBcs(flags=flags, vel=vel)
	s.step()

# this scene runs through with doubles, but not supported for tests so far
if getFloatSetting()==2:
	print("Note: test disabled for now for double precision!")
	exit(1)

# check final state
doTestGrid( sys.argv[0],"dens" , s, density , threshold=0.0001 , thresholdStrict=1e-10 )
doTestGrid( sys.argv[0],"vel"  , s, vel     , threshold=0.0001 , thresholdStrict=1e-10 )

