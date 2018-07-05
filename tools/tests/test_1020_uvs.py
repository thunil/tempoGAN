#
# Simple uv advection & init setup
# 
from manta import *
import os, shutil, math, sys
from helperInclude import *

# dimension two/three d
dim = 2

# how many grids of uv coordinates to use (more than 2 usually dont pay off here)
uvs = 3

# simulation resolution
res = 50
gs = vec3(res,int(1.5*res),res)
if (dim==2): gs.z = 1  # 2D

# setup low-res sim
sm = Solver(name='main', gridSize = gs, dim=dim)
sm.timestep = 0.5

# helper objects
source    = sm.create(Cylinder, center=gs*vec3(0.3,0.4,0.5), radius=res*0.10, z=gs*vec3(0.10, 0, 0))
sourceVel = sm.create(Cylinder, center=gs*vec3(0.3,0.4,0.5), radius=res*0.151, z=gs*vec3(0.151, 0, 0))


# init lower res solver & grids
flags = sm.create(FlagGrid)
flags.initDomain()
flags.fillGrid()

# create the array of uv grids
uv = []
for i in range(uvs):
	uv.append(i)
	uv[i] = sm.create(VecGrid) # empty name, bad... 
	resetUvGrid( uv[i] )

vel       = sm.create(MACGrid) 
density   = sm.create(RealGrid)
pressure  = sm.create(RealGrid)

if (0 and GUI):
	gui = Gui()
	gui.show(); gui.pause()

# compute one vel field
source.applyToGrid( grid=density , value=1. )
sourceVel.applyToGrid( grid=vel , value=vec3(5,0,0) )
setWallBcs(flags=flags, vel=vel)	
addBuoyancy(density=density, vel=vel, gravity=vec3(0,-1e-2,0), flags=flags)
solvePressure(flags=flags, vel=vel, pressure=pressure, cgMaxIterFac=2.0, cgAccuracy=1e-06 )
setWallBcs(flags=flags, vel=vel)

# main loop
for t in range(20):

	for i in range(uvs):
		advectSemiLagrange(flags=flags, vel=vel, grid=uv[i], order=1) 
		updateUvWeight( resetTime=11.0 , index=i, numUvs=uvs, uv=uv[i] ); 
	
	#gui.screenshot( 'out_%04d.png' % t );
	sm.step()

doTestGrid( sys.argv[0],"uv0" , sm, uv[0] , threshold=0.006 , thresholdStrict=1e-10 )
doTestGrid( sys.argv[0],"uv1" , sm, uv[1] , threshold=0.006 , thresholdStrict=1e-10 )
doTestGrid( sys.argv[0],"uv2" , sm, uv[2] , threshold=0.006 , thresholdStrict=1e-10 )

