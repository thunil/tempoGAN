#
# Testing shapes
# 

import sys
from manta import *
from helperInclude import *

# solver params
res = 42
gs  = vec3(res,res,res)
s   = Solver(name='main', gridSize = gs, dim=3)


# prepare grids

density = s.create(RealGrid)
vel     = s.create(MACGrid)

density.setConst( 0. )
vel.setConst( vec3(0, 0, 0) )

sph1 = s.create(Sphere, center=gs*vec3(0.3,0.4,0.5), radius=res*0.2)
sph1.applyToGrid(grid=density, value=0.302)

velVal = vec3(0.1,0.1,0.4)
sph2 = s.create(Sphere, center=gs*vec3(0.6,0.5,0.4), radius=res*0.25)
sph2.applyToGrid(grid=vel, value=velVal )

doTestGrid( sys.argv[0], "densSph" , s, density  , threshold=1e-07 , thresholdStrict=1e-10   )
doTestGrid( sys.argv[0], "velSph"  , s, vel      , threshold=1e-07 , thresholdStrict=1e-10   )


# sphere

density.setConst( 0. )
vel.setConst( vec3(0, 0, 0) )

box1 = s.create(Box, p0=gs*vec3(0.2,0.2,0.3), p1=gs*vec3(0.9,0.8,0.9) )
box1.applyToGrid(grid=density, value=0.812)

velVal = vec3(0.5,0.1,0.1)
box2 = s.create(Box, p0=gs*vec3(0.2,0.2,0.3), p1=gs*vec3(0.9,0.8,0.9) )
box2.applyToGrid(grid=vel, value=velVal )

doTestGrid( sys.argv[0], "densBox" , s, density  , threshold=1e-07 , thresholdStrict=1e-10   )
doTestGrid( sys.argv[0], "velBox"  , s, vel      , threshold=1e-07 , thresholdStrict=1e-10   )


# cylinder

density.setConst( 0. )
vel.setConst( vec3(0, 0, 0) )

cyl1 = s.create(Cylinder, center=gs*vec3(0.5,0.5,0.5), radius=res*0.2, z=gs*vec3(0, 0.3, 0))
cyl1.applyToGrid(grid=density, value=0.432)

velVal = vec3(0.4,0.3,0.2)
cyl2 = s.create(Cylinder, center=gs*vec3(0.5,0.5,0.5), radius=res*0.2, z=gs*vec3(0, 0.3, 0))
cyl2.applyToGrid(grid=vel, value=velVal )

doTestGrid( sys.argv[0], "densCyl" , s, density  , threshold=1e-07 , thresholdStrict=1e-10   )
doTestGrid( sys.argv[0], "velCyl"  , s, vel      , threshold=1e-07 , thresholdStrict=1e-10   )


