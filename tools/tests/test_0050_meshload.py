#
# Simple mesh load & sdf test
# 
import sys
print ("Running python "+sys.version)

from manta import *
from helperInclude import *

meshfile = 'test_0050_meshload.obj'

# resolution for level set / output mesh
res = 100 
gs = vec3(res,res,res)
s = Solver(name='main', gridSize = vec3(res,res,res) , dim=3)

# kernel radius for surface creation
radiusFactor = 2.5

# triangle scale relative to cell size
scale = 0.5

# prepare grids and particles
flags    = s.create(FlagGrid)
phi      = s.create(LevelsetGrid)
mesh     = s.create(Mesh)

# scene setup
flags.initDomain(boundaryWidth=0)
	
if 0 and (GUI):
	gui = Gui()
	gui.show()
	gui.pause()

#main 
mesh.load( meshfile )
mesh.scale( vec3(res/3.0) );
mesh.offset( gs*0.5 );
mesh.computeLevelset(phi, 2., -1.);

s.step()
	
doTestGrid( sys.argv[0], "phi" , s, phi  , threshold=1e-05 , thresholdStrict=5e-08 ) # higher double threshold for windows...

