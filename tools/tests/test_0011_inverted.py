#
# Inverted test, this should "fail", ie give a large difference...
# 

import sys
from manta import *
from helperInclude import *

# solver params
gs  = vec3(17, 177, 27)
s   = Solver(name='main', gridSize = gs, dim=3)


# prepare grids
density = s.create(RealGrid)
dummy   = s.create(RealGrid)

# set some value != 0
density.setConst( 25.01 )
dummy.setConst( -25.00 )

# verify , note - this should fail!
if (getGenRefFileSetting()==1):
	doTestGrid( sys.argv[0], "dens" , s, density )
else:
	doTestGrid( sys.argv[0], "dens" , s, dummy , threshold=50., thresholdStrict=50. , invertResult=True )


