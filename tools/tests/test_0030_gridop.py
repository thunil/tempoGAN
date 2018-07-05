#
# Basic test, grid operators
# 

import sys
print ("Running python "+sys.version)

from manta import *
from helperInclude import *


# solver params
gs  = vec3(10, 20, 30)
s   = Solver(name='main', gridSize = gs, dim=3)

# prepare grids
rlg1  = s.create(RealGrid)
rlg2  = s.create(RealGrid)
rlg3  = s.create(RealGrid)
vcg1  = s.create(MACGrid)
vcg2  = s.create(MACGrid)
vcg3  = s.create(MACGrid)
int1  = s.create(IntGrid)
int2  = s.create(IntGrid)
int3  = s.create(IntGrid)
vcgTmp= s.create(MACGrid)

genRefFiles = getGenRefFileSetting() 
if (genRefFiles==1):
	# manually init result
	rlg1.setConst( 1.1 )
	rlg2.setConst( 1.2 )
	rlg3.setConst( 2.9 )

	#vcg1.setConst( vec3(1.2, 1.2, 1.2) )
	#vcg2.setConst( vec3(0.5, 0.5, 0.5) )
	#vcg3.setConst( vec3(1.95, 1.95, 1.95) )

	vcg1.setConst( vec3(1.25, 1.25, 1.25) )
	vcg2.setConst( vec3(0.5, 0.5, 0.5) )
	vcg3.setConst( vec3(1.95, 1.95, 1.95) )

	int1.setConst( 125 )
	int2.setConst( 6 )
	int3.setConst( 143 )
else:	
# real test run, perform basic calculations

	rlg1.setConst( 1.0 )
	rlg2.setConst( 2.4 )
	rlg3.setConst( 9.6 )
	rlg1.addConst (0.1) # 1.1
	rlg2.multConst(0.5)  # 1.2
	rlg3.copyFrom( rlg1 )  # 1.1
	rlg3.add(rlg2)  # 2.3
	rlg3.addScaled(rlg2, 0.5) # 2.9
	#print "r1 %f , r2 %f , r3 %f " % ( rlg1.getMaxAbs() , rlg2.getMaxAbs() , rlg3.getMaxAbs() )

	vcg1.setConst( vec3(1.0, 1.0, 1.0) )
	vcg2.setConst( vec3(1.0, 1.0, 1.0) )
	vcg3.setConst( vec3(9.0, 9.0, 9.0) )

	vcg1.addConst ( vec3(0.25,0.25,0.25) ) # 1.25

	vcg2.multConst( vec3(0.5,0.5,0.5) ) # 0.5

	vcgTmp.setConst( vec3(1.2, 1.2, 1.2) )
	vcg3.copyFrom( vcgTmp )  # 1.2
	vcg3.add(vcg2) # 1.7
	vcg3.addScaled(vcg2, vec3(0.5, 0.5, 0.5) ) # 1.95
	#print "v1 %s , v2 %s , v3 %s " % ( vcg1.getMaxAbs() , vcg2.getMaxAbs(), vcg3.getMaxAbs() )

	int1.setConst( 123 )
	int2.setConst( 2 )
	int3.setConst( 9 )
	int1.addConst ( 2 ) # 125
	int2.multConst( 3 ) # 6
	int3.copyFrom( int1 ) # 125
	int3.add(int2)  # 131
	int3.addScaled(int2, 2) # 143
	#print "i1 %s , i2 %s , i3 %s " % ( int1.getMaxAbs() , int2.getMaxAbs() , int3.getMaxAbs() )


# verify

doTestGrid( sys.argv[0], "rlg1", s, rlg1 , threshold=1e-07 , thresholdStrict=1e-14  )
doTestGrid( sys.argv[0], "rlg2", s, rlg2 , threshold=1e-07 , thresholdStrict=1e-14  )
doTestGrid( sys.argv[0], "rlg3", s, rlg3 , threshold=1e-07 , thresholdStrict=1e-14  )

doTestGrid( sys.argv[0], "vcg1", s, vcg1 , threshold=5e-07 , thresholdStrict=1e-14  )
doTestGrid( sys.argv[0], "vcg2", s, vcg2 , threshold=5e-07 , thresholdStrict=1e-14  )
doTestGrid( sys.argv[0], "vcg3", s, vcg3 , threshold=5e-07 , thresholdStrict=1e-14  )

doTestGrid( sys.argv[0], "int1", s, int1 , threshold=1e-14 , thresholdStrict=1e-14  )
doTestGrid( sys.argv[0], "int2", s, int2 , threshold=1e-14 , thresholdStrict=1e-14  )
doTestGrid( sys.argv[0], "int3", s, int3 , threshold=1e-14 , thresholdStrict=1e-14  )

