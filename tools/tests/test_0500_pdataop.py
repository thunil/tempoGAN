#
# Basic test, grid operators
# 

import sys
print ("Running python "+sys.version)

from manta import *
from helperInclude import *


# solver params
gs  = vec3(12, 19, 31)
s   = Solver(name='main', gridSize = gs, dim=3)
pp  = s.create(BasicParticleSystem) 

# prepare grids
rlg1  = pp.create(PdataReal)
rlg2  = pp.create(PdataReal)
rlg3  = pp.create(PdataReal)
vcg1  = pp.create(PdataVec3)
vcg2  = pp.create(PdataVec3)
vcg3  = pp.create(PdataVec3)
int1  = pp.create(PdataInt)
int2  = pp.create(PdataInt)
int3  = pp.create(PdataInt)

genRefFiles = getGenRefFileSetting() 
if (genRefFiles==1):
	# generate and save particles
	addTestParts( pp , 10 )
	#pp.printParts(start=1, stop=10, printIndex=True)
	pp.save( referenceFilename( sys.argv[0], "parts" ) )

	# manually init result
	rlg1.setConst( 1.1 )
	rlg2.setConst( 1.2 )
	rlg3.setConst( 2.9 )
           
	vcg1.setConst( vec3(1.2, 1.2, 1.2) )
	vcg2.setConst( vec3(0.5, 0.5, 0.5) )
	vcg3.setConst( vec3(1.95, 1.95, 1.95) )
           
	int1.setConst( 125 )
	int2.setConst( 6 )
	int3.setConst( 143 )
else:	
	# real test run, perform basic calculations 
	pp.load( referenceFilename( sys.argv[0], "parts" ) )

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
	vcg1.addConst ( vec3(0.2,0.2,0.2) ) # 1.2
	vcg2.multConst( vec3(0.5,0.5,0.5) ) # 0.5
	vcg3.copyFrom( vcg1 )  # 1.2
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
	#int1.printPdata(start=2, stop=4, printIndex=True)
	#int2.printPdata(start=2, stop=4, printIndex=True)
	#int3.printPdata(start=2, stop=4, printIndex=True)
	#print "i1 %s , i2 %s , i3 %s " % ( int1.getMaxAbs() , int2.getMaxAbs() , int3.getMaxAbs() ) 

# verify

# note the strict/double prec threshold is ridiculously un-strict - problem is the float rounding in the uni/raw files, which makes this test pretty meaningless for doubles...
doTestGrid( sys.argv[0], "rlg1", pp, rlg1 , threshold=1e-07 , thresholdStrict=1e-06  )
doTestGrid( sys.argv[0], "rlg2", pp, rlg2 , threshold=1e-07 , thresholdStrict=1e-06  )
doTestGrid( sys.argv[0], "rlg3", pp, rlg3 , threshold=1e-07 , thresholdStrict=1e-06  )

doTestGrid( sys.argv[0], "vcg1", pp, vcg1 , threshold=5e-07 , thresholdStrict=1e-06  )
doTestGrid( sys.argv[0], "vcg2", pp, vcg2 , threshold=5e-07 , thresholdStrict=1e-06  )
doTestGrid( sys.argv[0], "vcg3", pp, vcg3 , threshold=5e-07 , thresholdStrict=1e-06  )

doTestGrid( sys.argv[0], "int1", pp, int1 , threshold=1e-14 , thresholdStrict=1e-14  )
doTestGrid( sys.argv[0], "int2", pp, int2 , threshold=1e-14 , thresholdStrict=1e-14  )
doTestGrid( sys.argv[0], "int3", pp, int3 , threshold=1e-14 , thresholdStrict=1e-14  )

