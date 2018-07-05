#
# 3D interpolation functions
# 

from manta import *
import sys
from helperInclude import *

showGui = False

# dimension two/three d
dim = 3
# how much to upres the XL sim?
upres = 2
# resolution
res = 40

# small
smgs = vec3(res/upres,res/upres,res/upres)
if (dim==2): smgs.z = 1  # 2D
smSolv = Solver(name='smallSolver', gridSize = smgs, dim=dim, fourthDim=smgs.x)
sm_flags    = smSolv.create(FlagGrid)
sm_flags.initDomain()
sm_flags.fillGrid()

# med
gs = vec3(res,res,res)
if (dim==2): gs.z = 1  # 2D
normSolv = Solver(name='mainSolver', gridSize = gs, dim=dim, fourthDim=gs.x)
flags    = normSolv.create(FlagGrid)
flags.initDomain()
flags.fillGrid()

# larger solver, recompute sizes...
xlgs = vec3(upres*gs.x,upres*gs.y,upres*gs.z)
if (dim==2): xlgs.z = 1  # 2D
xlSolv = Solver(name='largerSolver', gridSize = xlgs, dim=dim, fourthDim=xlgs.x)
xl_flags   = xlSolv.create(FlagGrid)
xl_flags.initDomain()
xl_flags.fillGrid()

# data

sm_density  = smSolv.create(  Grid4Real)
density     = normSolv.create(Grid4Real)
xl_density  = xlSolv.create(  Grid4Real)
density2    = normSolv.create(Grid4Real)
sm_density2 = smSolv.create(  Grid4Real)

sm_v3    = smSolv.create( Grid4Vec4)
v3       = normSolv.create(Grid4Vec4)
xl_v3    = xlSolv.create(Grid4Vec4)
v32      = normSolv.create(Grid4Vec4)
sm_v32   = smSolv.create(Grid4Vec4)

# for display only
sm_densDisp  = smSolv.create(  RealGrid)
densDisp     = normSolv.create(RealGrid)
xl_densDisp  = xlSolv.create(  RealGrid)
densDisp2    = normSolv.create(RealGrid)
sm_densDisp2 = smSolv.create(  RealGrid)

sm_velDisp    = smSolv.create( VecGrid)
velDisp       = normSolv.create(VecGrid)
xl_velDisp    = xlSolv.create(VecGrid)
velDisp2      = normSolv.create(VecGrid)
sm_velDisp2   = smSolv.create(VecGrid)

# sources

rs = smgs.x*0.3;
re = smgs.x*0.7;
rstart = Vec4(rs,rs,rs,rs);
rend   = Vec4(re,re,re,re);
# we have to use the setRegion function here, as regular shapes dont work in 4d so far

if showGui and (GUI):
	gui = Gui()
	gui.show()

# linear interpol only in 4d

for t in range(1):

	sm_density2.clear();
	sm_density.clear();
	density2.clear();
	density.clear();
	xl_density.clear();

	# small to large
	setRegion4d( sm_density, start=rstart, end=rend , value=1 );
	setRegion4dVec4( sm_v3 , start=rstart, end=rend , value=Vec4(1,1,1,1) );
	
	interpolateGrid4d( target=density, source=sm_density  )
	interpolateGrid4d( target=xl_density, source=density  )
	interpolateGrid4d( target=density2, source=xl_density )
	interpolateGrid4d( target=sm_density2, source=density2)
	#sm_density2.sub(sm_density)
	
	interpolateGrid4dVec( target=v3, source=sm_v3 )
	interpolateGrid4dVec( target=xl_v3, source=v3 )
	interpolateGrid4dVec( target=v32, source=xl_v3 )
	interpolateGrid4dVec( target=sm_v32, source=v32 )

	#for t in range(
	t=0.5;
	getSliceFrom4d(src=density, dst=densDisp, srct=int(gs.x*t) ); 
	getSliceFrom4d(src=density2, dst=densDisp2, srct=int(gs.x*t) ); 
	getSliceFrom4d(src=sm_density, dst=sm_densDisp, srct=int(smgs.x*t) ); 
	getSliceFrom4d(src=sm_density2, dst=sm_densDisp2, srct=int(smgs.x*t) ); 
	getSliceFrom4d(src=xl_density, dst=xl_densDisp, srct=int(xlgs.x*t) ); 

	getSliceFrom4dVec(src=v3, dst=velDisp, srct=int(gs.x*t) ); 
	getSliceFrom4dVec(src=v32, dst=velDisp2, srct=int(gs.x*t) ); 
	getSliceFrom4dVec(src=sm_v3, dst=sm_velDisp, srct=int(smgs.x*t) ); 
	getSliceFrom4dVec(src=sm_v32, dst=sm_velDisp2, srct=int(smgs.x*t) ); 
	getSliceFrom4dVec(src=xl_v3, dst=xl_velDisp, srct=int(xlgs.x*t) ); 

	smSolv.step()
	normSolv.step()
	xlSolv.step()    
	if showGui:
		gui.pause()

	for i in range(10):
		t=i*0.1;
		getSliceFrom4d(src=density, dst=densDisp, srct=int(gs.x*t) ); 
		getSliceFrom4d(src=density2, dst=densDisp2, srct=int(gs.x*t) ); 
		getSliceFrom4d(src=sm_density, dst=sm_densDisp, srct=int(smgs.x*t) ); 
		getSliceFrom4d(src=sm_density2, dst=sm_densDisp2, srct=int(smgs.x*t) ); 
		getSliceFrom4d(src=xl_density, dst=xl_densDisp, srct=int(xlgs.x*t) ); 
		smSolv.step()

	doTestGrid( sys.argv[0], "scalar1"  , normSolv, density   , threshold=1e-05 , thresholdStrict=1e-14  )
	doTestGrid( sys.argv[0], "scalar2"  , smSolv, sm_density  , threshold=1e-05 , thresholdStrict=1e-14  )
	doTestGrid( sys.argv[0], "scalar3"  , smSolv, sm_density2 , threshold=1e-05 , thresholdStrict=1e-14  )
	doTestGrid( sys.argv[0], "vec3t1"   , normSolv, v3        , threshold=1e-05 , thresholdStrict=1e-14  )
	doTestGrid( sys.argv[0], "vec3t2"   , smSolv  , sm_v3     , threshold=1e-05 , thresholdStrict=1e-14  )
	doTestGrid( sys.argv[0], "vec3t3"   , smSolv  , sm_v32    , threshold=1e-05 , thresholdStrict=1e-14  )



