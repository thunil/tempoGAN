#
# 2D interpolation functions
# 

from manta import *
import sys
from helperInclude import *

doGui = False;

# dimension two/three d
dim = 2
# how much to upres the XL sim?
upres = 5
# resolution
res = 60

# small
smgs = vec3(res/upres,int(1.0*res/upres),1)
if (dim==2): smgs.z = 1  # 2D
smSolv = Solver(name='smallSolver', gridSize = smgs, dim=dim)
sm_flags    = smSolv.create(FlagGrid)
sm_flags.initDomain()
sm_flags.fillGrid()

# med
gs = vec3(res,int(1.0*res),res)
if (dim==2): gs.z = 1  # 2D
normSolv = Solver(name='mainSolver', gridSize = gs, dim=dim)
flags    = normSolv.create(FlagGrid)
flags.initDomain()
flags.fillGrid()

# larger solver, recompute sizes...
xl_gs = vec3(upres*gs.x,upres*gs.y,upres*gs.z)
if (dim==2): xl_gs.z = 1  # 2D
xlSolv = Solver(name='largerSolver', gridSize = xl_gs, dim=dim)
xl_flags   = xlSolv.create(FlagGrid)
xl_flags.initDomain()
xl_flags.fillGrid()

# data

sm_density  = smSolv.create(RealGrid)
density     = normSolv.create(RealGrid)
xl_density  = xlSolv.create(RealGrid)
density2    = normSolv.create(RealGrid)
sm_density2 = smSolv.create(RealGrid)

sm_macvel   = smSolv.create(MACGrid)
macvel      = normSolv.create(MACGrid)
xl_macvel   = xlSolv.create(MACGrid)
macvel2     = normSolv.create(MACGrid)
sm_macvel2  = smSolv.create(MACGrid)

sm_v3    = smSolv.create(VecGrid)
v3       = normSolv.create(VecGrid)
xl_v3    = xlSolv.create(VecGrid)
v32      = normSolv.create(VecGrid)
sm_v32   = smSolv.create(VecGrid)

# sources

xsourceVel  = normSolv.create(Cylinder, center=xl_gs*vec3(0.5,0.5,0.5), radius=xl_gs.x*0.151, z=xl_gs*vec3(0.151, 0, 0))
smsourceVel = normSolv.create(Cylinder, center=smgs*vec3(0.5,0.5,0.5), radius=smgs.x*0.251, z=smgs*vec3(0.151, 0, 0))

if doGui and (GUI):
	gui = Gui()
	gui.show()

# linear interpol

ords = 1;

for t in range(1):

	sm_density2.clear();
	sm_density.clear();
	density2.clear();
	density.clear();
	xl_density.clear();

	# small to large
	smsourceVel.applyToGrid( grid=sm_density , value=1)
	smsourceVel.applyToGrid( grid=sm_v3      , value=vec3(1) )
	smsourceVel.applyToGrid( grid=sm_macvel     , value=vec3(1) )
	
	interpolateGrid( target=density, source=sm_density , orderSpace=ords )
	interpolateGrid( target=xl_density, source=density , orderSpace=ords )
	interpolateGrid( target=density2, source=xl_density , orderSpace=ords )
	interpolateGrid( target=sm_density2, source=density2 , orderSpace=ords )
	#sm_density2.sub(sm_density)
	
	interpolateGridVec3( target=v3, source=sm_v3 , orderSpace=ords )
	interpolateGridVec3( target=xl_v3, source=v3 , orderSpace=ords )
	interpolateGridVec3( target=v32, source=xl_v3 , orderSpace=ords )
	interpolateGridVec3( target=sm_v32, source=v32 , orderSpace=ords )
	
	interpolateMACGrid( target=macvel, source=sm_macvel , orderSpace=ords )
	interpolateMACGrid( target=xl_macvel, source=macvel , orderSpace=ords )
	interpolateMACGrid( target=macvel2, source=xl_macvel , orderSpace=ords )
	interpolateMACGrid( target=sm_macvel2, source=macvel2 , orderSpace=ords )

	smSolv.step()
	normSolv.step()
	xlSolv.step()    
	if doGui:
		gui.pause()

	doTestGrid( sys.argv[0], "scalar1"  , normSolv, density   , threshold=1e-06 , thresholdStrict=1e-14  )
	doTestGrid( sys.argv[0], "scalar2"  , smSolv, sm_density  , threshold=1e-06 , thresholdStrict=1e-14  )
	doTestGrid( sys.argv[0], "scalar3"  , smSolv, sm_density2 , threshold=1e-06 , thresholdStrict=1e-14  )
	doTestGrid( sys.argv[0], "vec3t1"   , normSolv, v3        , threshold=1e-06 , thresholdStrict=1e-14  )
	doTestGrid( sys.argv[0], "vec3t2"   , smSolv  , sm_v3     , threshold=1e-06 , thresholdStrict=1e-14  )
	doTestGrid( sys.argv[0], "vec3t3"   , smSolv  , sm_v32    , threshold=1e-06 , thresholdStrict=1e-14  )
	doTestGrid( sys.argv[0], "macvel1"  , normSolv, macvel    , threshold=1e-06 , thresholdStrict=1e-14  )
	doTestGrid( sys.argv[0], "macvel2"  , smSolv  , sm_macvel , threshold=1e-06 , thresholdStrict=1e-14  )
	doTestGrid( sys.argv[0], "macvel3"  , smSolv  , sm_macvel2, threshold=1e-06 , thresholdStrict=1e-14  )



# high order interpol

ords = 2;

for t in range(1):

	sm_density2.clear();
	sm_density.clear();
	density2.clear();
	density.clear();
	xl_density.clear();
	sm_v3.clear();
	sm_macvel.clear();

	# small to large
	smsourceVel.applyToGrid( grid=sm_density , value=1)
	smsourceVel.applyToGrid( grid=sm_v3      , value=vec3(1) )
	smsourceVel.applyToGrid( grid=sm_macvel     , value=vec3(1) )
	
	interpolateGrid( target=density, source=sm_density , orderSpace=ords )
	interpolateGrid( target=xl_density, source=density , orderSpace=ords )
	interpolateGrid( target=density2, source=xl_density , orderSpace=ords )
	interpolateGrid( target=sm_density2, source=density2 , orderSpace=ords )
	#sm_density2.sub(sm_density)
	
	interpolateGridVec3( target=v3, source=sm_v3 , orderSpace=ords )
	interpolateGridVec3( target=xl_v3, source=v3 , orderSpace=ords )
	interpolateGridVec3( target=v32, source=xl_v3 , orderSpace=ords )
	interpolateGridVec3( target=sm_v32, source=v32 , orderSpace=ords )
	
	interpolateMACGrid( target=macvel, source=sm_macvel , orderSpace=ords )
	interpolateMACGrid( target=xl_macvel, source=macvel , orderSpace=ords )
	interpolateMACGrid( target=macvel2, source=xl_macvel , orderSpace=ords )
	interpolateMACGrid( target=sm_macvel2, source=macvel2 , orderSpace=ords )

	smSolv.step()
	normSolv.step()
	xlSolv.step()    
	if doGui:
		gui.pause()

	doTestGrid( sys.argv[0], "hi_scalar1"  , normSolv, density   , threshold=1e-05 , thresholdStrict=1e-14  )
	doTestGrid( sys.argv[0], "hi_scalar2"  , smSolv, sm_density  , threshold=1e-05 , thresholdStrict=1e-14  )
	doTestGrid( sys.argv[0], "hi_scalar3"  , smSolv, sm_density2 , threshold=1e-05 , thresholdStrict=1e-14  )
	doTestGrid( sys.argv[0], "hi_vec3t1"   , normSolv, v3        , threshold=1e-05 , thresholdStrict=1e-14  )
	doTestGrid( sys.argv[0], "hi_vec3t2"   , smSolv  , sm_v3     , threshold=1e-05 , thresholdStrict=1e-14  )
	doTestGrid( sys.argv[0], "hi_vec3t3"   , smSolv  , sm_v32    , threshold=1e-05 , thresholdStrict=1e-14  )
	doTestGrid( sys.argv[0], "hi_macvel1"  , normSolv, macvel    , threshold=1e-05 , thresholdStrict=1e-14  )
	doTestGrid( sys.argv[0], "hi_macvel2"  , smSolv  , sm_macvel , threshold=1e-05 , thresholdStrict=1e-14  )
	doTestGrid( sys.argv[0], "hi_macvel3"  , smSolv  , sm_macvel2, threshold=1e-05 , thresholdStrict=1e-14  )


