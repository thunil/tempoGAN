#
# k-epsilon test case
# 
#
import sys
from manta import *
from helperInclude import *

# solver params
res    = 70
frames = 32

if getVisualSetting():
	# in visual mode
	res    = 102 * getVisualSetting()
	frames = 100 * getVisualSetting()

gs = vec3(res,res/2,res/2)
s = Solver(name='main', gridSize = gs)
s.timestep = 1.2

if getVisualSetting():
	s.timestep = 1.7 / getVisualSetting()

velInflow = vec3(0.52,0,0)

# prepare grids
flags = s.create(FlagGrid)
pressure = s.create(RealGrid, show=False)
vel = s.create(MACGrid)

k = s.create(RealGrid)
eps = s.create(RealGrid)
prod = s.create(RealGrid)
nuT= s.create(RealGrid)
strain= s.create(RealGrid)
vc=s.create(MACGrid)
temp=s.create(RealGrid)

# noise field
noise = s.create(NoiseField, loadFromFile=True )
noise.timeAnim = 0

# turbulence particles
turb = s.create(TurbulenceParticleSystem, noise=noise)

flags.initDomain()
flags.fillGrid()

# obstacle grid
for i in range(4):
	for j in range(4):
		obs = s.create(Sphere, center=gs*vec3(0.2,(i+1)/5.0,(j+1)/5.0), radius=res*0.025)
		obs.applyToGrid(grid=flags,value=FlagObstacle)

sdfgrad = obstacleGradient(flags)
sdf = obstacleLevelset(flags)
bgr = s.create(Mesh)
sdf.createMesh(bgr)

# particle inflow
box = s.create(Box, center = gs*vec3(0.05,0.43,0.6), size=gs*vec3(0.02,0.005,0.07))

# turbulence parameters
L0 = 0.01
mult = 0.1
intensity = 0.1
nu = 0.1
prodMult = 2.5
enableDiffuse = True

if 0 and (GUI):
	gui = Gui()
	gui.setBackgroundMesh(bgr)
	gui.show()

KEpsilonBcs(flags=flags,k=k,eps=eps,intensity=intensity,nu=nu,fillArea=True)

#main loop
for t in range(frames):
	turb.seed(box,500)
	turb.advectInGrid(flags=flags, vel=vel, integrationMode=IntRK4)
	turb.synthesize(flags=flags, octaves=1, k=k, switchLength=5, L0=L0, scale=mult, inflowBias=velInflow)
	#turb.projectOutside(sdfgrad)
	turb.deleteInObstacle(flags)

	KEpsilonBcs(flags=flags,k=k,eps=eps,intensity=intensity,nu=nu,fillArea=False)
	advectSemiLagrange(flags=flags, vel=vel, grid=k, order=1)
	advectSemiLagrange(flags=flags, vel=vel, grid=eps, order=1)
	KEpsilonBcs(flags=flags,k=k,eps=eps,intensity=intensity,nu=nu,fillArea=False)
	KEpsilonComputeProduction(vel=vel, k=k, eps=eps, prod=prod, nuT=nuT, strain=strain, pscale=prodMult) 
	KEpsilonSources(k=k, eps=eps, prod=prod)
	
	if enableDiffuse:
		KEpsilonGradientDiffusion(k=k, eps=eps, vel=vel, nuT=nuT, sigmaU=10.0);

	# base solver
	advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2, clampMode=1)
	setWallBcs(flags=flags, vel=vel)
	setInflowBcs(vel=vel,dir='xXyYzZ',value=velInflow)
	solvePressure(flags=flags, vel=vel, pressure=pressure, cgMaxIterFac=0.5)
	setWallBcs(flags=flags, vel=vel)
	setInflowBcs(vel=vel,dir='xXyYzZ',value=velInflow)
	
	s.step()

	if 1 and getVisualSetting() and (t%getVisualSetting()==0):
		#maxv = (k.getMax()); print "Max k %f \n"%(maxv)
		projectPpmFull( k, '%s_%04d.ppm' % (sys.argv[0],t/getVisualSetting()) , 0, 20.0 );
   
# check final state
doTestGrid( sys.argv[0],"k"    , s, k    , threshold=0.00001 , thresholdStrict=1e-10 )
doTestGrid( sys.argv[0],"eps"  , s, eps  , threshold=0.00001 , thresholdStrict=1e-10 )
doTestGrid( sys.argv[0],"vel"  , s, vel  , threshold=0.0001  , thresholdStrict=1e-10 )

