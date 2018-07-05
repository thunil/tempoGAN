#
# Test advection modes (and basic noise field)
# 

import sys
from manta import *
from helperInclude import *

res = 36
gs = vec3(res,res,res)
s = Solver(name='main', gridSize = gs, dim=3)
# use non standard timestep...
s.timestep = 1.2

density  = s.create(RealGrid)
vgrid      = s.create(VecGrid)
mgrid      = s.create(MACGrid)

flags    = s.create(FlagGrid)
vel      = s.create(MACGrid)
pressure = s.create(RealGrid)

flags.initDomain()
flags.fillGrid()

if 0 and (GUI):
	gui = Gui()
	gui.show()
	gui.pause()

velSource1 = s.create(Box, p0=gs*vec3(0.25,0.30,0.35), p1=gs*vec3(0.45,0.50,0.55) )
velSource2 = s.create(Box, p0=gs*vec3(0.75,0.70,0.65), p1=gs*vec3(0.90,0.85,0.85) )
dSource    = s.create(Box, p0=gs*vec3(0.1), p1=gs*vec3(0.9) )

# simple noise field
noise = s.create(NoiseField, loadFromFile=True )
noise.posScale = vec3(40)
noise.valScale = 2
noise.valOffset = -0.5

# ============================

# get a somewhat interesting flow field
vel.setConst( vec3(0,0,0) )
velSource1.applyToGrid(grid=vel, value=vec3(0.1,  2, 0.2) ) 
velSource2.applyToGrid(grid=vel, value=vec3(-0.1,-2,-0.2) ) 
setWallBcs(flags=flags, vel=vel) 
solvePressure(flags=flags, vel=vel, pressure=pressure, cgMaxIterFac=99, cgAccuracy=1e-04, zeroPressureFixing=False)


doTestGrid( sys.argv[0], "dens0init" , s, pressure , threshold=1e-04, thresholdStrict=1e-10)

def initGrids(s,v,m):
	s.setConst(0.)
	densityInflow(flags=flags, density=s, noise=noise, shape=dSource, scale=1, sigma=0.5)
	v.setConst(Vec3(0.))
	setComponent(s, v, 0)
	setComponent(s, v, 1)
	setComponent(s, v, 2)
	m.setConst(Vec3(0.))
	m.copyFrom(v)

# first order
initGrids(density, vgrid, mgrid)
for i in range(10):
	advectSemiLagrange(flags=flags, vel=vel, grid=density, order=1)
	advectSemiLagrange(flags=flags, vel=vel, grid=vgrid  , order=1)
	advectSemiLagrange(flags=flags, vel=vel, grid=mgrid  , order=1)
	s.step()

doTestGrid( sys.argv[0], "dens1" , s, pressure , threshold=1e-04, thresholdStrict=1e-10)
doTestGrid( sys.argv[0], "vgrid1" , s, pressure , threshold=1e-04, thresholdStrict=1e-10)
doTestGrid( sys.argv[0], "mgrid1" , s, pressure , threshold=1e-04, thresholdStrict=1e-10)

# macCormack
initGrids(density, vgrid, mgrid)
for i in range(10):
	advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2, clampMode=1)
	advectSemiLagrange(flags=flags, vel=vel, grid=vgrid  , order=2, clampMode=1)
	advectSemiLagrange(flags=flags, vel=vel, grid=mgrid  , order=2, clampMode=1)
	s.step()

doTestGrid( sys.argv[0], "dens2" , s, pressure , threshold=1e-04, thresholdStrict=1e-10)
doTestGrid( sys.argv[0], "vgrid2" , s, pressure , threshold=1e-04, thresholdStrict=1e-10)
doTestGrid( sys.argv[0], "mgrid2" , s, pressure , threshold=1e-04, thresholdStrict=1e-10)

# macCormack, less aggressive clamping
initGrids(density, vgrid, mgrid)
for i in range(10):
	advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2, clampMode=2)
	advectSemiLagrange(flags=flags, vel=vel, grid=vgrid  , order=2, clampMode=2)
	advectSemiLagrange(flags=flags, vel=vel, grid=mgrid  , order=2, clampMode=2)
	s.step()

doTestGrid( sys.argv[0], "dens3" , s, pressure , threshold=1e-04, thresholdStrict=1e-10)
doTestGrid( sys.argv[0], "vgrid3" , s, pressure , threshold=1e-04, thresholdStrict=1e-10)
doTestGrid( sys.argv[0], "mgrid3" , s, pressure , threshold=1e-04, thresholdStrict=1e-10)


