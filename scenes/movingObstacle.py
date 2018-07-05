#
# Simple example scene with a moving obstacle
#
from manta import *

dim = 2
res = 50
gs  = vec3(res,res, 1 if dim==2 else res )
s   = Solver(name='main', gridSize = gs, dim=dim)

# allocate grids
flags   = s.create(FlagGrid)
vel     = s.create(MACGrid)
density = s.create(RealGrid)
pressure = s.create(RealGrid)
obsVel  = s.create(MACGrid)

bWidth = 1
flags.initDomain(boundaryWidth=bWidth) 
flags.fillGrid()
setOpenBound(flags, bWidth,'yY',FlagOutflow|FlagEmpty) 

source = Box( parent=s, p0=gs*vec3(0.45,0.1,0.1), p1=gs*vec3(0.55,0.9,0.9))
source.applyToGrid(grid=density, value=1)

# init obstacle properties
obsPos = vec3(0.2,0.4,0.5)
obsVelVec = vec3(0.6,0.2,0.0) * (1./100.) * float(res) # velocity in grid units for 100 steps
obsSize = 0.1
obsVel.setConst(obsVelVec)
obsVel.setBound(value=Vec3(0.), boundaryWidth=bWidth+1) # make sure walls are static
obs = "dummy"; phiObs = "dummy2"

if (GUI):
	gui = Gui()
	gui.show( True ) 
	#gui.pause()
	
#main loop
for t in range(400):
	mantaMsg('\nFrame %i' % (s.frame))
		
	advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2) 
	advectSemiLagrange(flags=flags, vel=vel, grid=vel,     order=2, openBounds=True, boundaryWidth=bWidth)
	resetOutflow(flags=flags,real=density) 

	if t<=100:
		flags.initDomain(boundaryWidth=bWidth) 
		flags.fillGrid()
		setOpenBound(flags, bWidth,'yY',FlagOutflow|FlagEmpty) 

		del obs, phiObs
		# use sphere or box?
		if 1:
			obs = Sphere( parent=s, center=gs*obsPos + float(t) * obsVelVec, radius=res*obsSize)
		else:
			obs = Box( parent=s, p0=gs*vec3(0.15-obsSize*0.5,0.2,0.4) + float(t) * obsVelVec, \
			                     p1=gs*vec3(0.15+obsSize*0.5,0.5,0.6) + float(t) * obsVelVec)
		phiObs = obs.computeLevelset()

		setObstacleFlags(flags=flags, phiObs=phiObs) 
		flags.fillGrid()

		obs.applyToGrid(grid=density, value=0.) # clear smoke inside
	elif t==101:
		# stop moving
		obsVel.setConst(Vec3(0.))

	setWallBcs(flags=flags, vel=vel, phiObs=phiObs, obvel=obsVel)
	solvePressure(flags=flags, vel=vel, pressure=pressure)

	s.step()
