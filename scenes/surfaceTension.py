#
# Level-set surface tracking with surface tension
# 
from manta import *

# solver params
surfaceTension = 0.1
dim = 3
res = 40
gs = Vec3(res,res,res)
if (dim==2):
	gs.z=1
s = Solver(name='main', gridSize = gs, dim=dim)
s.timestep  = 0.25

# scene file params
ghostFluid  = True
accuracy    = 5e-4
# using fast marching is more accurate, esp. important for surface tension
useMarching = True
lsOrder     = 1

# prepare grids and particles
curv     = s.create(RealGrid)
flags = s.create(FlagGrid)
vel = s.create(MACGrid)
pressure = s.create(RealGrid)
mesh = s.create(Mesh)
phiBackup = s.create(LevelsetGrid)

# scene setup
bWidth=1
flags.initDomain(boundaryWidth=bWidth)
if 0:
	basin = Box( parent=s, p0=gs*Vec3(0,0,0), p1=gs*Vec3(1,0.2,1))
	drop  = Sphere( parent=s , center=gs*Vec3(0.5,0.5,0.5), radius=res*0.125)
	phi = basin.computeLevelset()
	phi.join(drop.computeLevelset())
else:
	fluidbox = Box( parent=s, p0=gs*vec3(0.25,0.25,0.25), p1=gs*vec3(0.75,0.75,0.75)) # centered falling block
	phi = fluidbox.computeLevelset()
flags.updateFromLevelset(phi)

		
if (GUI):
	gui = Gui()
	gui.show()
	gui.pause()
	

#main loop
for t in range(1000):
	mantaMsg('\nFrame %i, simulation time %f' % (s.frame, s.timeTotal))
	
	# update and advect levelset
	if useMarching:
		phi.reinitMarching(flags=flags, velTransport=vel) 
	else:
		extrapolateLsSimple(phi=phi, distance=5, inside=False)
		extrapolateLsSimple(phi=phi, distance=5, inside=True )
		extrapolateMACSimple( flags=flags, vel=vel, distance=5 )

	advectSemiLagrange(flags=flags, vel=vel, grid=phi, order=lsOrder) 
	phi.setBoundNeumann(bWidth)
	flags.updateFromLevelset(phi)

	advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2)
	
	setWallBcs(flags=flags, vel=vel)

	getLaplacian(laplacian=curv, grid=phi)
	solvePressure(flags=flags, vel=vel, pressure=pressure, phi=phi, curv=curv, surfTens=surfaceTension, cgAccuracy=accuracy)
	
	if (dim==3):
		phi.createMesh(mesh)
		#mesh.save('phi%04d.bobj.gz' % t)
	
	s.step()



