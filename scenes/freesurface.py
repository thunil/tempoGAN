#
# Simple example of a free-surface simulation with level-set
# (Optionally, demos outflow boundaries)
# 
from manta import *

# solver params
dim = 3
res = 64
#res = 128 
gs = Vec3(res,res,res)
if (dim==2):
	gs.z=1
s = Solver(name='main', gridSize = gs, dim=dim)
s.timestep  = 0.15 

# scene file params
ghostFluid  = True
doOpen      = False
accuracy    = 5e-4
# using fast marching is more accurate, but currently causes asymmetries
useMarching = False

# prepare grids and particles
phi = s.create(LevelsetGrid)
flags = s.create(FlagGrid)
vel = s.create(MACGrid)
pressure = s.create(RealGrid)
mesh = s.create(Mesh)

# scene setup
bWidth=1
flags.initDomain(boundaryWidth=bWidth)
basin = Box( parent=s, p0=gs*Vec3(0,0,0), p1=gs*Vec3(1,0.2,1))
#basin = Box( parent=s, p0=gs*Vec3(0,0,0), p1=gs*Vec3(0.2,1,1)) # Y
drop  = Sphere( parent=s , center=gs*Vec3(0.5,0.5,0.5), radius=res*0.125)
phi.setConst(1e10)
phi.join(basin.computeLevelset())
phi.join(drop.computeLevelset())
flags.updateFromLevelset(phi)

# optionally, enable open boundaries here and below...
if doOpen:
	setOpenBound(flags,bWidth,'xXzZ',FlagOutflow|FlagEmpty) 
		
if (GUI):
	gui = Gui()
	gui.show()
	gui.pause()
	

#main loop
for t in range(1000):
	mantaMsg('\nFrame %i, simulation time %f' % (s.frame, s.timeTotal))
	
	# update and advect levelset
	if useMarching:
		phi.reinitMarching(flags=flags, velTransport=vel, ignoreWalls=True) 
	else:
		extrapolateLsSimple(phi=phi, distance=5, inside=False)
		extrapolateLsSimple(phi=phi, distance=5, inside=True )
		extrapolateMACSimple( flags=flags, vel=vel, distance=5 )

	advectSemiLagrange(flags=flags, vel=vel, grid=phi, order=2, clampMode=2) 

	phi.setBound(bWidth, 1.) # enforce outside values at border
	if doOpen:
		resetOutflow(flags=flags,phi=phi) # open boundaries
	flags.updateFromLevelset(phi)
	
	# velocity self-advection
	advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2, openBounds=doOpen, boundaryWidth=bWidth )
	addGravity(flags=flags, vel=vel, gravity=Vec3(0,-0.025,0))
	#addGravity(flags=flags, vel=vel, gravity=Vec3(-0.025,0,0)) # Y
	
	# pressure solve
	setWallBcs(flags=flags, vel=vel)
	if ghostFluid:
		solvePressure(flags=flags, vel=vel, pressure=pressure, cgMaxIterFac=0.5, cgAccuracy=accuracy, phi=phi )
	else:
		solvePressure(flags=flags, vel=vel, pressure=pressure, cgMaxIterFac=0.5, cgAccuracy=accuracy)
	
	# note: these meshes are created by fast marching only, should smooth
	#       geometry and normals before rendering (only in 3D for now)
	if (dim==3):
		phi.createMesh(mesh)
		#mesh.save('phi%04d.bobj.gz' % t)
	
	s.step()
	#gui.screenshot( 'freesurface_%04d.png' % t );


