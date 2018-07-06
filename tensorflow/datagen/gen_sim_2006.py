#
# tempoGAN: A Temporally Coherent, Volumetric GAN for Super-resolution Fluid Flow
# Copyright 2018 You Xie, Erik Franz, Mengyu Chu, Nils Thuerey
#
# Plume data generation, 2D
# 
from manta import *
import os, shutil, math, sys
import numpy as np
sys.path.append("../tools")
import paramhelpers as ph

simId = 2006
simPath = '../2ddata_sim/'
simPath,simId = ph.getNextSimPath(simId, simPath)

# how much to reduce target sim size
targetFac = 0.25
savenpz = 1

# source solver params
dim = 2
#res = 128
res = 512
gs  = vec3(res,int(1.0*res),res)
if (dim==2): gs.z = 1  # 2D

sm = Solver(name='main', gridSize = gs, dim=dim)
sm.timestep = 1.5 
sm.timestep = 0.75 

# inflow noise field
noise = NoiseField( parent=sm, fixedSeed=265, loadFromFile=True)
noise.posScale = vec3(24)
noise.clamp = True
noise.clampNeg = 0
noise.clampPos = 2
noise.valScale = 1
noise.valOffset = 0.075
noise.timeAnim = 0.3
noise.timeAnim = 0.5

cylWidth = 0.13
source = Cylinder(parent=sm, center=gs*vec3(0.5,0.1,0.5), radius=res*cylWidth, z=gs*vec3(0, 0.04, 0))


# target solver, recompute sizes...

target_gs = vec3(targetFac*gs.x,targetFac*gs.y,targetFac*gs.z)
if (dim==2): target_gs.z = 1  # 2D
targs = Solver(name='target', gridSize = target_gs, dim=dim)
targs.timestep = sm.timestep 

dummy = targs.create(MACGrid)
target_flags   = targs.create(FlagGrid)
target_vel     = targs.create(MACGrid)
target_density = targs.create(RealGrid)

target_flags.initDomain()
target_flags.fillGrid()

target_source = Cylinder(parent=targs, center=target_gs*vec3(0.5,0.1,0.5), radius=res*targetFac*cylWidth, z=target_gs*vec3(0, 0.04, 0))

if savenpz:
	arR = np.zeros([int(gs.z), int(gs.y), int(gs.x), 1])
	arV = np.zeros([int(gs.z), int(gs.y), int(gs.x), 3])
	target_arR = np.zeros([int(target_gs.z), int(target_gs.y), int(target_gs.x), 1])
	target_arV = np.zeros([int(target_gs.z), int(target_gs.y), int(target_gs.x), 3])

# allocate other grids
flags    = sm.create(FlagGrid)
vel      = sm.create(MACGrid)
density  = sm.create(RealGrid)
pressure = sm.create(RealGrid)
blurden  = sm.create(RealGrid)
blurvel  = sm.create(MACGrid)

bWidth=0
flags.initDomain(boundaryWidth=bWidth)
flags.fillGrid() 
setOpenBound(flags,bWidth,'yY',FlagOutflow|FlagEmpty) 


if (GUI):
	gui = Gui()
	gui.setCamPos(0., 0., -1.3)
	gui.show()

# main loop
for t in range(400):
	mantaMsg('\nFrame %i, simulation time %f' % (sm.frame, sm.timeTotal))
		
	advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2, clampMode=2 )    
	advectSemiLagrange(flags=flags, vel=vel, grid=vel,     order=2, clampMode=2 , openBounds=True, boundaryWidth=bWidth )
	
	applyInflow=False
	if (sm.timeTotal>=0 and sm.timeTotal<150.):
		densityInflow( flags=flags, density=density, noise=noise, shape=source, scale=1, sigma=0.5 )
		applyInflow=True
	
	setWallBcs(flags=flags, vel=vel)    
	#addBuoyancy(density=density, vel=vel, gravity=vec3(0,-1e-3,0), flags=flags)
	addBuoyancy(density=density, vel=vel, gravity=vec3(0,-2e-4,0), flags=flags)

	#vorticityConfinement( vel=vel, flags=flags, strength=0.1 )
	
	solvePressure(flags=flags, vel=vel, pressure=pressure , cgMaxIterFac=1.0, cgAccuracy=0.01 )
	setWallBcs(flags=flags, vel=vel)
	
	sm.step()
	
	# copy to target
	if 1:
		blurSig = float(1./targetFac) / 3.544908 # 3.544908 = 2 * sqrt( PI )
		blurRealGrid( density, blurden, blurSig)
		interpolateGrid( target=target_density, source=blurden )

		blurMacGrid( vel, blurvel, blurSig)
		interpolateMACGrid( target=target_vel, source=blurvel )
		target_vel.multConst( vec3(targetFac) )

	# save
	if 0 and t%2==0: 
		frameNr = t / 2
		framedir = "frame_%04d" % frameNr
		os.mkdir( framedir )

		target_vel.save("%s/vel_low_%04d_%04d.uni" % (framedir,simId,frameNr) )
		target_density.save("%s/density_low_%04d_%04d.uni" % (framedir,simId,frameNr) )
		density.save("%s/density_high_%04d_%04d.uni" % (framedir,simId,frameNr) )
	
		#gui.screenshot( 'plume_%04d.png' % frameNr );

	if savenpz and t%2==0: 
		tf = t / 2
		print("Writing NPZs for frame %d"%tf)
		copyGridToArrayReal( target=target_arR, source=target_density )
		np.savez_compressed( simPath + 'density_low_%04d.npz' % (tf), target_arR )
		copyGridToArrayVec3( target=target_arV, source=target_vel )
		np.savez_compressed( simPath + 'velocity_low_%04d.npz' % (tf), target_arV )
		copyGridToArrayReal( target=arR, source=density )
		np.savez_compressed( simPath + 'density_high_%04d.npz' % (tf), arR )
		copyGridToArrayVec3( target=arV, source=vel )
		np.savez_compressed( simPath + 'velocity_high_%04d.npz' % (tf), arV )


	targs.step()    


