#
# tempoGAN: A Temporally Coherent, Volumetric GAN for Super-resolution Fluid Flow
# Copyright 2018 You Xie, Erik Franz, Mengyu Chu, Nils Thuerey
#
# Plume data generation, 3D
#

from manta import *
import os, shutil, math, sys, time
from datetime import datetime
sys.path.append("../tools")
import paramhelpers as ph
import uniio

outputpath = "../3ddata_sim_test/" # default target, but actual path set by getNextSimPath below 
simNo           =     3006
outputpath,simNo = ph.getNextSimPath(simNo, outputpath)
basePath        =     ph.getParam( "basePath",        "../data/"        )
savedata		= int(ph.getParam( "savedata",        1))>0

def cubeSetup( len, setupFactor, dim ):
	res = int(setupFactor*len)
	gs = vec3(res,res,res)
	if (dim==2): gs.z = 1  # 2D test

	mainDt = 0.5
	s = FluidSolver(name='main', gridSize = gs, dim = dim)
	s.timestep = mainDt/int(setupFactor)

	vel	 = s.create(MACGrid)
	flags   = s.create(FlagGrid)
	density = s.create(RealGrid)
	pressure= s.create(RealGrid)
	vorticity= s.create(Vec3Grid)#vorticity
	norm     = s.create(RealGrid)

	bWidth=1 * setupFactor
	flags.initDomain(boundaryWidth=bWidth)
	flags.fillGrid()
	setOpenBound(flags, bWidth,'xXyYzZ',FlagOutflow|FlagEmpty) 

	# noise field
	noise = s.create(NoiseField, fixedSeed=265, loadFromFile=True)
	noise.posScale = vec3(15)
	noise.clamp = True
	noise.clampNeg = -0.5
	noise.clampPos = 1.5
	noise.valScale = 2
	noise.valOffset = 0.75
	noise.timeAnim  = 0.4 * setupFactor

	orgs = vec3 (int(gs.x / setupFactor), int(gs.y / setupFactor), int(gs.z / setupFactor))
	upz  = vec3(0, 1, 0) * setupFactor
	cpos = vec3(0.5,0.1,0.5) * orgs * setupFactor
	radius = 0.14 * orgs.x * setupFactor

	source   = s.create(Cylinder, center=cpos, radius=radius, z=upz)
	#sourceVel= s.create(Cylinder, center=gs*cpos, radius=res*0.18, z=gs*vec3(0, 0.028, 0))
	#velInflow 	= vec3(0.002, 0.002, 0)
	
	return s, noise, source, vel, flags, density, pressure, bWidth, vorticity, norm
	

def cubeStep( vel, flags, density, pressure, bWidth ):
	advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2, clampMode=2, orderSpace = 1, openBounds=True, boundaryWidth=bWidth, strength=0.8)
	advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2, clampMode=2, openBounds=True, boundaryWidth=bWidth, strength=0.8)

	setWallBcs(flags=flags, vel=vel)
	addBuoyancy(density=density, vel=vel, gravity=vec3(0,-7e-4,0), flags=flags)
	#vorticityConfinement( vel=vel, flags=flags, strength=0.001 )
	solvePressure(flags=flags, vel=vel, pressure=pressure ,  cgMaxIterFac=1.0, cgAccuracy=0.01 )
	setWallBcs(flags=flags, vel=vel)
	#dissipate(density, 0.05, 0.005);


simId = 3006
showGui = 0
dim = 3 # only supports 3D

baseRes = 64
setupFactor = 4.0

sm, noise, source, vel, flags, density, pressure, bWidth, _, _ = cubeSetup( baseRes, setupFactor, 3 )
targs, target_noise, target_source, target_vel, target_flags, target_density, target_pressure, target_bWidth, target_vorticity, target_norm = cubeSetup( baseRes, 1.0, 3 )
# targs.timestep = sm.timestep ?

dummy = targs.create(MACGrid)
blurden  = sm.create(RealGrid)
blurvel  = sm.create(MACGrid)

if (showGui and GUI):
	gui = Gui()
	gui.setCamPos(0., 0., -1.3)
	gui.show()

# main loop

if savedata: 
	folderNo = simNo
	outputpath,simNo = ph.getNextSimPath(simNo, basePath)

	# add some more info for json file
	ph.paramDict["res"]   = baseRes
	ph.paramDict["fac"]   = setupFactor
	ph.paramDict["simNo"] = simNo
	ph.paramDict["dim"]   = dim
	ph.paramDict["type"]  = "smoke"
	ph.paramDict["name"]  = "plume3006"
	ph.paramDict["version"] = printBuildInfo()
	ph.paramDict["creation_date"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S") 
	ph.writeParams(outputpath + "description.json") # export sim parameters 

	sys.stdout = ph.Logger(outputpath)
	print("Called on machine '"+ os.uname()[1] +"' with: " + str(" ".join(sys.argv) ) )
	print("Saving to "+outputpath+", "+str(simNo))
	uniio.backupFile(__file__, outputpath)

#if(not os.path.exists(outputpath)):
#	os.mkdir( outputpath )
#nowtstr = time.strftime('%Y%m%d_%H%M%S')
#logfilepath = '%sscene_%s.py' % (outputpath, nowtstr)
#shutil.copyfile( sys.argv[0], logfilepath )

for t in range(400):

	mantaMsg('\nFrame %i, simulation time %f' % (targs.frame, targs.timeTotal))

	for substep in range ( int(targs.timestep / sm.timestep) ):
		cubeStep( vel, flags, density, pressure, bWidth )
		if( t < 160 ):
			densityInflow( flags=flags, density=density, noise=noise, shape=source, scale=1, sigma=0.5 )
		density.setBound(0, 5)
		sm.step()
	
	# copy to target
	if 1:
		blurSig = float(setupFactor) / 3.544908 # 3.544908 = 2 * sqrt( PI )
		blurRealGrid( density, blurden, blurSig)
		interpolateGrid( target=target_density, source=blurden )

		blurMacGrid( vel, blurvel, blurSig)
		interpolateMACGrid( target=target_vel, source=blurvel )
		target_vel.multConst( vec3(targs.timestep / sm.timestep / setupFactor) )
	else:
		cubeStep( target_vel, target_flags, target_density, target_pressure, target_bWidth )
		if( t < 320 ):
			densityInflow( flags=target_flags, density=target_density, noise=target_noise, shape=target_source, scale=1, sigma=0.5 )
		target_density.setBound(0, 5)
	# save
	if savedata and t%2==0: 
		frameNr = t / 2
		framedir = "frame_%04d" % frameNr
		#os.mkdir( outputpath + framedir )
		computeVorticity(vel = target_vel,vorticity = target_vorticity,norm = target_norm)#vorticity
		target_vorticity.save('%s/vorticity_low_%04d.uni' % (outputpath,frameNr))
		target_vel.save("%s/velocity_low_%04d.uni" % (outputpath,frameNr) )
		target_density.save("%s/density_low_%04d.uni" % (outputpath,frameNr) )
		density.save("%s/density_high_%04d.uni" % (outputpath,frameNr) )
		if(0):
			projectPpmFull( density, '%ssource_high_%04d.ppm' % (outputpath,frameNr), 0, 3.0 )
			projectPpmFull( target_density, '%ssource_low_%04d.ppm' % (outputpath, frameNr), 0, 3.0 )
	if (showGui and GUI):
		gui.screenshot( 'plume_%04d.png' % frameNr );
	
	targs.step()


