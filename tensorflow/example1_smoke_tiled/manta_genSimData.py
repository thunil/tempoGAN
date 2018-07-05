#******************************************************************************
#
# Double sim data generation
# (note, no blur when transferring hi->lo at the moment)
#
#******************************************************************************
from manta import *
import os, shutil, math, sys, time
import numpy as np
sys.path.append("../tools")
import paramhelpers as ph

# Main params  ----------------------------------------------------------------------#
steps    = 200
savedata = True
saveppm  = False
simNo    = 1000  # start ID
showGui  = 0
basePath = '../data/'
npSeedstr   = "-1"
res      = 64

# debugging
#steps = 50       # shorter test
#savedata = False # debug , dont write...
#showGui  = 1

basePath        =     ph.getParam( "basePath",        basePath  )
npSeedstr       =     ph.getParam( "npSeed"  ,        npSeedstr )
simNo           =     int(ph.getParam( "simNo" ,      simNo     ))
res             =     int(ph.getParam( "res"   ,      res       ))
steps           =     int(ph.getParam( "steps" ,      steps     ))
npSeed          =     int(npSeedstr)
ph.checkUnusedParams()

# Scene settings  ---------------------------------------------------------------------#
setDebugLevel(1)

# Solver params  ----------------------------------------------------------------------#
dim    = 2 
offset = 20
interval = 1

scaleFactor = 4

sm_gs = vec3(res,res,res)
xl_gs = sm_gs * float(scaleFactor)
if (dim==2):  xl_gs.z = sm_gs.z = 1  # 2D

#buoy    = vec3(0,-9e-4,0)
buoy    = vec3(0,-1e-3,0)
xl_buoy = buoy * vec3(1./scaleFactor)
velOffset    = vec3(0.)
xl_velOffset = vec3(0.)

# wlt Turbulence input fluid
sm = Solver(name='smaller', gridSize = sm_gs, dim=dim)
sm.timestep = 0.5

# wlt Turbulence output fluid
xl = Solver(name='larger', gridSize = xl_gs, dim=dim)
xl.timestep = sm.timestep 

# Simulation Grids  -------------------------------------------------------------------#
flags    = sm.create(FlagGrid)
vel      = sm.create(MACGrid)
velTmp   = sm.create(MACGrid)
density  = sm.create(RealGrid)
pressure = sm.create(RealGrid)

xl_flags   = xl.create(FlagGrid)
xl_vel     = xl.create(MACGrid)
xl_velTmp  = xl.create(MACGrid)
xl_blurvel = xl.create(MACGrid)
xl_density = xl.create(RealGrid)
xl_blurden = xl.create(RealGrid)
xl_pressure= xl.create(RealGrid)

# open boundaries
bWidth=1
flags.initDomain(boundaryWidth=bWidth)
flags.fillGrid()
xl_flags.initDomain(boundaryWidth=bWidth)
xl_flags.fillGrid()

setOpenBound(flags,    bWidth,'yY',FlagOutflow|FlagEmpty) 
setOpenBound(xl_flags, bWidth,'yY',FlagOutflow|FlagEmpty) 

# inflow sources ----------------------------------------------------------------------#
if(npSeed>0): np.random.seed(npSeed)

# init random density
noise    = []
sources  = []

noiseN = 12
nseeds = np.random.randint(10000,size=noiseN)

cpos = vec3(0.5,0.5,0.5)

randoms = np.random.rand(noiseN, 8)
for nI in range(noiseN):
	noise.append( sm.create(NoiseField, fixedSeed= int(nseeds[nI]), loadFromFile=True) )
	noise[nI].posScale = vec3( res * 0.1 * (randoms[nI][7] + 1) )
	noise[nI].clamp = True
	noise[nI].clampNeg = 0
	noise[nI].clampPos = 1.0
	noise[nI].valScale = 1.0
	noise[nI].valOffset = -0.01 # some gap
	noise[nI].timeAnim = 0.3
	noise[nI].posOffset = vec3(1.5)
	
	# random offsets
	coff = vec3(0.4) * (vec3( randoms[nI][0], randoms[nI][1], randoms[nI][2] ) - vec3(0.5))
	radius_rand = 0.035 + 0.035 * randoms[nI][3]
	upz = vec3(0.95)+ vec3(0.1) * vec3( randoms[nI][4], randoms[nI][5], randoms[nI][6] )
	if(dim == 2): 
		coff.z = 0.0
		upz.z = 1.0
	if( nI%2 == 0 ):
		sources.append(xl.create(Cylinder, center=xl_gs*(cpos+coff), radius=xl_gs.x*radius_rand, \
			z=xl_gs*radius_rand*upz))
	else:
		sources.append(xl.create(Sphere, center=xl_gs*(cpos+coff), radius=xl_gs.x*radius_rand, scale=upz))
		
	print (nI, "centre", xl_gs*(cpos+coff), "radius", xl_gs.x*radius_rand, "other", upz )
	
	densityInflow( flags=xl_flags, density=xl_density, noise=noise[nI], shape=sources[nI], scale=1.0, sigma=1.0 )

# init random velocity
Vrandom = np.random.rand(3)
v1pos = vec3(0.7 + 0.4 *(Vrandom[0] - 0.5) ) #range(0.5,0.9) 
v2pos = vec3(0.3 + 0.4 *(Vrandom[1] - 0.5) ) #range(0.1,0.5)
vtheta = Vrandom[2] * math.pi * 0.5
velInflow = 0.04 * vec3(math.sin(vtheta), math.cos(vtheta), 0)

if(dim == 2):
	v1pos.z = v2pos.z = 0.5
	xl_sourcV1 = xl.create(Sphere, center=xl_gs*v1pos, radius=xl_gs.x*0.1, scale=vec3(1))
	xl_sourcV2 = xl.create(Sphere, center=xl_gs*v2pos, radius=xl_gs.x*0.1, scale=vec3(1))
	xl_sourcV1.applyToGrid( grid=xl_vel , value=(-velInflow*float(xl_gs.x)) )
	xl_sourcV2.applyToGrid( grid=xl_vel , value=( velInflow*float(xl_gs.x)) )
elif(dim == 3):
	VrandomMore = np.random.rand(3)
	vtheta2 = VrandomMore[0] * math.pi * 0.5
	vtheta3 = VrandomMore[1] * math.pi * 0.5
	vtheta4 = VrandomMore[2] * math.pi * 0.5
	for dz in range(1,10,1):
		v1pos.z = v2pos.z = (0.1*dz)
		vtheta_xy = vtheta *(1.0 - 0.1*dz ) + vtheta2 * (0.1*dz)
		vtheta_z  = vtheta3 *(1.0 - 0.1*dz ) + vtheta4 * (0.1*dz)
		velInflow = 0.04 * vec3( math.cos(vtheta_z) * math.sin(vtheta_xy), math.cos(vtheta_z) * math.cos(vtheta_xy),  math.sin(vtheta_z))
		xl_sourcV1 = xl.create(Sphere, center=xl_gs*v1pos, radius=xl_gs.x*0.1, scale=vec3(1))
		xl_sourcV2 = xl.create(Sphere, center=xl_gs*v2pos, radius=xl_gs.x*0.1, scale=vec3(1))
		xl_sourcV1.applyToGrid( grid=xl_vel , value=(-velInflow*float(xl_gs.x)) )
		xl_sourcV2.applyToGrid( grid=xl_vel , value=( velInflow*float(xl_gs.x)) )

blurSig = float(scaleFactor) / 3.544908 # 3.544908 = 2 * sqrt( PI )
xl_blurden.copyFrom( xl_density )
# todo blurRealGrid( xl_density, xl_blurden, blurSig)
interpolateGrid( target=density, source=xl_blurden )

xl_blurvel.copyFrom( xl_vel )
# todo blurMacGrid( xl_vel, xl_blurvel, blurSig)
interpolateMACGrid( target=vel, source=xl_blurvel )
vel.multConst( vec3(1./scaleFactor) )

printBuildInfo()

# Setup UI ---------------------------------------------------------------------#
if (showGui and GUI):
	gui=Gui()
	gui.show()
	gui.pause()

if savedata:
	simPath, simNo = ph.getNextSimPath(simNo, basePath)
	sys.stdout = ph.Logger(simPath)

t = 0
resetN = 20

# main loop --------------------------------------------------------------------#
while t < steps+offset:
	curt = t * sm.timestep
	mantaMsg( "Current time t: " + str(curt) +" \n" )
	
	newCentre = calcCenterOfMass(xl_density)
	xl_velOffset = xl_gs*float(0.5) - newCentre
	xl_velOffset = xl_velOffset * (1./ xl.timestep)
	velOffset = xl_velOffset * (1./ float(scaleFactor))
	
	#velOffset = xl_velOffset = vec3(0.0) # re-centering off
	if(dim == 2):
		xl_velOffset.z = velOffset.z = 0.0
	
	# high res fluid
	advectSemiLagrange(flags=xl_flags, vel=xl_velTmp, grid=xl_vel, order=2, openBounds=True, boundaryWidth=bWidth)
	setWallBcs(flags=xl_flags, vel=xl_vel)
	addBuoyancy(density=xl_density, vel=xl_vel, gravity=buoy , flags=xl_flags)
	if 1 and ( t< offset ): 
		vorticityConfinement( vel=xl_vel, flags=xl_flags, strength=0.05 )
	solvePressure(flags=xl_flags, vel=xl_vel, pressure=xl_pressure ,  cgMaxIterFac=10.0, cgAccuracy=0.0001 )
	setWallBcs(flags=xl_flags, vel=xl_vel)
	xl_velTmp.copyFrom( xl_vel )
	xl_velTmp.addConst( xl_velOffset )
	if( dim == 2 ):
		xl_vel.multConst( vec3(1.0,1.0,0.0) )
		xl_velTmp.multConst( vec3(1.0,1.0,0.0) )
	advectSemiLagrange(flags=xl_flags, vel=xl_velTmp, grid=xl_density, order=2, openBounds=True, boundaryWidth=bWidth)
	xl_density.clamp(0.0, 2.0)

	# low res fluid, velocity
	if( t % resetN == 0) :
		xl_blurvel.copyFrom( xl_vel )
		# optional blurMacGrid( xl_vel, xl_blurvel, blurSig)
		interpolateMACGrid( target=vel, source=xl_blurvel )
		vel.multConst( vec3(1./scaleFactor) )
	else:
		advectSemiLagrange(flags=flags, vel=velTmp, grid=vel, order=2, openBounds=True, boundaryWidth=bWidth)
		setWallBcs(flags=flags, vel=vel)
		addBuoyancy(density=density, vel=vel, gravity=xl_buoy , flags=flags)
		if 1 and ( t< offset ): 
			vorticityConfinement( vel=vel, flags=flags, strength=0.05/scaleFactor )
		solvePressure(flags=flags, vel=vel, pressure=pressure , cgMaxIterFac=10.0, cgAccuracy=0.0001 )
		setWallBcs(flags=flags, vel=vel)

	velTmp.copyFrom(vel)
	velTmp.addConst( velOffset )

	# low res fluid, density
	if( t % resetN == 0) :
		xl_blurden.copyFrom( xl_density )
		# optional blurRealGrid( xl_density, xl_blurden, blurSig)
		interpolateGrid( target=density, source=xl_blurden )
	else:
		advectSemiLagrange(flags=flags, vel=velTmp, grid=density, order=2, openBounds=True, boundaryWidth=bWidth)
		density.clamp(0.0, 2.0)

	# save low and high res
	# save all frames
	if savedata and t>=offset and (t-offset)%interval==0:
		tf = (t-offset)/interval
		framePath = simPath + 'frame_%04d/' % tf
		os.makedirs(framePath)
		density.save(framePath + 'density_low_%04d_%04d.uni' % (simNo, tf))
		vel.save(framePath + 'vel_low_%04d_%04d.uni' % (simNo, tf))
		xl_density.save(framePath + 'density_high_%04d_%04d.uni' % (simNo, tf))
		if(saveppm):
			projectPpmFull( xl_density, simPath + 'density_high_%04d_%04d.ppm' % (simNo, tf), 0, 1.0 )
			projectPpmFull( density, simPath + 'density_low_%04d_%04d.ppm' % (simNo, tf), 0, 1.0 )

	sm.step()
	#gui.screenshot( 'outLibt1_%04d.png' % t )

	xl.step()
	t = t+1

