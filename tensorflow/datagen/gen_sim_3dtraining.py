#******************************************************************************
#
# Varying density data gen, 2d/3d
#
# 2 modes:
# - double sim mode (mode==1)
# 		Runs hi-res, then downsamples to coarse ("sm") sim in intervals, by default
# 		every frame
# - wavelet turbulence mode (mode==2)
# 		Runs low-res, then upsamples and adds WLT
#
#******************************************************************************

from manta import *
import os, shutil, math, sys, time
from datetime import datetime
import numpy as np
sys.path.append("../tools")
import paramhelpers as ph
import uniio # for backup

# Main params  ----------------------------------------------------------------------#
steps    = 120
savedata = False
saveppm  = False
simNo    = 1000  # start ID
showGui  = 1
basePath = '../data/'
npSeedstr   = "-1"
dim         = 2 
doRecenter  = False   # re-center densities

# Solver params  
res         = 64
scaleFactor = 4 
resetN      = 1

# cmd line
basePath        =     ph.getParam( "basePath",        basePath        )
npSeedstr       =     ph.getParam( "seed"    ,        npSeedstr       )
npSeed          =     int(npSeedstr)
resetN			= int(ph.getParam( "reset"   ,        resetN))
dim   			= int(ph.getParam( "dim"     ,        dim))
simMode			= int(ph.getParam( "mode"    ,        1))  # 1 = double sim, 2 = wlt
savedata		= int(ph.getParam( "savedata",        1 if savedata else 0))>0
saveppm 		= int(ph.getParam( "saveppm" ,        1 if saveppm  else 0))>0
showGui 		= int(ph.getParam( "gui"     ,        showGui))
res     		= int(ph.getParam( "res"     ,        res))
scaleFactor 	= int(ph.getParam( "fac"     ,        scaleFactor))
steps     		= int(ph.getParam( "steps"   ,        steps))
timeOffset   	= int(ph.getParam( "warmup" ,        20))    # skip certain no of steps 
ph.checkUnusedParams()

setDebugLevel(1)

if savedata: 
	folderNo = simNo
	simPath,simNo = ph.getNextSimPath(simNo, basePath)

	# add some more info for json file
	ph.paramDict["simNo"] = simNo
	ph.paramDict["type"] = "smoke"
	ph.paramDict["name"] = "gen6combined"
	ph.paramDict["version"] = printBuildInfo()
	ph.paramDict["creation_date"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S") 
	ph.writeParams(simPath + "description.json") # export sim parameters 

	sys.stdout = ph.Logger(simPath)
	print("Called on machine '"+ os.uname()[1] +"' with: " + str(" ".join(sys.argv) ) )
	print("Saving to "+simPath+", "+str(simNo))
	uniio.backupFile(__file__, simPath)

if(npSeed<0): 
	npSeed = np.random.randint(0, 2**32-1 )
print("Random seed %d" % npSeed)
np.random.seed(npSeed)

# Init solvers -------------------------------------------------------------------#
sm_gs = vec3(res,res,res) 
xl_gs = sm_gs * float(scaleFactor)
if (dim==2):  xl_gs.z = sm_gs.z = 1  # 2D

#buoy    = vec3(0,-9e-4,0)
buoyFac = 0.1 + np.random.rand()
if buoyFac < 0.4:
	buoyFac = 0. # make sure we have some sims without buoyancy
buoy    = vec3(0,-3e-3,0) * buoyFac
print("Buoyancy: " + str(buoyFac) +", factor " + str(buoyFac))
xl_buoy = buoy * vec3(1./scaleFactor)

# wlt Turbulence input fluid
sm = Solver(name='smaller', gridSize = sm_gs, dim=dim)
sm.timestep = 0.5

# wlt Turbulence output fluid
xl = Solver(name='larger', gridSize = xl_gs, dim=dim)
xl.timestep = sm.timestep 

timings = Timings()

# Simulation Grids  -------------------------------------------------------------------#
flags    = sm.create(FlagGrid)
vel      = sm.create(MACGrid)
density  = sm.create(RealGrid)
#velTmp   = sm.create(MACGrid)
tmp      = sm.create(RealGrid)
vorticityTmp= sm.create(Vec3Grid)#vorticity
norm     = sm.create(RealGrid)

xl_flags   = xl.create(FlagGrid)
xl_vel     = xl.create(MACGrid)
xl_density = xl.create(RealGrid)
xl_velTmp  = xl.create(MACGrid)
xl_tmp     = xl.create(RealGrid)

# for domain centering
velRecenter    = sm.create(MACGrid)
xl_velRecenter = xl.create(MACGrid)

# open boundaries
bWidth=1
flags.initDomain(boundaryWidth=bWidth)
flags.fillGrid()
xl_flags.initDomain(boundaryWidth=bWidth)
xl_flags.fillGrid()

setOpenBound(flags,    bWidth,'yY',FlagOutflow|FlagEmpty) 
setOpenBound(xl_flags, bWidth,'yY',FlagOutflow|FlagEmpty) 

# wavelet turbulence octaves

wltnoise = NoiseField( parent=xl, loadFromFile=True)
# scale according to lowres sim , smaller numbers mean larger vortices
wltnoise.posScale = vec3( int(1.0*sm_gs.x) ) * 0.5
wltnoise.timeAnim = 0.1

wltnoise2 = NoiseField( parent=xl, loadFromFile=True)
wltnoise2.posScale = wltnoise.posScale * 2.0
wltnoise2.timeAnim = 0.1

wltnoise3 = NoiseField( parent=xl, loadFromFile=True)
wltnoise3.posScale = wltnoise2.posScale * 2.0
wltnoise3.timeAnim = 0.1

# inflow sources ----------------------------------------------------------------------#

# init random density
sources  = []
noise    = []  # xl
sourSm   = []
noiSm    = []  # sm
inflowSrc = [] # list of IDs to use as continuous density inflows

noiseN = 12
#noiseN = 1
nseeds = np.random.randint(10000,size=noiseN)

cpos = vec3(0.5,0.5,0.5)

randoms = np.random.rand(noiseN, 10)
for nI in range(noiseN):
	#noise.append( sm.create(NoiseField, fixedSeed= int(nseeds[nI]), loadFromFile=True) )
	noise.append( xl.create(NoiseField, fixedSeed= int(nseeds[nI]), loadFromFile=True) )
	noise[nI].posScale = vec3( res * 0.1 * (randoms[nI][7] + 1) ) * ( float(scaleFactor))
	noise[nI].clamp = True
	noise[nI].clampNeg = 0
	noise[nI].clampPos = 1.0
	noise[nI].valScale = 1.0
	noise[nI].valOffset = 0.5 * randoms[nI][9]
	noise[nI].timeAnim = 0.3
	noise[nI].posOffset = vec3(1.5)
	
	noiSm.append( sm.create(NoiseField, fixedSeed= int(nseeds[nI]), loadFromFile=True) )
	#noiSm[nI].posScale = vec3( res * 0.1 * (randoms[nI][7] + 1) ) * ( float(scaleFactor))
	noiSm[nI].timeAnim = 0.3 / ( float(scaleFactor))
	noiSm[nI].posOffset = vec3(1.5) # * ( 1./float(scaleFactor))
	noiSm[nI].posScale = noise[nI].posScale
	noiSm[nI].clamp    = noise[nI].clamp    
	noiSm[nI].clampNeg = noise[nI].clampNeg 
	noiSm[nI].clampPos = noise[nI].clampPos 
	noiSm[nI].valScale = noise[nI].valScale 
	noiSm[nI].valOffset= noise[nI].valOffset
	
	# random offsets
	coff = vec3(0.4) * (vec3( randoms[nI][0], randoms[nI][1], randoms[nI][2] ) - vec3(0.5))
	radius_rand = 0.035 + 0.035 * randoms[nI][3]
	upz = vec3(0.95)+ vec3(0.1) * vec3( randoms[nI][4], randoms[nI][5], randoms[nI][6] )

	if 1 and randoms[nI][8] > 0.5: # turn into inflow?
		if coff.y > -0.2:
			coff.y += -0.4
		coff.y *= 0.5
		inflowSrc.append(nI)

	if(dim == 2): 
		coff.z = 0.0
		upz.z = 1.0
	#if( nI%2 == 0 ):
		#sources.append(xl.create(Cylinder, center=xl_gs*(cpos+coff), radius=xl_gs.x*radius_rand, z=xl_gs*radius_rand*upz))
	#else:
	sources.append(xl.create(Sphere, center=xl_gs*(cpos+coff), radius=xl_gs.x*radius_rand, scale=upz))
	sourSm.append( sm.create(Sphere, center=sm_gs*(cpos+coff), radius=sm_gs.x*radius_rand, scale=upz))
		
	print (nI, "centre", xl_gs*(cpos+coff), "radius", xl_gs.x*radius_rand, "other", upz )
	
	densityInflow( flags=xl_flags, density=xl_density, noise=noise[nI], shape=sources[nI], scale=1.0, sigma=1.0 )
	densityInflow( flags=flags, density=density, noise=noiSm[nI], shape=sourSm[nI], scale=1.0, sigma=1.0 )

# init random velocities

inivel_sources = []
inivel_vels = []
inivel_sourcesSm = []
inivel_velsSm = []
if 1: # from fluidnet
	c = 3 + np.random.randint(3) # "sub" mode
	xgs = xl_gs
	if 1:
		# 3..5 - ini vel sources
		if c==3: numFac = 1; sizeFac = 0.9;
		if c==4: numFac = 3; sizeFac = 0.7;
		if c==5: numFac = 5; sizeFac = 0.6;
		numNs = int( numFac * float(dim) )  
		for ns in range(numNs):
			p = [0.5,0.5,0.5]
			Vrand = np.random.rand(10) 
			for i in range(3):
				p[i] += (Vrand[0+i]-0.5) * 0.6
			p = Vec3(p[0],p[1],p[2])

			# org: size = 0.05+0.1*Vrand[3]
			size = ( 0.05 + 0.1*Vrand[3] ) * sizeFac

			v = [0.,0.,0.]
			for i in range(3):
				v[i] -= (Vrand[0+i]-0.5) * 0.6 * 2. # invert pos offset , towards middle
				v[i] += (Vrand[4+i]-0.5) * 0.3      # randomize a bit, parametrized for 64 base
			v = Vec3(v[0],v[1],v[2])
			v = v*0.9 # tweaking
			v = v*(1. + 0.5*Vrand[7] ) # increase by up to 50% 
			v *= float(scaleFactor)

			print( "IniVel Pos " + format(p) + ", size " + format(size) + ", vel " + format(v) )
			sourceV = xl.create(Sphere, center=xgs*p, radius=xgs.x*size, scale=vec3(1))
			inivel_sources.append(sourceV)
			inivel_vels.append(v)
			sourceVsm = sm.create(Sphere, center=sm_gs*p, radius=sm_gs.x*size, scale=vec3(1))
			inivel_sourcesSm.append(sourceVsm)
			inivel_velsSm.append(v* (1./scaleFactor))
		#print("ini velsrc 3..5 "+str(c))

# init low-res vel from hi-res
blurSig = float(scaleFactor) / 3.544908 # 3.544908 = 2 * sqrt( PI )
xl_velTmp.copyFrom( xl_vel )
blurMacGrid( xl_vel, xl_velTmp, blurSig)
interpolateMACGrid( target=vel, source=xl_velTmp )
vel.multConst( vec3(1./scaleFactor) )

printBuildInfo()

# wlt params ---------------------------------------------------------------------#

if simMode==2:
	wltStrength = 0.8
	if resetN==1:
		print("Warning!!!!!!!!!!!!!! Using resetN=1 for WLT doesnt make much sense, resetting to never")
		resetN = 99999

def calcCOM(dens):
	if doRecenter:
		newCentre = calcCenterOfMass(xl_density)
		#mantaMsg( "Current moff "+str(newCentre) )
		xl_velOffset = xl_gs*float(0.5) - newCentre
		xl_velOffset = xl_velOffset * (1./ xl.timestep) 
		velOffset = xl_velOffset * (1./ float(scaleFactor)) 
		if(dim == 2):
			xl_velOffset.z = velOffset.z = 0.0 
	else: 
		velOffset = xl_velOffset = vec3(0.0) # re-centering off

	return velOffset, xl_velOffset

# Setup UI ---------------------------------------------------------------------#
if (showGui and GUI):
	gui=Gui()
	gui.show()
	gui.pause()

t = 0
doPrinttime = False

# main loop --------------------------------------------------------------------#
while t < steps+timeOffset:
	curt = t * sm.timestep
	sys.stdout.write( "Current sim time t: " + str(curt) +" \n" )
	#density.setConst(0.); xl_density.setConst(0.); # debug reset

	if doPrinttime:
		starttime = time.time()
		print("starttime: %2f" % starttime)	

	# --------------------------------------------------------------------#
	if simMode==1: 
		velOffset , xl_velOffset = calcCOM(xl_density)

		if 1 and len(inflowSrc)>0:
			# note - the density inflows currently move with the offsets!
			for nI in inflowSrc:
				densityInflow( flags=xl_flags, density=xl_density, noise=noise[nI], shape=sources[nI], scale=1.0, sigma=1.0 )
				densityInflow( flags=flags, density=density, noise=noiSm[nI], shape=sourSm[nI], scale=1.0, sigma=1.0 )
		
		# high res fluid
		advectSemiLagrange(flags=xl_flags, vel=xl_velRecenter, grid=xl_vel, order=2, clampMode=2, openBounds=True, boundaryWidth=bWidth)
		setWallBcs(flags=xl_flags, vel=xl_vel)
		addBuoyancy(density=xl_density, vel=xl_vel, gravity=buoy , flags=xl_flags)
		if 1:
			for i in range(len(inivel_sources)):
				inivel_sources[i].applyToGrid( grid=xl_vel , value=inivel_vels[i] )
		if 1 and ( t< timeOffset ): 
			vorticityConfinement( vel=xl_vel, flags=xl_flags, strength=0.05 )

		solvePressure(flags=xl_flags, vel=xl_vel, pressure=xl_tmp ,  cgMaxIterFac=2.0, cgAccuracy=0.001, preconditioner=PcMGStatic )
		setWallBcs(flags=xl_flags, vel=xl_vel)
		xl_velRecenter.copyFrom( xl_vel )
		xl_velRecenter.addConst( xl_velOffset )
		if( dim == 2 ):
			xl_vel.multConst( vec3(1.0,1.0,0.0) )
			xl_velRecenter.multConst( vec3(1.0,1.0,0.0) )
		advectSemiLagrange(flags=xl_flags, vel=xl_velRecenter, grid=xl_density, order=2, clampMode=2, openBounds=True, boundaryWidth=bWidth)

		# low res fluid, velocity
		if( t % resetN == 0) :
			xl_velTmp.copyFrom( xl_vel )
			blurMacGrid( xl_vel, xl_velTmp, blurSig)
			interpolateMACGrid( target=vel, source=xl_velTmp )
			vel.multConst( vec3(1./scaleFactor) )
		else:
			advectSemiLagrange(flags=flags, vel=velRecenter, grid=vel, order=2, clampMode=2, openBounds=True, boundaryWidth=bWidth)
			setWallBcs(flags=flags, vel=vel)
			addBuoyancy(density=density, vel=vel, gravity=xl_buoy , flags=flags)
			if 1:
				for i in range(len(inivel_sourcesSm)):
					inivel_sourcesSm[i].applyToGrid( grid=vel , value=inivel_velsSm[i] )
			if 1 and ( t< timeOffset ): 
				vorticityConfinement( vel=vel, flags=flags, strength=0.05/scaleFactor )
			solvePressure(flags=flags, vel=vel, pressure=tmp , cgMaxIterFac=2.0, cgAccuracy=0.001, preconditioner=PcMGStatic )
			setWallBcs(flags=flags, vel=vel)

		velRecenter.copyFrom(vel)
		velRecenter.addConst( velOffset )

		# low res fluid, density
		if( t % resetN == 0) :
			xl_tmp.copyFrom( xl_density )
			blurRealGrid( xl_density, xl_tmp, blurSig)
			interpolateGrid( target=density, source=xl_tmp )
		else:
			advectSemiLagrange(flags=flags, vel=velRecenter, grid=density, order=2, clampMode=2, openBounds=True, boundaryWidth=bWidth)

	# --------------------------------------------------------------------#
	elif simMode==2: 
		# low res fluid, density
		if( t % resetN == 0) :
			xl_tmp.copyFrom( xl_density )
			blurRealGrid( xl_density, xl_tmp, blurSig)
			interpolateGrid( target=density, source=xl_tmp )
		
		advectSemiLagrange(flags=flags, vel=velRecenter, grid=density, order=2, clampMode=2)    
		if t<=1: velRecenter.copyFrom(vel); # special , center only density once, leave vel untouched 
		advectSemiLagrange(flags=flags, vel=velRecenter, grid=vel,     order=2, clampMode=2, openBounds=True, boundaryWidth=bWidth )
		
		if 1 and len(inflowSrc)>0:
			# note - the density inflows currently move with the offsets!
			for nI in inflowSrc:
				densityInflow( flags=xl_flags, density=xl_density, noise=noise[nI], shape=sources[nI], scale=1.0, sigma=1.0 )
				densityInflow( flags=flags, density=density, noise=noiSm[nI], shape=sourSm[nI], scale=1.0, sigma=1.0 )
		
		setWallBcs(flags=flags, vel=vel)    
		#addBuoyancy(density=density, vel=vel, gravity=vec3(0,-1e-4,0), flags=flags)
		addBuoyancy(density=density, vel=vel, gravity=buoy , flags=flags)
		if 1:
			for i in range(len(inivel_sourcesSm)):
				inivel_sourcesSm[i].applyToGrid( grid=vel , value=inivel_velsSm[i] )

		vorticityConfinement( vel=vel, flags=flags, strength=0.1 ) 
		#applyNoiseVec3( flags=flags, target=vel, noise=noise, scale=1 ) # just to test, add everywhere...
		
		solvePressure(flags=flags, vel=vel, pressure=tmp , cgMaxIterFac=3.0, cgAccuracy=0.001, preconditioner=PcMGStatic )
		setWallBcs(flags=flags, vel=vel)
		
		computeEnergy(flags=flags, vel=vel, energy=tmp)
		# standard weights for wavelet turbulence
		computeWaveletCoeffs(tmp)
		
		# xl solver, update up-res'ed grids ...

		# new centre of mass , from XL density
		velOffset , xl_velOffset = calcCOM(xl_density)
		xl_velOffset = velOffset  # note - hires advection does "scaleFac" substeps below! -> same offset

		if 1 and len(inflowSrc)>0:
			velOffset *= 0.5;  xl_velOffset *= 0.5;  # re-centering reduced

		# high res sim
		
		interpolateGrid( target=xl_tmp, source=tmp )
		interpolateMACGrid( source=vel, target=xl_vel )
		
		applyNoiseVec3( flags=xl_flags, target=xl_vel, noise=wltnoise, scale=wltStrength*1.0 , weight=xl_tmp)
		# manually weight and apply further octaves
		applyNoiseVec3( flags=xl_flags, target=xl_vel, noise=wltnoise2, scale=wltStrength*0.8 , weight=xl_tmp)
		applyNoiseVec3( flags=xl_flags, target=xl_vel, noise=wltnoise3, scale=wltStrength*0.8*0.8 , weight=xl_tmp)
		
		xl_velRecenter.copyFrom( xl_vel )
		xl_velRecenter.addConst( xl_velOffset )
		if( dim == 2 ):
			xl_velRecenter.multConst( vec3(1.0,1.0,0.0) )

		for substep in range(scaleFactor): 
			advectSemiLagrange(flags=xl_flags, vel=xl_velRecenter, grid=xl_density, order=2, clampMode=2)    

		velRecenter.copyFrom(vel)
		velRecenter.addConst( velOffset )
		if( dim == 2 ):
			#vel.multConst( vec3(1.0,1.0,0.0) )
			velRecenter.multConst( vec3(1.0,1.0,0.0) )
	else:
		print("Unknown sim mode!"); exit(1)

	if doPrinttime:
		endtime = time.time()
		print("endtime: %2f" % endtime)
		print("runtime: %2f" % (endtime-starttime))

	# --------------------------------------------------------------------#

	# save low and high res
	# save all frames
	if savedata and t>=timeOffset:
		tf = t-timeOffset
		density.save(simPath + 'density_low_%04d.uni' % (tf))
		velRecenter.save(simPath + 'velocity_low_%04d.uni' % (tf))
		computeVorticity(vel = velRecenter,vorticity = vorticityTmp,norm = norm)#vorticity
		vorticityTmp.save(simPath + 'vorticity_low_%04d.uni' % (tf))
		#vel.save(simPath + 'vel_low_%04d.uni' % (tf)) # optional
		xl_density.save(simPath + 'density_high_%04d.uni' % (tf))
		if(saveppm):
			projectPpmFull( xl_density, simPath + 'density_high_%04d.ppm' % (tf), 0, 5.0 )
			projectPpmFull( density, simPath + 'density_low_%04d.ppm' % (tf), 0, 5.0 )
	if not savedata and t>=timeOffset:
		if(saveppm):
			tf = t-timeOffset; simPath = "./"
			projectPpmFull( xl_density, simPath + 'density_high_%04d.ppm' % (tf), 0, 5.0 )
			projectPpmFull( density, simPath + 'density_low_%04d.ppm' % (tf), 0, 5.0 )

	sm.step()
	xl.step()
	#gui.screenshot( 'outt2_%04d.jpg' % t ) 
	#timings.display() 
	t = t+1


