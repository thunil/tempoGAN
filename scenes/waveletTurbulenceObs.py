#
# Complexer smoke simulation with wavlet turbulence plus
# - obstacle handling
# - and uv coordinates
# 

from manta import *
import os, shutil, math, sys

# how much to upres the XL sim?
# set to zero to disable the second one completely
upres = 4
 
# overall wavelet noise strength
wltStrength = 0.3 

# how many grids of uv coordinates to use (more than 2 usually dont pay off here)
uvs = 1

# and how many octaves of wavelet turbulence? could be set manually, but ideally given by upres factor (round)
octaves = 0
if(upres>0):
	octaves = int( math.log(upres)/ math.log(2.0) + 0.5 )

# simulation resolution 
dim = 2
res = 80
gs = vec3(res,int(1.5*res),res)
if (dim==2): gs.z = 1 

# setup low-res sim
sm = Solver(name='main', gridSize = gs, dim=dim)
sm.timestep = 1.5
# to do larger timesteps, without adaptive time steps, we need to set the length of a frame (1 is default)
sm.frameLength = sm.timestep 
timings = Timings()

# note - world space velocity, convert to grid space later
velInflow = vec3(0.015, 0, 0)

# inflow noise field
noise = NoiseField( parent=sm, fixedSeed=265, loadFromFile=True)
noise.posScale = vec3(20) # note, this is normalized to the grid size...
noise.clamp = True
noise.clampNeg = 0
noise.clampPos = 2
noise.valScale = 1
noise.valOffset = 0.075
noise.timeAnim = 0.3

# helper objects: inflow region, and obstacle
source    = Cylinder( parent=sm, center=gs*vec3(0.3,0.2,0.5), radius=res*0.081, z=gs*vec3(0.081, 0, 0))
sourceVel = Cylinder( parent=sm, center=gs*vec3(0.3,0.2,0.5), radius=res*0.15 , z=gs*vec3(0.15 , 0, 0))
obs       = Sphere( parent=sm,   center=gs*vec3(0.5,0.5,0.5), radius=res*0.15)

# larger solver, recompute sizes...
if(upres>0):
	xl_gs = vec3(upres*gs.x,upres*gs.y,upres*gs.z)
	if (dim==2): xl_gs.z = 1  # 2D
	xl = Solver(name='larger', gridSize = xl_gs, dim=dim)
	xl.timestep    = sm.timestep 
	xl.frameLength = xl.timestep 

	xl_flags   = xl.create(FlagGrid)
	xl_vel     = xl.create(MACGrid)
	xl_density = xl.create(RealGrid)

	xl_flags.initDomain()
	xl_flags.fillGrid()

	xl_source = Cylinder( parent=xl, center=xl_gs*vec3(0.3,0.2,0.5), radius=xl_gs.x*0.081, z=xl_gs*vec3(0.081, 0, 0))
	xl_obs    = Sphere(   parent=xl, center=xl_gs*vec3(0.5,0.5,0.5), radius=xl_gs.x*0.15)
	xl_obs.applyToGrid(grid=xl_flags, value=FlagObstacle)

	xl_noise = NoiseField( parent=xl, fixedSeed=265, loadFromFile=True)
	xl_noise.posScale = noise.posScale
	xl_noise.clamp    = noise.clamp
	xl_noise.clampNeg = noise.clampNeg
	xl_noise.clampPos = noise.clampPos
	xl_noise.valScale = noise.valScale
	xl_noise.valOffset = noise.valOffset
	xl_noise.timeAnim  = noise.timeAnim * upres


# init lower res solver & grids
bWidth=1
flags = sm.create(FlagGrid)
flags.initDomain(boundaryWidth=bWidth)
flags.fillGrid()
setOpenBound(flags, bWidth,'yY', FlagOutflow|FlagEmpty) 
obs.applyToGrid(grid=flags, value=FlagObstacle)

# create the array of uv grids
uv = []
for i in range(uvs):
	uvGrid = sm.create(VecGrid)
	uv.append(uvGrid)
	resetUvGrid( uv[i] )

vel       = sm.create(MACGrid) 
density   = sm.create(RealGrid)
pressure  = sm.create(RealGrid)
energy    = sm.create(RealGrid)
tempFlag  = sm.create(FlagGrid)

# wavelet turbulence noise field
xl_wltnoise = NoiseField( parent=xl, loadFromFile=True)
# scale according to lowres sim , smaller numbers mean larger vortices
# note - this noise is parented to xl solver, thus will automatically rescale
xl_wltnoise.posScale = vec3( int(1.0*gs.x) ) * 0.5
xl_wltnoise.timeAnim = 0.1

# setup user interface
if (GUI):
	gui = Gui()
	sliderStr = gui.addControl(Slider, text='turb. strength', val=wltStrength, min=0, max=2)
	gui.show()
	#gui.pause()

#printBuildInfo() 

# main loop
for t in range(200):
	mantaMsg('\nFrame %i, simulation time %f' % (sm.frame, sm.timeTotal))

	if (GUI):
		wltStrength = sliderStr.get()
	
	advectSemiLagrange(flags=flags, vel=vel, grid=density,  order=2)	
	advectSemiLagrange(flags=flags, vel=vel, grid=vel,      order=2, openBounds=True, boundaryWidth=bWidth )

	for i in range(uvs):
		advectSemiLagrange(flags=flags, vel=vel, grid=uv[i], order=2) 
		# now we have to update the weights of the different uv channels
		# note: we have a timestep of 1.5 in this setup! so the value of 16.5 means reset every 11 steps
		updateUvWeight( resetTime=16.5 , index=i, numUvs=uvs, uv=uv[i] ); 
		# also note, we have to update the weight after the advection 
		# as it is stored at uv[i](0,0,0) , the advection overwrites this...
		
	applyInflow=False
	if (sm.timeTotal>=0 and sm.timeTotal<50.):
		densityInflow( flags=flags, density=density, noise=noise, shape=source, scale=1, sigma=0.5 )
		sourceVel.applyToGrid( grid=vel , value=(velInflow*float(res)) )
		applyInflow=True
	
	setWallBcs(flags=flags, vel=vel)	
	addBuoyancy(density=density, vel=vel, gravity=vec3(0,-1e-3,0), flags=flags)

	vorticityConfinement( vel=vel, flags=flags, strength=0.4 )
	
	solvePressure(flags=flags, vel=vel, pressure=pressure , cgMaxIterFac=1.0, cgAccuracy=0.01 )
	setWallBcs(flags=flags, vel=vel)
	
	# determine weighting
	computeEnergy(flags=flags, vel=vel, energy=energy)

	# mark outer obstacle region by extrapolating flags for 2 layers
	tempFlag.copyFrom(flags)
	extrapolateSimpleFlags( flags=flags, val=tempFlag, distance=2, flagFrom=FlagObstacle, flagTo=FlagFluid );
	# now extrapolate energy weights into obstacles to fix boundary layer
	extrapolateSimpleFlags( flags=tempFlag, val=energy, distance=6, flagFrom=FlagFluid, flagTo=FlagObstacle ); 
	computeWaveletCoeffs(energy)

	#density.save('densitySm_%04d.vol' % t)
	sm.step()
	
	# xl ...
	if(upres>0):
		interpolateMACGrid( source=vel, target=xl_vel )
		
		# add all necessary octaves
		sStr = 1.0 * wltStrength 
		sPos = 2.0 
		for o in range(octaves):
			# add wavelet noise for each grid of uv coordinates 
			#xl_vel.clear() # debug , show only noise eval 
			for i in range(uvs):
				uvWeight = getUvWeight(uv[i]) 
				applyNoiseVec3( flags=xl_flags, target=xl_vel, noise=xl_wltnoise, scale=sStr * uvWeight, scaleSpatial=sPos , 
					weight=energy, uv=uv[i] )
			#mantaMsg( "Octave "+str(o)+", ss="+str(sStr)+" sp="+str(sPos)+" uvs="+str(uvs) ) # debug output 

			# update octave parameters for next iteration
			sStr *= 0.06 # magic kolmogorov factor
			sPos *= 2.0 
		
		# now advect
		for substep in range(upres): 
			advectSemiLagrange(flags=xl_flags, vel=xl_vel, grid=xl_density, order=2)    
		
		# manually recreate inflow
		if (applyInflow):
			 densityInflow( flags=xl_flags, density=xl_density, noise=xl_noise, shape=xl_source, scale=1, sigma=0.5 )
		
		#xl_density.save('densityXl_%04d.vol' % t) 
		xl.step()    

	#timings.display()
	# small and xl grid update done
	#gui.screenshot( 'wltObs_%04d.png' % t );

