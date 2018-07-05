#
# Testing symmetry of advection
# checks 2D & 3D , for 4 and 6 directions along axes (ech scalar and mac grid advection)
#
import sys
from manta import *
from helperInclude import *


showGui = False
if showGui and (GUI):
	gui = Gui()
	gui.show()
	gui.pause()

# print max error for each step?
showErrs = False

# globals
dirsSymm = [ 0,2,   1,2,   1,0 ]
dirsVel  = [ vec3(0,2,0), vec3(0,-2,0),
			 vec3(2,0,0), vec3(-2,0,0),
			 vec3(0,0,2), vec3(0,0,-2) ]

steps = 5
outputScale = 1e05
errThresh   = 1e-05
if (DOUBLEPRECISION):
	outputScale = 1e14

# loop over both 2d and 3d
for dim in range(2,4):

	# solver params
	res = 34

	gs = vec3(res,res,res)
	if (dim==2):
		gs.z=1
	s = Solver(name='main', gridSize = gs, dim=dim)
	s.timestep = 1.0
	accuracy   = 1e-3

	# error calculation
	errR1 = s.create(RealGrid)
	errV1 = s.create(RealGrid)
	# for 3d, we have 2 directions to check...
	errR2 = s.create(RealGrid)
	errV2 = s.create(RealGrid)

	# other grids
	flags = s.create(FlagGrid)
	vel   = s.create(MACGrid)
	rhs   = s.create(RealGrid)
	phi   = s.create(LevelsetGrid)
	pressure = s.create(RealGrid)

	# init scalar grid to advect... doesnt matter too much
	drop  = s.create(Sphere, center=gs*vec3(0.5,0.5,0.5), radius=res*0.25)




	for symms in range(2*dim):
	#if 1:
		#symms = 3; # single run

		# init & reset fields
		flags.initDomain(boundaryWidth=0)
		vel.setConst( vec3(0,0,0) )
		errR1.setConst(0)
		errV1.setConst(0)
		pressure.setConst(0)
		# note - ui problem with on-the-fly created LS grids, dont do "phi = drop.computeLs"
		rhs.setConst(0)
		phi.setConst(1e10)
		phi.join( drop.computeLevelset() )

		# setup velocity field

		# increase geometry along z for 2d
		fvOffsetZ = 0.0
		if (dim==2):
			fvOffsetZ = 1.25

		# init 2 swirl field
		flags.fillGrid()
		vel.setConst(vec3(0,0,0))

		dir1   = dirsSymm[symms-(symms%2)+0]
		dir2   = dirsSymm[symms-(symms%2)+1]
		velDir = dirsVel[symms]

		fluidVel = s.create(Box, p0=gs*vec3(0.30,0.30,0.30-fvOffsetZ), p1=gs*vec3(0.70,0.70,0.70+fvOffsetZ))
		fluidVel.applyToGrid( grid=vel , value=velDir )
		solvePressure(flags=flags, vel=vel, pressure=pressure, cgMaxIterFac=99., cgAccuracy=accuracy, retRhs = rhs)
				
		# symmetrize (note pressure checks here are just for debugging...)
		checkSymmetry    (a=pressure, err=errR1, axis=dir1)
		checkSymmetryVec3(a=vel     , err=errV1, axis=dir1)

		maxErrR = outputScale * errR1.getMax()
		maxErrV = outputScale * errV1.getMax()
		#print( "Initial symmetry err %f , %f" % (maxErrR, maxErrV) )
		#errV1.printGrid(zSlice=0)
		#s.step()

		# re-symmetrize
		checkSymmetry    (a=pressure, symmetrize=True, axis=dir1)
		checkSymmetryVec3(a=vel     , symmetrize=True, axis=dir1)

		# double check initial symmetry
		checkSymmetry    (a=pressure, err=errR1, axis=dir1)
		checkSymmetryVec3(a=vel     , err=errV1, axis=dir1)

		maxErrR = outputScale * errR1.getMax()
		maxErrV = outputScale * errV1.getMax()
		if(showErrs):
			print( "Initial symmetry check %f , %f " % (maxErrR, maxErrV) )

		if(dim==3): 
			checkSymmetry    (a=pressure, symmetrize=True, axis=dir2)
			checkSymmetryVec3(a=vel     , symmetrize=True, axis=dir2)

			checkSymmetry    (a=pressure, err=errR2, axis=dir2)
			checkSymmetryVec3(a=vel     , err=errV2, axis=dir2)

			maxErrR = outputScale * errR2.getMax()
			maxErrV = outputScale * errV2.getMax()
			if(showErrs):
				print( "Initial symmetry check %f , %f " % (maxErrR, maxErrV) )

		#errV1.printGrid(zSlice=0)
		#errR1.printGrid(zSlice=0) 
		#phi.printGrid(zSlice=0)

		#flags.updateFromLevelset(phi)

		# vel field done ...


		# add obstacle?
		if 1:
			obsBox = s.create(Box, p0=gs*vec3(0.4,0.4,0.4-fvOffsetZ), p1=gs*vec3(0.6,0.6,0.6+fvOffsetZ)) 
			obsBox.applyToGrid(grid=flags, value=(FlagObstacle) )


		print( "Checking symmetry, dirs " +str(dir1)+","+str(dir2)+ "  velocity " + str(velDir) )

		# part1
		print( "Checking scalar advection " )
		for t in range(steps):

			# re-symmetrize		
			checkSymmetry(a=phi, symmetrize=True, axis=dir1)
			if(dim==3):
				checkSymmetry(a=phi, symmetrize=True, axis=dir2)
			phi.setBoundNeumann(0)

			advectSemiLagrange(flags=flags, vel=vel, grid=phi, order=2, clampMode=1)

			checkSymmetry(a=phi, err=errR1, axis=dir1)
			maxErrR = outputScale * errR1.getMax()
			if(showErrs):
				print( "Max symmetry err1 real: %f " % (maxErrR) )

			if(dim==3):
				checkSymmetry(a=phi, err=errR2, axis=dir2)
				maxErrR = outputScale * errR2.getMax()
				if(showErrs):
					print( "Max symmetry err2 real: %f " % (maxErrR) )

			#errR1.printGrid(zSlice=0, printIndex=True) 
			s.step()

		# check result
		doTestGrid( sys.argv[0], ("errr1-%d-%d"%(dim,symms)) , s, errR1 , threshold=errThresh , thresholdStrict=1e-13 )
		if(dim==3):
			doTestGrid( sys.argv[0], ("errr2-%d-%d"%(dim,symms)) , s, errR2 , threshold=errThresh , thresholdStrict=1e-13 )


		# part2 
		print( "Checking vec3 advection " )
		for t in range(steps):
			phi.setBoundNeumann(0)

			# re-symmetrize		
			checkSymmetryVec3(a=vel, symmetrize=True, axis=dir1)
			if(dim==3):
				checkSymmetryVec3(a=vel, symmetrize=True, axis=dir2)

			advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2, clampMode=1) 
			#vel.printGrid(zSlice=0, printIndex=True)

			checkSymmetryVec3(a=vel, err=errV1, axis=dir1) 
			maxErrV = outputScale * errV1.getMax()
			if(showErrs):
				print( "Max symmetry err1 vec3: %f " % (maxErrV) )
			if(dim==3):
				checkSymmetryVec3(a=vel, err=errV2, axis=dir2) 
				maxErrV = outputScale * errV2.getMax()
				if(showErrs):
					print( "Max symmetry err2 vec3: %f " % (maxErrV) )

			# debug version, check symmetry only for a single component of the MAC grid
			#checkSymmetryVec3(a=vel, err=errV1, axis=dir1, disable=6); maxErrV = outputScale * errV1.getMax(); print( "Max err run X  %f " % (maxErrV))
			#checkSymmetryVec3(a=vel, err=errV1, axis=dir1, disable=5); maxErrV = outputScale * errV1.getMax(); print( "Max err run Y  %f " % (maxErrV))
			#checkSymmetryVec3(a=vel, err=errV1, axis=dir1, disable=3); maxErrV = outputScale * errV1.getMax(); print( "Max err run Z  %f " % (maxErrV))
			#errV1.printGrid(zSlice=0, printIndex=True)

			s.step()

		# check result
		doTestGrid( sys.argv[0], ("errv1-%d-%d"%(dim,symms)) , s, errV1 , threshold=errThresh , thresholdStrict=1e-12 )
		if(dim==3):
			doTestGrid( sys.argv[0], ("errv2-%d-%d"%(dim,symms)) , s, errV2 , threshold=errThresh , thresholdStrict=1e-12 )

	# prevent UI crash
	if showGui:
		print( "\nGUI right now does not support switching solvers & dimensions... exiting (choose 2d/3d beforehand)\n" )
		exit(1)


