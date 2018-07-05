#
# Helper functions for test runs in mantaflow
# 

from manta import *
import os
import shutil
import re
from helperGeneric import *


# ------------------------------------------------------------------------------------------
# test result checking


def checkResult( name, result, resultRel , thresh, threshStrict, invertResult=False ):
	curr_thresh = thresh

	# enable strict thresholds for double prec tests
	if(getFloatSetting()==2):
		curr_thresh = threshStrict
	print ("Checking '%s', result=%f , thresh=%f" % ( name , result , curr_thresh) )

	if   ( ( result > 0.) and (result < 1e-04) ):
		print ("Note: small difference: %f (output scaled by 1e5)" % ( result * 1e05 ) ) # debugging...
	elif ( ( result > 0.) and (result < 1e-08) ):
		print ("Note: small difference: %f (output scaled by 1e9)" % ( result * 1e09 ) ) # debugging...
	#elif ( result == 0.0):
		#print ("Result is really zero...")

	allGood = 0
	if ( result <= curr_thresh) :
		allGood = 1

	# for checks that should fail
	if ( invertResult == True) :
		if ( allGood == 0) :
			allGood = 1
		else:
			allGood = 0

	# now react on outcome...
	if ( allGood == 1 ):
		print("OK! Results for "+name+" match...")
		return 0
	else:
		print("FAIL! Allowed "+name+" threshold "+str(curr_thresh)+", results differ by "+str(result) +" (abs) , and by "+str(resultRel)+" (rel)" )
		return 1



# global var to print manta version once per test
printVersion = 1

# compare a grid, in generation mode (MANTA_GEN_TEST_DATA=1) it
# creates the data on disk, otherwise it loads the disk data,
# computes the largest per cell error, and checks whether it matches
# the allowed thresholds
#
# note, there are two thresholds:
# 	- the "normal" one is intended for comparing single precision calculations across different compilers
#	- the "strict" one for double precision compiles (detected automatically)
#   - the "grid" object can be either a Grid<T>, or a ParticleDataImpl<T> ; parent is either FluidSolver or ParticleSystem
#
def doTestGrid( file , name, parent , grid, threshold=0, thresholdStrict=0, invertResult=False, debugShowDifference=False ):
	global printVersion

	# both always have to given together (if not default)
	if ( threshold!=0 and thresholdStrict==0 ):
		print( "Error doTestGrid - give both thresholds at the same time...")
		return 1
	if ( threshold==0 and thresholdStrict!=0 ):
		print( "Error doTestGrid - give both thresholds at the same time...")
		return 1
	if getVisualSetting():
		print( "Visual mode, skipping data file checks" )
		return 0

	# handle grid types that need conversion
	#print( "doTestGrid, incoming grid type :" + type(grid).__name__ + " class:"+grid._class+ " T:"+grid._T )
	if ( type(grid).__name__ == "MACGrid" ):
		gridTmpMac = parent.create(VecGrid)
		copyMacToVec3(grid , gridTmpMac )
		ret = doTestGrid( file, name, parent, gridTmpMac , threshold, thresholdStrict, invertResult, debugShowDifference)
		if debugShowDifference: grid.copyFrom( gridTmpMac )
		return ret
	if ( type(grid).__name__ == "LevelsetGrid" ):
		gridTmpLs = parent.create(RealGrid)
		copyLevelsetToReal(grid , gridTmpLs )
		ret = doTestGrid( file, name, parent, gridTmpLs  , threshold, thresholdStrict, invertResult, debugShowDifference)
		if debugShowDifference: grid.copyFrom( gridTmpLs )
		return ret

	# now we should only have real & vec3 grids

	# sanity check data type & parent
	if ( grid._class == "Grid" and parent._class != "FluidSolver" ):
		print( "Error doTestGrid - pass fluid solver as parent for grids, is '"+ parent._class +"'" );
		return 1
	if ( grid._class == "ParticleDataImpl" and parent._class != "BasicParticleSystem" ):
		print( "Error doTestGrid - pass particle system as parent for pdata" );
		return 1

	# create temp grid (parent can be either fluid solver or particle system)
	if ( grid._class == "Grid" and grid._T == "Real" ):
		compareTmpGrid = parent.create(RealGrid)
	elif ( grid._class == "Grid" and grid._T == "Vec3" ):
		compareTmpGrid = parent.create(VecGrid)
	elif ( grid._class == "Grid" and grid._T == "int" ):
		compareTmpGrid = parent.create(IntGrid)
	elif ( grid._class == "ParticleDataImpl" and grid._T == "Real" ):
		compareTmpGrid = parent.create(PdataReal)
	elif ( grid._class == "ParticleDataImpl" and grid._T == "Vec3" ):
		compareTmpGrid = parent.create(PdataVec3)
	elif ( grid._class == "ParticleDataImpl" and grid._T == "int" ):
		compareTmpGrid = parent.create(PdataInt)
	elif ( grid._class == "Grid4d" and grid._T == "Real" ):
		compareTmpGrid = parent.create(Grid4Real)
	elif ( grid._class == "Grid4d" and grid._T == "int" ):
		compareTmpGrid = parent.create(Grid4Int)
	elif ( grid._class == "Grid4d" and grid._T == "Vec3" ):
		compareTmpGrid = parent.create(Grid4Vec3)
	elif ( grid._class == "Grid4d" and grid._T == "Vec4" ):
		compareTmpGrid = parent.create(Grid4Vec4)
	else:
		print( "Error doTestGrid - unknown grid type " + type(grid).__name__+ " class:"+grid._class+ " T:"+grid._T  )
		return 1

	genRefFiles = getGenRefFileSetting()
	fname = referenceFilename( file, name )

	if (genRefFiles==1):
		grid.save( fname )
		print( "OK! Generated reference file '" + fname + "'")

		# test data generation log
		if 1:
			infofilename = dataDirectory(file)+"/test_data_info.txt"
			text_file = open(dataDirectory(file)+"/test_data_info.txt", "a");
			if printVersion:
				printVersion = 0
				text_file.write( "\n%s, %s\n" % (file, str(printBuildInfo())) );
			text_file.write( "    %s\n" % ( fname ) );
			text_file.close();
		return 0
	else:
		# give error if file doesnt exist
		if( not os.path.isfile( fname ) ):
			print( "Error - unable to load test file %s" % referenceFilename( file, name ) )
			print("FAIL! Reference data missing..." );
			return 1

		compareTmpGrid.load( fname )

		errVal = 1e10
		if ( grid._class == "Grid" and grid._T == "Real" ):
			errVal = gridMaxDiff    ( grid, compareTmpGrid )
		elif ( grid._class == "Grid" and grid._T == "Vec3" ):
			errVal = gridMaxDiffVec3( grid, compareTmpGrid )
		elif ( grid._class == "Grid" and grid._T == "int" ):
			errVal = gridMaxDiffInt ( grid, compareTmpGrid )
		elif ( grid._class == "ParticleDataImpl" ):
			errVal = pdataMaxDiff ( grid, compareTmpGrid )
		elif ( grid._class == "Grid4d" and grid._T == "Real" ):
			errVal = grid4dMaxDiff    ( grid, compareTmpGrid )
		elif ( grid._class == "Grid4d" and grid._T == "int" ):
			errVal = grid4dMaxDiffInt ( grid, compareTmpGrid )
		elif ( grid._class == "Grid4d" and grid._T == "Vec3" ):
			errVal = grid4dMaxDiffVec3( grid, compareTmpGrid )
		elif ( grid._class == "Grid4d" and grid._T == "Vec4" ):
			errVal = grid4dMaxDiffVec4( grid, compareTmpGrid )
		else:
			print( "Error doTestGrid - error calculation missing" )
			return 1

		# debug mode to return difference in source grid, warning - no error measurements possible anymore
		if debugShowDifference: 
			print("Warning debugShowDifference active, test data invalidated for UI display")
			grid.sub( compareTmpGrid )
			return 0

		# debug info , print min/max
		if 0:
			minVal1 = grid.getMin()
			maxVal1 = grid.getMax()
			minVal2 = compareTmpGrid.getMin()
			maxVal2 = compareTmpGrid.getMax()
			print( "Test "+name+" min/max curr "+str(minVal1)+" to "+str(maxVal1)+" min/max ref "+str(minVal2)+" to "+str(maxVal2) );

		maxVal = grid.getMaxAbs() + 1e-15
		errValRel = errVal/maxVal

		# finally, compare max error to allowed threshold, and return result
		return checkResult( name, errVal , errValRel, threshold , thresholdStrict, invertResult )


# ------------------------------------------------------------------------------------------
# smaller helpers (directories, global settings)


# for xl test, load test data afterwards to keep sims in sync
def doTestDataLoad( file , name, solver , grid ):
	genRefFiles = getGenRefFileSetting()

	if (genRefFiles!=1):
		print( "Loading %s" % referenceFilename( file, name ) )
		grid.load( referenceFilename( file, name ) )

# reset and generate info file with version string when in data gen mode
def doResetInfoFile( file ):
	if(getGenRefFileSetting()==1):
		infofilename = dataDirectory(file)+"/test_data_info.txt"
		print( "Resetting test data info file "+infofilename )
		text_file = open(dataDirectory(file)+"/test_data_info.txt", "w");
		text_file.write( "\n%s\n\n" % (str(printBuildInfo())) );
		text_file.close();

# read test data

# try to load uni file if it exists
def tryToGetSize( basename, suffix, number , appendNumber ):
	if(appendNumber==True):
		suffix = suffix+("_%04d" % number )
	rfile = referenceFilename( basename, suffix ) 
	#print("Trying to get grid size from " + rfile)
	size = vec3(0,0,0)
	if(os.path.isfile(rfile)):
		size = getUniFileSize(rfile) 
		#print("Found " + str(size) )
	return size

# configure input filenames

# try to load uni file if it exists
def tryToLoad( grid, basename, suffix, number , appendNumber , buildInfo ):
	if(appendNumber==True):
		suffix = suffix+("_%04d" % number )
	rfile = referenceFilename( basename, suffix ) 
	print("Trying to load " + rfile)
	if(os.path.isfile(rfile)):
		grid.load(rfile)
		if(buildInfo==1):
			printUniFileInfoString(rfile) # more detailed build info
	else:
		grid.clear()
	return 1


