#
# Helper functions independent of mantaflow 
# (exception: getFloatSetting, see below)
# 

import os
import re
import shutil


# ------------------------------------------------------------------------------------------
# smaller helpers (filenames etc.)

def outputFilename( file, gridname ):
	return file +"_"+ gridname + "_out.uni" 

# original, simpler...
def referenceFilename_old( file, gridname ):
	return file +"_"+ gridname + "_ref.uni" 


def getGenRefFileSetting( ):
	# check env var for generate data setting
	ret = int(os.getenv('MANTA_GEN_TEST_DATA', 0))
	# print("Gen-data-setting: " + str(ret))
	if(ret>0):
		return 1
	return 0

def getStrictSetting( ):
	print("Warning - deprecated, do not use! Strict thresholds are used automatically for double precision versions. ")
	# check env var whether strict mode enabled
	ret = int(os.getenv('MANTA_TEST_STRICT', 0))
	#print("Strict-test-setting: " + str(ret))
	if(ret>0):
		return 1
	return 0

# visual mode on? returns multiplier
def getVisualSetting( ):
	ret = int(os.getenv('MANTA_VISUAL', 0))
	if(ret>0):
		return ret
	return 0

# new version, extract directory & basename...
def referenceFilename( file, gridname ):
	(name,ext) = os.path.splitext( os.path.basename(file) )
	ddir = dataDirectory(file)
	suffix = "uni" 
	# double prec mode uses raw files (uni is always single prec float!)
	# uni files can be used to test IO , but strict threshold will cause "FAILs" then
	if getFloatSetting()==2: 
		suffix = "raw"
		#suffix = "uni" 
	return ddir+"/"+ name +"_"+ gridname + "." + suffix

def dataDirectory( file ):
	# extract path from script call
	basename = os.path.basename(file)
	basedir  = os.path.dirname (file)
	if len(basedir)==0:
		basedir = "."
	# different data sets for single/double
	dataDir = "testdata"
	if getFloatSetting()==2:
		dataDir = "testdataDouble"
	return basedir +"/"+ "../" + dataDir


# ------------------------------------------------------------------------------------------
# floating point accuracy , used throughout


# floating point accuracy setting - using settings in executable
# note, this "can" use manta python functions (if extCall is disabled)
globalFpAcc = -1
def getFloatSetting(extCall = False, platform = -1, mantaExe = ""):
	global globalFpAcc
	if globalFpAcc>0:
		return globalFpAcc

	# check build info string
	reCnt = 1
	if not extCall:
		from manta import printBuildInfo 
		buildInfo = printBuildInfo()
	else:
		# note - this call actually prints the info twice, once from the manta call, then from the buildInfo call...
		# hence, check for count of 2 here
		bifile = "helperBuildInfo.py"
		if   platform == 0: # unix
			buildInfo = os.popen(mantaExe + " " + bifile).read() 
		else: # windows
			buildInfo = os.popen('"'+ mantaExe + '" '+ bifile).read() 
		reCnt = 2
	
	fp1 = re.findall(r" fp1 ", buildInfo)
	fp2 = re.findall(r" fp2 ", buildInfo)
	#print( "FP '"+buildInfo+"' strings: " + format(fp1) + " " + format(fp2) + " | len: " + str(len(fp1)) + " " + str(len(fp2)) ) # debug
	if len(fp1) >= reCnt:
		globalFpAcc = 1
	elif len(fp2) >= reCnt:
		print("Double precision build detected")
		globalFpAcc = 2
	else:
		print("Error: Unable to determine floating point accuracy from running executable; Output: \n"+buildInfo +"\n")
		exit(1);
	return globalFpAcc


