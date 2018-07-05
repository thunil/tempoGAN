#!/usr/bin/python
#
# Simple script to run all base tests
# The test scripts are named test_XXXX_description.py
# They are assumed to be in ascending order of complexity.
#
# Note - the visual mode (MANTA_VISUAL=1) requires ImageMagick and gnuplot!
#
# Rough ordering:
# 0xxx tests are very basic (mostly single operator calls)
# 1xxx is for 2d sims
# 2xxx for "real" 3d sims
# 204x ff. are liquids
# 30xx ff. are for non-fluid related things 
# 
import datetime
import platform
import os
import shutil
import sys
import re
from subprocess import *
from helperGeneric import *

# note - todo, right now this script assumes the test_X.py files are
# in the current working directory... this should be changed at some 
# point.

# determine platform
platform=platform.system()
if   platform == 'Linux' or platform == 'Darwin':
	platform = 0; # unix
elif platform == 'Windows' and os.name=='nt':
	# windows specifics:
	#    - some commands depend on the operating system - check platform
	#    - cygwin under windows doesn't behave like linux nor like windows at my computer and cannot be distinguished
	#      from the normal windows console -> set platform manually to 2 if needed
	#      cygwin should be detected by os.name==posix but it doesn't work at my computer
	#    - some print messages needed to set default codepage of the console to UTF-8 with the command "chcp 65001"
	platform = 1; # windows
elif platform == 'Windows' and os.name == 'posix':
	platform = 2; # probably cygwin
else:
	platform = 0; 
	print("OS cannot be detected (%s), assuming unix" % platform);
 

# debugging, print outputs from all manta calls
printAllOutpus = 0
filePrefix = "test_"

# by default, disable UI (if compiled)
os.environ["MANTA_DISABLE_UI"] = "1"

if(len(sys.argv)<2):
	print("Usage runTests.py <manta-executable>")
	exit(1)

manta = sys.argv[1]
print ("Using mantaflow executable '" + manta + "' " )


# extract path from script call
basedir  = os.path.dirname (sys.argv[0])
#os.path.splitext(base)
if len(basedir)==0:
	basedir = "."

# store test data in separate directory
getFloatSetting(extCall = True, platform=platform, mantaExe=manta ) # check from external call
datadir = dataDirectory(sys.argv[0])


#unix only: currdate = os.popen("date \"+%y%m%d%H%M\"").read() 
if getGenRefFileSetting():
	print ("\nNote - generating test data for all tests!")
	print ("Tests results will not be evaluated...\n")

	# create directory if necessary
	if not os.path.exists( datadir ):
		os.makedirs( datadir )    

currdate = '{:%Y%m%d_%H%M%S}'.format(datetime.datetime.now())
currdate = str(currdate)[:-1]

# in visual mode, also track runtimes
visModeTrashDir = basedir+"/trash"
outpngdir       = basedir+"/result_%s"%(currdate)
if getVisualSetting():
	dirname = basedir+"/runtimes"
	if not os.path.exists( dirname ):
		os.makedirs( dirname )    
	# make sure no previous files are left
	if not os.path.exists( visModeTrashDir ):
		os.makedirs( visModeTrashDir )    
	os.popen( "mv -f ./test_*.ppm %s"%(visModeTrashDir) )

	if not os.path.exists( outpngdir ):
		os.makedirs( outpngdir )	
	print("Note - running in visual test mode...");

# limit the runs for debugging
visModeDebugCount = 0

# old: files = os.popen("ls "+basedir+"/"+str(filePrefix)+"????_*.py").read() 
if   platform==0:
	files = os.popen("ls "+basedir+"/"+str(filePrefix)+"????_*.py").read() 
elif platform==1:
	files = os.popen("dir /a-d /b "+basedir+"\\"+str(filePrefix)+"????_*.py").read() 
elif platform==2:
	files = os.popen("ls "+basedir+str(filePrefix)+"????_*.py").read() # some cygwin under windows
#print ("Debug - using test scene files: "+files)

# for graphs in visual mode only
gnuplotExe = "/usr/bin/gnuplot"
if len( os.getenv('MANTA_GNUPLOT', "") )>0:
	gnuplotExe = os.getenv('MANTA_GNUPLOT', "")
	print("Using gnuplot at %s" % gnuplotExe);

# ready to go...
currdate = os.popen("date \"+%y%m%d%H%M\"").read() 
currdate = str(currdate)[:-1]

# in visual mode, also track runtimes
visModeTrashDir = basedir+"/trash"
outpngdir       = basedir+"/result_%s"%(currdate)
if getVisualSetting():
	dirname = basedir+"/runtimes"
	if not os.path.exists( dirname ):
		os.makedirs( dirname )	
	# make sure no previous files are left
	if not os.path.exists( visModeTrashDir ):
		os.makedirs( visModeTrashDir )	
	os.popen( "mv -f ./test_*.ppm %s"%(visModeTrashDir) )

	if not os.path.exists( outpngdir ):
		os.makedirs( outpngdir )	
	print("Note - running in visual test mode...");

# limit the runs for debugging
visModeDebugCount = 0

files = os.popen("ls "+basedir+"/"+str(filePrefix)+"????_*.py").read() 
#print "Debug - using test scene files: "+files


# ready to go...

num = 0
numOks = 0
numFail = 0
failedTests = ""

files = files.split('\n')
for file in files:
	if ( len(file) < 1) or (visModeDebugCount>1):
		continue

	(utime1, stime1, cutime1, cstime1, elapsed_time1) = os.times() 

	num += 1
	print("Running '" + file + "' ");
	# result = os.popen(manta + " " + file + " 2>&1 ").read() 
	if   platform == 0:
		result = os.popen(manta + " " + file + " 2>&1 ").read() 
	elif platform == 1:
		result = os.popen('"'+ manta + '" '+ file).read() 
	elif platform == 2:
		result = os.popen('"'+ manta + '" '+ file + " 2>&1 ").read() 

	(utime2, stime2, cutime2, cstime2, elapsed_time2) = os.times() 

	(utime2, stime2, cutime2, cstime2, elapsed_time2) = os.times() 

	oks = re.findall(r"\nOK!", result)
	#print oks
	numOks += len(oks)

	fails = re.findall(r"\nFAIL!", result)
	# also check for "Errors"
	if (len(fails)==0) :
		# note - if there are errors, the overall count of tests won't be valid anymore... 
		fails = re.findall(r"Error", result) 
	
	if (len(fails)!=0) :
		# we had failed tests! add to list...
		failedTests = failedTests+" "+file 
	#print fails
	numFail += len(fails)

	if (len(fails)>0) or (printAllOutpus==1):
		print(" ");
		print("Full output: " + result);
		print(" ");

	# store benchmarking results (if theres any output) , and generate plots
	timefile = "%s/runtimes/%s_v%d" % (basedir, os.path.basename(file), getVisualSetting()) 
	if getVisualSetting() and ( os.path.isfile( "%s_0001.ppm"%(file) ) or os.path.isfile(timefile+".time") ) and os.path.isfile( gnuplotExe ):
		runtime = elapsed_time2-elapsed_time1 
		if runtime>0.0:
			text_file = open(timefile+".time", "a");
			text_file.write( "%s %f \n" % (currdate,runtime) );
			text_file.close();
		else:
			print("Zero runtime! Something went wrong...");
		
		if os.path.isfile(gnuplotExe):
			plot = Popen(gnuplotExe, stdin=PIPE)
			plot.stdin.write("unset key\n")
			plot.stdin.write("set terminal png\n")
			plot.stdin.write("set output '%s.png'\n" % (timefile) )
			plot.stdin.write("set terminal png size 1000, 700\n")
			plot.stdin.write("set ylabel \"Time [s]\"\n")
			plot.stdin.write("set xlabel \"%s, Date\"\n" % (timefile) )
			plot.stdin.write("plot '%s.time' using 1:2 with lines\n" % (timefile) )
			plot.stdin.write("quit\n")

		# for debugging, only execute a few files
		#visModeDebugCount += 1

		# generate log file for visual run, same dir as output pngs
		logfilename = "%s/%s.log" % (outpngdir, os.path.basename(file) )
		log_file = open(logfilename, "w");
		log_file.write(result); log_file.close();


if getGenRefFileSetting():
	print("Test data generated");
	exit(0)

if getVisualSetting():
	print("Viusal data generated");
	# now convert & remove all ppms"
	os.popen("for i in ./test_*.ppm ; do convert $i $(basename $i .ppm).png; done")
	os.popen( "mv ./test_*.png %s"%(outpngdir)  )
	os.popen( "mv -f ./test_*.ppm %s"%(visModeTrashDir) )
	exit(0)

print("\n ============================================= \n");
print("Test summary: "  + str(numOks) + " passed, " + str(numFail) + " failed.   (from "+str(num) + " files) ");

if (numFail==0) and (numOks==0):
	print("Failure! Perhaps manta executable didnt work, or runTests not called in test directory \n");
	exit(1)
elif (numFail==0) and (numOks>0):
	print("All good :) \n");
	exit(0)
else:
	print("Oh no :( the following tests failed: %s \n" % failedTests);
	exit(2)

