#!/usr/bin/python
import os
import shutil
import sys
import re
doDebug = False

# helper function to write output file
def writeHeader( filename, content ):
	try:
		outfile = open(filename, "w")
		outfile.write( content )
		outfile.close()
		if(doDebug):
			print( "Wrote '" + filename +"' " )
	except IOError:
		print("Warning, unable to write file '"+filename+"' ")
		exit(1)

dummyContent = "\n// could not determine git version\n\n"

# params

if len(sys.argv)<2:
	print("Usage makeHgVersion.py <out-file> <optional: path-to-git-exe> ")
	print("Warning, the target file <out-file> will be overwritten! ")
	exit(1)

# target file
outname = sys.argv[1]

# path to git executable, try a few options
# note /opt/local/bin/xxx is a double entry, can be overwritten by command line arg
exenames = [ "--replace--", "--replace--", "/opt/local/bin/git", "/usr/local/bin/git" ]
# check default
exenames[1] = os.popen("which git").read() 
exenames[1] = exenames[1].rstrip('\n')
# optionally, make argument
if len(sys.argv)>2:
	exenames[0] = sys.argv[2]

exename = ""
for nameCheck in exenames:
	#print "exe entry '"+nameCheck+"' "  # debug
	if( os.path.isfile(nameCheck) ):
		exename = nameCheck

# write empty file if no exe found
if(exename == ""):
	writeHeader( outname, dummyContent )
	print("Warning, no exe found - writing dummy header")
	exit(0); # dont throw an error for make, we can still continue...

if(doDebug):
	print("Params: outname '"+outname+"' , exename '"+exename)

# read old contents
oldContent = ""
doWrite    = True
try:
	infile = open(outname, "r")
	oldContent = infile.read()
	infile.close()
	if(doDebug):
		print("\n Old file content '"+oldContent+"' end \n")
except IOError:
	if(doDebug):
		print("Old file not found...")


# get git version
#gitVersion = os.popen(exename+" id").read() 
# get gid id
gitVersion = os.popen(exename+" log -1 ").read() 
# remove newlines...
if len(gitVersion)>0:
	gitVersion = gitVersion.splitlines()[0]
	gitVersion = gitVersion.rstrip('\n')
else:
	# probably git error, no repo or so; init default
	writeHeader( outname, dummyContent )
	print("Warning, git info failed - writing dummy header")
	exit(0)

if(doDebug):
	print( "Got git info: '" + gitVersion +"' " )

# matches old?
newContent = "\n\n#define MANTA_GIT_VERSION \"" + gitVersion + "\" \n\n" 

if(newContent == oldContent):
	if(doDebug):
		print("MATCHES! No write")
	doWrite = False
else:
	if(doDebug):
		print("Old info different, writing")

# write temp file
if(doWrite):
	writeHeader( outname, newContent )
	print( "Updated repository info header , "+gitVersion )


