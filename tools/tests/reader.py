#
# Read in generated data, display & compare it
# Needs manual file paths and grid dimension settings for now...
# 
import sys
from manta import *
from helperInclude import *

# use numbered files? or just base name?
appendNumber = False
# input file prefixes
basename1 = "test01.py"
basename2 = "test02.py"
# sim range
startFrame = 1
endFrame   = 150

# id string for grids to load (real, vec3, particles)
# leave string blank to skip loading
nameScalar = "phi" 
nameVec3   = "" # "vel"
nameParts  = "" #  "parts"

# setup

# try to load first file
gs = vec3(0,0,0)
if 1:
	for t in range(startFrame, endFrame):
		# optionally, also search for other grid types?
		gs = tryToGetSize( basename1, nameScalar , t, appendNumber )
		if(gs.x!=0):
			break;
		gs = tryToGetSize( basename2, nameScalar , t, appendNumber )
		if(gs.x!=0):
			break;
		gs = tryToGetSize( basename1, nameVec3 , t, appendNumber )
		if(gs.x!=0):
			break;
		gs = tryToGetSize( basename2, nameVec3 , t, appendNumber )
		if(gs.x!=0):
			break;
		if(appendNumber==False):
			break;

# solver params
dim = 3
if (gs.z==1):
	dim=2
print("Using grid size " + str(gs) + " , dim "+str(dim) )

# print info about running build, and those used to create data files?
buildInfo=1

# solver setup
s = Solver(name='main', gridSize = gs, dim=dim)
flags    = s.create(FlagGrid)

realErr  = s.create(RealGrid)
real1    = s.create(RealGrid)
real2    = s.create(RealGrid)

#ls1      = s.create(LevelsetGrid)
#ls2      = s.create(LevelsetGrid)
#lsErr    = s.create(RealGrid)

macErr   = s.create(VecGrid)
mac1     = s.create(VecGrid)
mac2     = s.create(VecGrid)

parts1    = s.create(BasicParticleSystem) 
pDens1    = parts1.create(PdataReal) 
pTest1    = parts1.create(PdataReal) 

parts2    = s.create(BasicParticleSystem) 
pDens2    = parts2.create(PdataReal) 

flags.initDomain(boundaryWidth=0)

if 1 and (GUI):
	gui = Gui()
	gui.show()
	gui.pause()
    
if(buildInfo==1):
	printBuildInfo() # more detailed build info , about what's running

# to be initialized later on...
realErrMax = 0
macErrMax = 0
lsErrMax = 0
partErrMax = 0

#main loop
for t in range(startFrame, endFrame):

	if(nameScalar != ""):
		tryToLoad( real1, basename1, nameScalar , t, appendNumber, buildInfo )
		tryToLoad( real2, basename2, nameScalar , t, appendNumber, buildInfo )
		realErr.copyFrom(real1);
		realErr.sub(real2);
		realErrMax = gridMaxDiff(real1, real2)
	
		#realErr.print(zSlice=15) 
		print("Max difference in step " +str(t) + " = "+ str(realErrMax) )

	if(nameVec3 != ""):
		tryToLoad( mac1, basename1, nameVec3  , t, appendNumber, buildInfo )
		tryToLoad( mac2, basename2, nameVec3  , t, appendNumber, buildInfo )
		macErr.copyFrom(mac1);
		macErr.sub(mac2);
		macErrMax = gridMaxDiffVec3(mac1, mac2)
	
		#macErr.print(zSlice=15) 
		print("Max vec3 difference in step " +str(t) + " = "+ str(macErrMax) )

	# load particles
	if(nameParts != ""):
		tryToLoad( parts1 , basename1, "parts"  , t, appendNumber, buildInfo )
		tryToLoad( parts2 , basename2, "parts"  , t, appendNumber, buildInfo )

	s.step()
    
