#
# Helpers for handling command line parameters and the like
# example: path = getParam("path", "path.uni")
#
import sys, os, shutil, json

# global for checking used params
paramUsed = []
# additionally store parameters and values pairs
paramDict = {}

# ======================================================================================================================
# read parameters

#! check for a specific parameter, note returns strings, no conversion; not case sensitive! all converted to lower case
def getParam(name, default):
	global paramUsed
	v = default
	while( len(paramUsed)<len(sys.argv) ):
		paramUsed.append(0)
	for iter in range(1, len(sys.argv)):
		#if(iter <  len(sys.argv)-1): print("Param %s , used %d, val %s" %( sys.argv[iter].lower(), paramUsed[iter] , sys.argv[iter+1]) ); # debug
		if(sys.argv[iter].lower() == name.lower()) and (iter+1<len(paramUsed)):
			paramUsed[iter] = paramUsed[iter+1] = 1
			v = sys.argv[iter+1]
	paramDict[name] = v
	return v

def checkUnusedParams():
	global paramUsed
	err = False
	for iter in range(1, len(sys.argv)):
		if(paramUsed[iter]==0):
			print("Error: param %d '%s' not used!" % (iter,sys.argv[iter]) )
			err = True
	if err:
		exit(1)


# write param data to json file
def writeParams(filename="params.json", data=None):
	if data is None:
		data = paramDict
	with open(filename, 'w') as f:
		json.dump(data, f, indent=4)

# read param data from file
def readParams(filename="params.json"):
	with open(filename, 'r') as f:
		data = json.load(f)
	return data

def paramsToString():
	s = ""
	for keys,values in paramDict.items():
		s = s + "\t{}: {}\n".format(keys,values)
	return s


# ======================================================================================================================
# others / directory handling


# search & create next output dir
def getNextGenericPath(dirPrefix, folder_no = 1, basePath="../data/"):
	test_path_addition = '%s_%04d/' % (dirPrefix, folder_no)
	while os.path.exists(basePath + test_path_addition):
		folder_no += 1
		test_path_addition = '%s_%04d/' % (dirPrefix, folder_no)
		test_folder_no = folder_no
	test_path = basePath + test_path_addition
	print("Using %s dir '%s'" % (dirPrefix, test_path) )
	os.makedirs(test_path)
	return (test_path, folder_no)

def getNextTestPath(folder_no = 1, basePath="../data/"):
	return getNextGenericPath("test", folder_no, basePath)

def getNextSimPath(folder_no = 1, basePath="../data/"):
	return getNextGenericPath("sim", folder_no, basePath)

# custom Logger to write Log to file
class Logger(object):
	def __init__(self, test_path):
		self.terminal = sys.stdout
		self.log = open(test_path + "logfile.log", "a")

	def write(self, message):
		self.terminal.write(message)
		self.log.write(message)

	def flush(self): 
		# to avoid errormsg, " AttributeError: 'Logger' object has no attribute 'flush' "
		pass

