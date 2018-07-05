import os
from subprocess import call

scenefile = "manta_genSimData.py" 
mantaexe = "../../build/manta"      # windows: "../../build/Release/manta"
my_env = os.environ.copy()
my_env["MANTA_DISABLE_UI"]="1"

for i in range(0,10):
	call([mantaexe, scenefile, "npSeed","%d"%i, "basePath","../data/"], env=my_env)

