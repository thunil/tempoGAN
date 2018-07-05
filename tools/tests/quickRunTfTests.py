import os

# assume it's called from tools/tests
os.chdir("../../tensorflow")

pyExe = "python"
mantaExe = "/Users/sinithue/devel/manta/buildMaster/manta"
dataPath = "./data"

runCommands = True
def doExe(cmd):
    if runCommands:
        os.system(cmd)
    else:
        print(cmd)

if os.path.exists(dataPath):
    print("Error - %s dir exists, stopping"%dataPath)
    exit(1)
else:
    if(runCommands):
	    os.makedirs(dataPath)

# ---
if 1: # ex0
    os.chdir("example0_simple")
    # gen 3 data sets
    doExe("{} manta_genSimSimple.py".format(mantaExe))
    doExe("{} manta_genSimSimple.py".format(mantaExe))
    doExe("{} manta_genSimSimple.py".format(mantaExe))
    doExe("{} tf_simple.py".format(pyExe))
    # TODO move output
    os.chdir("..")
    #exit(1)

# ---
if 1: # ex1
    os.chdir("example1_smoke_tiled")
    # gen 2 data sets
    doExe("{} manta_genSimData.py simNo 3000 res 48 ".format(mantaExe))
    doExe("{} manta_genSimData.py simNo 3000 res 96 steps 20 ".format(mantaExe))
    doExe("./runTests.sh")
    os.chdir("..")
    #exit(1)

if 1: # ex2
    os.chdir("example2_liquid")
    doExe("{} manta_flip.py ".format(mantaExe))
    doExe("{} manta_gendata.py ".format(mantaExe))
    doExe("{} tf_train.py ../data/manta-flip/training_data/ -s 1000 ".format(pyExe))
    doExe("{} manta_mlflip.py --t_end 0.5 ".format(mantaExe))  # note - will look funny, just for testing...
    os.chdir("..")
    #exit(1)

# ---
if 1: # check files
    os.chdir("refData")
    fns = os.popen("find .").read().split('\n')
    os.chdir("../data")
    missCnt = 0
    foundCnt = 0
    for f in fns:
        if len(f)<1: continue
        if f.find("events.out.tfevents.")!=-1: continue # skip TF event files

        #print("checking {}".format(f))
        if not os.path.exists(f):
            print("    Missing file {}".format(f))
            missCnt += 1
        else:
            foundCnt += 1
    if missCnt==0:
        print("All {} files found, nice!".format(foundCnt))
    else:
        print("Error - not all files generated, see above")


