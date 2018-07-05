#!/bin/bash
# 
# note, uses & deletes test indices 7xx (models), and 8xx (outputs)
# to really clean all tests, remove test_07?? and test_08?? dirs, plus tiles in sim input dirs (uncomment below)
#
# by default, expects training data in sim 3000 (res 48), and applies to 3001 (res 96)
#

BPATH=../data/
SAFEDEL='rm -r'
SCENE_TRAIN=3000  # assumed to be resolution 48
SCENE_NP=3001     # assumed to be resolution 96
SHOW_INPUTS=1
EPOS=1000

echo Using data path ${BPATH}

# --- numpy npz file tests ---

# optionally, clear all tiles
#find ${BPATH}/sim_${SCENE_TRAIN}/ -iname "tiles*" -exec rm -fr \{\} \;
#find ${BPATH}/sim_${SCENE_NP}/    -iname "tiles*" -exec rm -fr \{\} \;

# remove old output dir
echo
echo "************** Test 1 (npz) **************"
${SAFEDEL} ${BPATH}test_0700
${SAFEDEL} ${BPATH}test_0800

# train a model
python tf_train.py out 0 basePath ${BPATH}  useVelocities 0  trainingEpochs ${EPOS}   alwaysSave 1  testPathStartNo 700 fromSim ${SCENE_TRAIN}        simSizeLow 48

# and apply to small data set
python tf_train.py out 1  basePath ${BPATH}  useVelocities 0  testPathStartNo 800  fromSim ${SCENE_NP}  loadModelTest 700     simSizeLow 96  outInputs ${SHOW_INPUTS}


# same for a model using velocities
echo
echo "************** Test 2 (npz,vel) **************"
${SAFEDEL} ${BPATH}test_0710
${SAFEDEL} ${BPATH}test_0810

python tf_train.py out 0 basePath ${BPATH}  useVelocities 1  trainingEpochs ${EPOS}  alwaysSave 1  testPathStartNo 710 fromSim ${SCENE_TRAIN}         simSizeLow 48
python tf_train.py out 1 basePath ${BPATH}  useVelocities 1  testPathStartNo 810  fromSim ${SCENE_NP}  loadModelTest 710     simSizeLow 96   outInputs ${SHOW_INPUTS}

# --- keras version tests ---

echo
echo "************** Test 3 (keras) **************"
${SAFEDEL} ${BPATH}test_0720
${SAFEDEL} ${BPATH}test_0820

python tf_train_keras.py out 0 basePath ${BPATH}  useVelocities 0  trainingEpochs ${EPOS}  testPathStartNo 720 fromSim ${SCENE_TRAIN}         simSizeLow 48
python tf_train_keras.py out 1 basePath ${BPATH}  useVelocities 0  testPathStartNo 820  fromSim ${SCENE_NP}  loadModelTest 720     simSizeLow 96   outInputs ${SHOW_INPUTS}

