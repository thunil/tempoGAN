#!/bin/bash
# 
# Run liquid (MLFLIP) test
#

BPATH=../data/
MANTAEXE=../build/manta
SAFEDEL='rm -r'  # optionally replace with sth safer...

echo Using data path \'${BPATH}\', and manta \'${MANTAEXE}\'

echo
echo "************** MLFLIP Test **************"

# optionally, remove old output dir
#${SAFEDEL} ${BPATH}manta-flip
#${SAFEDEL} ${BPATH}mlflip-tf

${MANTAEXE} manta_flip.py

${MANTAEXE} manta_gendata.py

python tf_train.py ${BPATH}manta-flip/training_data/ 

${MANTAEXE} manta_mlflip.py 

