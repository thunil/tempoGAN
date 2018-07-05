#******************************************************************************
#
# MantaFlow fluid solver framework
# Copyright 2017 Nils Thuerey
#
# This program is free software, distributed under the terms of the
# Apache License, Version 2.0 
# http://www.apache.org/licenses/LICENSE-2.0
#
# As-simple-as-possible manta- & tensor-flow example
# Hard coded for 2d 64*64 density data
#
#******************************************************************************

import time
import os
import shutil
import sys
import math
import random

import tensorflow as tf
import numpy as np
import scipy.misc
np.random.seed(13)
tf.set_random_seed(13)

sys.path.append("../tools")
import uniio

# path to sim data, trained models and output are also saved here
basePath = '../data/'

trainingEpochs = 2500
batchSize      = 10
inSize         = 64 * 64 * 1 # warning - hard coded to scalar values 64^2

# load data
densities = []

# start reading simSimple 1000 ff.
for sim in range(1000,2000): 
	if os.path.exists( "%s/simSimple_%04d" % (basePath, sim) ):
		for i in range(0,100): 
			filename = "%s/simSimple_%04d/density_%04d.uni" 
			uniPath = filename % (basePath, sim, i)  # 100 files per sim
			header, content = uniio.readUni(uniPath) # returns [Z,Y,X,C] np array
			h = header['dimX']
			w  = header['dimY']
			arr = content[:, ::-1, :, :] # reverse order of Y axis
			arr = np.reshape(arr, [w, h, 1]) # discard Z
			densities.append( arr )

loadNum = len(densities)
if loadNum<200:
	print("Error - use at least two full sims, generate data by running 'manta ./manta_genSimSimple.py' a few times..."); exit(1)

densities = np.reshape( densities, (len(densities), 64,64,1) )

print("Read uni files, total data " + format(densities.shape) )
valiSize = max(100, int(loadNum * 0.1)) # at least 1 full sim...
valiData = densities[loadNum-valiSize:loadNum,:] 
densities = densities[0:loadNum-valiSize,:]
print("Split into %d training and %d validation samples" % (densities.shape[0], valiData.shape[0]) )
loadNum = densities.shape[0]



# set up the network

x = tf.placeholder(tf.float32, shape=[None, 64,64, 1])
y = tf.placeholder(tf.float32, shape=[None, 64,64, 1])

xIn = tf.reshape(x, shape=[-1, inSize ]) # flatten
fc_1w = tf.Variable(tf.random_normal([inSize, 50], stddev=0.01))
fc_1b   = tf.Variable(tf.random_normal([50], stddev=0.01))

fc1 = tf.add(tf.matmul(xIn, fc_1w), fc_1b)
fc1 = tf.nn.tanh(fc1)
fc1 = tf.nn.dropout(fc1, 0.5) # plenty of dropout...

fc_2w = tf.Variable(tf.random_normal([50, inSize], stddev=0.01))  # back to input size
fc_2b = tf.Variable(tf.random_normal([inSize], stddev=0.01))

y_pred = tf.add(tf.matmul(fc1, fc_2w), fc_2b)
y_pred = tf.reshape( y_pred, shape=[-1, 64, 64, 1])

cost = tf.nn.l2_loss(y - y_pred) 
opt  = tf.train.AdamOptimizer(0.0001).minimize(cost)



# now we can start training...

print("Starting training...")
sess = tf.InteractiveSession()
sess.run(tf.global_variables_initializer())

for epoch in range(trainingEpochs):
	c = (epoch * batchSize) % densities.shape[0]
	batch = []
	for currNo in range(0, batchSize):
		r = random.randint(0, loadNum-1) 
		batch.append( densities[r] )

	_ , currentCost = sess.run([opt, cost], feed_dict={x: batch, y: batch})
	#print("Epoch %d/%d: cost %f " % (epoch, trainingEpochs, currentCost) ) # debug, always output cost
	
	if epoch%10==9 or epoch==trainingEpochs-1:
		[valiCost,vout] = sess.run([cost, y_pred], feed_dict={x: valiData, y: valiData})
		print("Epoch %d/%d: cost %f , validation cost %f " % (epoch, trainingEpochs, currentCost, valiCost) )

		if epoch==trainingEpochs-1:
			outDir = "%s/test_simple" % basePath
			if not os.path.exists(outDir): os.makedirs(outDir)
			print("\n Training done. Writing %d images from validation data to directory %s..." % (len(valiData),outDir) )
			for i in range(len(valiData)):
				scipy.misc.toimage( np.reshape(valiData[i], [64, 64]) , cmin=0.0, cmax=1.0).save("%s/in_%d.png" % (outDir,i))
				scipy.misc.toimage( np.reshape(vout[i]    , [64, 64]) , cmin=0.0, cmax=1.0).save("%s/out_%d.png" % (outDir,i))



print("Done")

