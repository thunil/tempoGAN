#******************************************************************************
#
# MantaFlow fluid solver framework
# Copyright 2017 Daniel Hook, Nils Thuerey
#
# This program is free software, distributed under the terms of the
# Apache License, Version 2.0 
# http://www.apache.org/licenses/LICENSE-2.0
#
# Manta & tensor flow example with tiles
#
#******************************************************************************

import time
import os
import shutil
import sys
import math

import tensorflow as tf
import numpy as np

# load manta tools
sys.path.append("../tools")
import tilecreator as tiCr
import uniio
import paramhelpers as ph

# path to sim data, trained models and output are also saved here
basePath = '../data/'

# main mode switch:
outputOnly   = True  # apply model, or run full training?
outputInputs = False

simSizeLow   = 64
tileSizeLow  = 16
upRes        = 4

# dont use for training! only for applying model, add overlap here if necessary (i.e., cropOverlap>0) 
# note:  cropTileSizeLow + (cropOverlap * 2) = tileSizeLow
cropOverlap     = 0
cropTileSizeLow = tileSizeLow - 2*cropOverlap

emptyTileValue  = 0.01
learningRate    = 0.00005
trainingEpochs  = 10000 # for large values, stop manualy with ctrl-c...
dropout         = 0.1   # slight
batchSize       = 96
testInterval    = 200
saveInterval    = 1000
fromSim = toSim = -1    # range of sim data directories to use
keepAll         = False
numTests        = 10    # evaluate on 10 data points from test data
randSeed        = 1
brightenOutput  = -1    # multiplied with output to brighten it up

# optional, add velocity as additional channels to input
useDensity      = 1     # default, only density
useVelocities   = 0

# load pressure as reference data (instead of density by default)
outputPressure  = 0
outputDataName  = '' # name of data to be regressed; by default, does nothing (density), e.g. if output data is pressure set to "pressure"
bWidth          = -1 # special: boundaryWidth to be cut away, in line with "bWidth" in manta scene files

# previous test to load (id X for test_X dir)
loadModelTest = -1
loadModelNo   = -1 # specific saved model to load, searches for latest by default
testPathStartNo = 1

# ---------------------------------------------

# command line params, explanations mostly above with variables
outputOnly      = int(ph.getParam( "out",             outputOnly ))>0
trainingEpochs  = int(ph.getParam( "trainingEpochs",  trainingEpochs ))
loadModelTest   = int(ph.getParam( "loadModelTest",   loadModelTest))  
loadModelNo     = int(ph.getParam( "loadModelNo",     loadModelNo))    
basePath        =     ph.getParam( "basePath",        basePath)
useDensity		= int(ph.getParam( "useDensity",      useDensity  ))
useVelocities   = int(ph.getParam( "useVelocities",   useVelocities  ))
outputPressure	= int(ph.getParam( "outputPressure",  outputPressure  ))
testPathStartNo = int(ph.getParam( "testPathStartNo", testPathStartNo  ))
fromSim         = int(ph.getParam( "fromSim",         fromSim  )) 
toSim           = int(ph.getParam( "toSim",           toSim  ))
alwaysSave      = int(ph.getParam( "alwaysSave",      False  )) # by default, only save checkpoint when cost is lower, can be turned off here
randSeed        = int(ph.getParam( "randSeed",        randSeed )) 
simSizeLow      = int(ph.getParam( "simSizeLow",      simSizeLow )) 
upRes           = int(ph.getParam( "upRes",           upRes ))
outputInputs    = int(ph.getParam( "outInputs",       outputInputs)) # for debugging, write images for input data
brightenOutput  = int(ph.getParam( "brightenOutput",  brightenOutput)) 
outputDataName  =    (ph.getParam( "outName",         outputDataName))
bWidth			= int(ph.getParam( "bWidth",          bWidth))
ph.checkUnusedParams()

# initialize
tiCr.setBasePath(basePath) 
np.random.seed(randSeed)
tf.set_random_seed(randSeed)

simSizeHigh  = simSizeLow   * upRes
tileSizeHigh = tileSizeLow  * upRes

if toSim==-1:
	toSim = fromSim

# debug helper, copy sim data to different ID
#tiCr.copySimData( fromSim, toSim ); exit(1);  # uncomment to run...

if not outputOnly:
	# run training!
	loadModelTest = -1
	if fromSim==-1:
		fromSim = toSim = 1000 # short, use single sim

	if cropOverlap>0:
		print("ERROR: dont use cropOverlap != 0 for training...")
		exit(1)

else:
	keepAll = True
	emptyTileValue = -1. # dont discard any
	# dont train, just apply to input seq, by default use plume (2004)
	if fromSim==-1:
		fromSim = toSim = 3000

# check if batchsize is multiple of tilesInImage
tileSizeHiCrop = upRes * cropTileSizeLow
tilesPerImg = (simSizeHigh // tileSizeHiCrop) ** 2

# ---------------------------------------------

n_input  = tileSizeLow  ** 2 
n_output = tileSizeHigh ** 2
n_inputChannels = 0

if useDensity:
	n_inputChannels += 1
if useVelocities:
	n_inputChannels += 3
if n_inputChannels == 0:
	print("ERROR: No inputs set. Use useDensity, useVelocities, ...")
	exit()
n_input *= n_inputChannels

# create model loading path
if not loadModelTest == -1:
	if not os.path.exists(basePath + 'test_%04d/' % loadModelTest):
		print('ERROR: Test to load does not exist.')
	# search for newest model if no loadModelNo is given
	if loadModelNo == -1:
		for currModel in range(0, 200):
			if os.path.isfile(basePath + 'test_%04d/model_%04d.ckpt.index' % (loadModelTest, currModel)):
				loadModelNo = currModel
		if loadModelNo == -1:
			print('ERROR: Model with id below 200 does not exist. Please specify model id as "loadModelNo".')
			exit()
		# print('Latest model: %d.' % loadModelNo)

	load_path = basePath + 'test_%04d/model_%04d.ckpt' % (loadModelTest, loadModelNo)

(test_path,test_folder_no) = ph.getNextTestPath(testPathStartNo, basePath)
if not outputOnly: uniio.backupFile(__file__, test_path)

sys.stdout = ph.Logger(test_path)

# print Variables to log
def print_variables():
	print('\nUsing variables:')
	print('fromSim: {}'.format(fromSim))
	print('toSim: {}'.format(toSim))
	print('simSizeLow: {}'.format(simSizeLow))
	print('tileSizeLow: {}'.format(tileSizeLow))
	print('cropOverlap: {}'.format(cropOverlap))
	print('cropTileSizeLow: {}'.format(cropTileSizeLow))
	print('upRes: {}'.format(upRes))
	print('emptyTileValue: {}'.format(emptyTileValue))
	print('learningRate: {}'.format(learningRate))
	print('trainingEpochs: {}'.format(trainingEpochs))
	print('dropout: {}'.format(dropout))
	print('batchSize: {}'.format(batchSize))
	print('\n')

print_variables()

# ---------------------------------------------
# TENSORFLOW SETUP
x      = tf.placeholder(tf.float32, shape=[None, 1, tileSizeLow, tileSizeLow, n_inputChannels])  # fixed to 2D for now
y_true = tf.placeholder(tf.float32, shape=[None, 1, tileSizeHigh,tileSizeHigh, 1             ])  # also fixed, always output 1 channel
training = tf.placeholder(tf.bool)

hl1_size = 512 
xIn = tf.reshape(x, shape=[-1, n_input ]) # flatten
fc_1_weight = tf.Variable(tf.random_normal([n_input, hl1_size], stddev=0.01))
fc_1_bias   = tf.Variable(tf.random_normal([hl1_size], stddev=0.01))

fc1 = tf.add(tf.matmul(xIn, fc_1_weight), fc_1_bias)
fc1 = tf.nn.tanh(fc1)
fc1 = tf.layers.dropout(fc1, rate=dropout, training=training)

fc_out_weight = tf.Variable(tf.random_normal([hl1_size, tileSizeHigh * tileSizeHigh], stddev=0.01))
fc_out_bias   = tf.Variable(tf.random_normal([tileSizeHigh * tileSizeHigh], stddev=0.01))

y_pred = tf.add(tf.matmul(fc1, fc_out_weight), fc_out_bias)
y_pred = tf.reshape( y_pred, shape=[-1, 1, tileSizeHigh, tileSizeHigh, 1])

costFunc = tf.nn.l2_loss(y_true - y_pred) 
optimizer = tf.train.AdamOptimizer(learningRate).minimize(costFunc)

# create session and saver
sess = tf.InteractiveSession()
saver = tf.train.Saver()

# init vars or load model
if loadModelTest == -1:
	sess.run(tf.global_variables_initializer())
else:
	saver.restore(sess, load_path)
	print("Model restored from %s." % load_path)

if outputPressure:
	outputDataName = "pressure"

# load test data
tiCr.loadTestDataNpz(fromSim, toSim, emptyTileValue, cropTileSizeLow, cropOverlap, 0.95, 0.05, load_vel=useVelocities, low_res_size=simSizeLow, upres=upRes, keepAll=keepAll, special_output_type=outputDataName, bWidth=bWidth)

if useVelocities and not useDensity:
	tiCr.reduceInputsToVelocity(dimensions=3)
	tiCr.splitTileData(0.95, 0.05)
	
# create a summary to monitor cost tensor
lossTrain  = tf.summary.scalar("loss",     costFunc)
lossTest   = tf.summary.scalar("lossTest", costFunc)
merged_summary_op = tf.summary.merge_all() 
summary_writer    = tf.summary.FileWriter(test_path, sess.graph)

# ---------------------------------------------
# START TRAINING
training_duration = 0.0
cost = 0.0
save_no = 0

if not outputOnly:
	try:
		print('\n*****TRAINING STARTED*****\n')
		print('(stop with ctrl-c)')
		avgCost = 0
		startTime = time.time()
		epochTime = startTime
		lastSave = 1
		lastCost = 1e10
		for epoch in range(trainingEpochs):
			batch_xs, batch_ys = tiCr.selectRandomTiles(batchSize)
			_, cost, summary = sess.run([optimizer, costFunc, lossTrain], feed_dict={x: batch_xs, y_true: batch_ys, training:True})

			# save model
			if ((cost < lastCost) or alwaysSave) and (lastSave >= saveInterval):
				saver.save(sess, test_path + 'model_%04d.ckpt' % save_no)
				save_no += 1
				lastSave = 1
				lastCost = cost
				print('Saved Model with cost %f.' % cost)
			else:
				lastSave += 1

			# display error
			avgCost += cost
			if (epoch + 1) % testInterval == 0:
				avgCostTest = 0.0
				for curr in range(numTests):
					batch_xs, batch_ys = tiCr.selectRandomTiles(batchSize, isTraining=False)
					cost_test, summary_test = sess.run([costFunc, lossTest], feed_dict={x: batch_xs, y_true: batch_ys, training:False}) 
					avgCostTest += cost_test
				avgCostTest /= numTests

				avgCost /= testInterval
				print('\nEpoch {:04d}/{:04d} - Avg. cost= {:.9f} - Avg. validation cost= {:.9f}'.format((epoch + 1), trainingEpochs, avgCost, avgCostTest))
				print('%d epochs took %.02f seconds.' % (testInterval, (time.time() - epochTime)))
				#print('Estimated time: %.02f minutes.' % ((trainingEpochs - epoch) / testInterval * (time.time() - epochTime) / 60.0))
				epochTime = time.time()
				summary_writer.add_summary(summary, epoch)
				summary_writer.add_summary(summary_test, epoch)
				avgCost = 0

	except KeyboardInterrupt:
		pass

	print('\n*****TRAINING FINISHED*****')
	training_duration = (time.time() - startTime) / 60.0
	print('Training needed %.02f minutes.' % (training_duration))
	print('To apply the trained model, set "outputOnly" to True, add "out 1 loadModelTest %d" to script call ' % test_folder_no)

else: 
	# ---------------------------------------------
	# outputOnly: apply to a full data set, and re-create full outputs from tiles

	batch_xs, batch_ys = tiCr.tile_inputs_all_complete, tiCr.tile_outputs_all_complete
	print('Creating %d tile outputs...' % (len(batch_xs)) )

	tileSizeHiCrop = upRes * cropTileSizeLow
	batchSize = tilesPerImg = (simSizeHigh // tileSizeHiCrop) ** 2 # note - override batchsize for output
	outrange = int( len(tiCr.tile_inputs_all_complete) / tilesPerImg )

	for currOut in range(outrange): 
		batch_xs = []
		batch_ys = []
		batch_velocity = [] # for optional output of velocity input pngs

		# use batchSize if its a multiple of tilesPerImage. Important for conv_trans networks.
		stopOutput = False
		for curr_tile in range(batchSize):
			idx = currOut * tilesPerImg + curr_tile
			if idx > len(tiCr.tile_inputs_all_complete) - 1:
				print('Warning - not all tiles used for output')
				stopOutput = True
				break
			batch_xs.append(tiCr.tile_inputs_all_complete[idx])
			batch_ys.append(tiCr.tile_outputs_all_complete[idx])

			# to output velocity inputs
			if useVelocities and outputInputs:
				batch_velocity.append(tiCr.tile_inputs_all_complete[idx])

		if stopOutput:
			break
		resultTiles = y_pred.eval(feed_dict={x: batch_xs, y_true: batch_ys, training:False}) 

		if brightenOutput > 0:
			for i in range(len(resultTiles)):
				resultTiles[i] *= brightenOutput
			for i in range(len(batch_ys)):
				batch_ys[i] *= brightenOutput
		tiCr.debugOutputPngsCrop(resultTiles, tileSizeHigh, simSizeHigh, test_path, imageCounter=currOut, cut_output_to=tileSizeHiCrop, tiles_in_image=tilesPerImg, name='output')
		tiCr.debugOutputPngsCrop(batch_ys,    tileSizeHigh, simSizeHigh, test_path, imageCounter=currOut, cut_output_to=tileSizeHiCrop, tiles_in_image=tilesPerImg, name='expected_out')

		if outputInputs:
			if not useVelocities:
				tiCr.debugOutputPngsSingle(batch_xs,         tileSizeLow, simSizeLow, test_path, imageCounter=currOut, name='input')
			else: 
				tiCr.debugOutputPngsSingle(batch_velocity, tileSizeLow, simSizeLow, test_path, imageCounter=currOut, name='in_vel_x', channel=1)
				tiCr.debugOutputPngsSingle(batch_velocity, tileSizeLow, simSizeLow, test_path, imageCounter=currOut, name='in_vel_y', channel=2)


# write summary to test overview
loaded_model = ''
if not loadModelTest == -1:
	loaded_model = ', Loaded %04d, %d' % (loadModelTest , loadModelNo)
with open(basePath + 'test_overview.txt', "a") as text_file:
	text_file.write(test_path[-10:-1] + ': {:.2f} min, {} Epochs, cost {:.4f}, {}'.format(training_duration, trainingEpochs, cost, " ") + loaded_model + '\n')


