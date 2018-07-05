#******************************************************************************
#
# MantaFlow fluid solver framework
# Copyright 2017 Daniel Hook, Nils Thuerey
#
# This program is free software, distributed under the terms of the
# Apache License, Version 2.0 
# http://www.apache.org/licenses/LICENSE-2.0
#
# Manta & tensor flow example with tiles using Keras
# (note, this is a slightly simplified version of the tf_train script, with a deeper network and better regularization)
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

import keras

# path to sim data, trained models and output are also saved here
basePath = '../data/'

# main mode switch:
outputOnly = True  # apply model, or run full training?
outputInputs = False							 

simSizeLow   = 64
tileSizeLow  = 16
upRes        = 4
useVelocities= 0       # optional, add velocity (x,y,z) as additional channels to input?

kerasChunk      = 200     # run this many iterations per keras fit call
emptyTileValue  = 0.01
learningRate    = 0.002
trainingEpochs  = 10000 # for large values, stop manualy with ctrl-c...
dropout         = 0.9   # slight...
batchSize       = 96
saveInterval    = 10
fromSim = toSim = -1
keepAll         = False
numTests        = 10      # evaluate on 10 data points from test data
randSeed        = 1


# dont use for training! for applying model, add overlap here if necessary (i.e., cropOverlap>0) 
# note:  cropTileSizeLow + (cropOverlap * 2) = tileSizeLow
cropOverlap     = 0
cropTileSizeLow = tileSizeLow - 2*cropOverlap



# ---------------------------------------------

# load an existing model when load_ values > -1
# when training , manually abort when it's good enough
# then enter test_XXXX id and model checkpoint ID below to load

loadModelTest = -1
loadModelNo   = -1
testPathStartNo = 1

# command line params, explanations mostly above with variables
outputOnly      = int(ph.getParam( "out",             outputOnly ))>0
trainingEpochs  = int(ph.getParam( "trainingEpochs",  trainingEpochs ))
loadModelTest   = int(ph.getParam( "loadModelTest",   loadModelTest))
loadModelNo     = int(ph.getParam( "loadModelNo",     loadModelNo))
basePath        =     ph.getParam( "basePath",        basePath        )
useVelocities   = int(ph.getParam( "useVelocities",   useVelocities  ))
testPathStartNo = int(ph.getParam( "testPathStartNo", testPathStartNo  ))
fromSim         = int(ph.getParam( "fromSim",         fromSim  )) 
toSim           = int(ph.getParam( "toSim",           toSim  ))
randSeed        = int(ph.getParam( "randSeed",        randSeed )) 
simSizeLow      = int(ph.getParam( "simSizeLow",      simSizeLow )) 
upRes           = int(ph.getParam( "upRes",           upRes ))
outputInputs    = int(ph.getParam( "outInputs",       outputInputs)) 
ph.checkUnusedParams()

# initialize
simSizeHigh  = simSizeLow   * upRes
tileSizeHigh = tileSizeLow  * upRes
if outputOnly: # dont discard
	emptyTileValue = -1.

if toSim==-1:
	toSim = fromSim
tiCr.setBasePath(basePath)

np.random.seed(randSeed)
tf.set_random_seed(randSeed)

if not outputOnly:
	# run train!
	loadModelTest = -1
	if fromSim==-1:
		fromSim = toSim  = 1000 # short default, use single sim

	if cropOverlap>0:
		print("Error - dont use cropOverlap != 0 for training...")
		exit(1)

else:
	keepAll = True

# ---------------------------------------------

# create model loading path
if not loadModelTest == -1:
	if not os.path.exists(basePath + 'test_%04d/' % loadModelTest):
		print('ERROR: Test to load does not exist.')
	# search for newest model if no loadModelNo is given
	if loadModelNo == -1:
		for currModel in range(0, 999):
			if os.path.isfile(basePath + 'test_%04d/model_%04d.kkpt' % (loadModelTest, currModel)):
				loadModelNo = currModel
		if loadModelNo == -1:
			print('ERROR: Model with id below 200 does not exist. Please specify model id as "loadModelNo".')
			exit()
		# print('Latest model: %d.' % loadModelNo)

load_path = basePath + 'test_%04d/model_%04d.kkpt' % (loadModelTest, loadModelNo)

(test_path,test_folder_no) = ph.getNextTestPath(testPathStartNo, basePath)
if not outputOnly: uniio.backupFile(__file__, test_path)

sys.stdout = ph.Logger(test_path)

print("Call: " + str(" ".join(sys.argv) ) )

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

n_input  = tileSizeLow  ** 2 
n_output = tileSizeHigh ** 2
n_inputChannels = 1 
if useVelocities:
	n_inputChannels = 4
n_input *= n_inputChannels

clFMs = int(8 / n_inputChannels)
#print( "inputs " + format(len( tiCr.tile_data['inputs_train']) ))

# ---------------------------------------------
# model

model = keras.models.Sequential()
model.add( keras.layers.Conv2D(clFMs/2, (2,2), activation='relu', strides=(2,2), input_shape=(tileSizeLow,tileSizeLow,n_inputChannels), padding='same' ) )

model.add( keras.layers.Conv2D(clFMs  , (2,2), activation='relu', strides=(2,2), padding='same' ) )
model.add( keras.layers.BatchNormalization() )  
model.add( keras.layers.Flatten() ) # not really needed

model.add( keras.layers.Dense(4*4*clFMs, activation='relu') )
model.add( keras.layers.BatchNormalization() )  
model.add( keras.layers.Dropout(0.25) )  
model.add( keras.layers.Reshape( (4,4,clFMs) ) ) # for flatten, not really needed

model.add( keras.layers.convolutional.Conv2DTranspose(clFMs,   (2,2), activation='relu', strides=(2,2), padding='same' ) )
model.add( keras.layers.convolutional.Conv2DTranspose(clFMs/2, (2,2), activation='relu', strides=(2,2), padding='same' ) )
model.add( keras.layers.convolutional.Conv2DTranspose(clFMs/4, (2,2), activation='relu', strides=(2,2), padding='same' ) )
model.add( keras.layers.convolutional.Conv2DTranspose(1,       (4,4), activation='relu', strides=(2,2), padding='same' , name="out") )

model.compile( loss='mse', optimizer=keras.optimizers.adam(lr=learningRate) )

if 1: # count DOFs?
	DOFs = 0
	for l in model.layers:
		for i in l.get_weights():
			DOFs += i.size
	print("Total DOFs: " + format(DOFs) )

if not loadModelTest == -1:
	model.load_weights( load_path )
	print("Model restored from %s." % load_path)

#print( format( model.layers[ len(model.layers)-1 ].output )) # print output shape

# load test data, note no split into train & test/validation set necessary here, done by keras during fit
tiCr.loadTestDataNpz(fromSim, toSim, emptyTileValue, cropTileSizeLow, cropOverlap, 1.0, 0.0, load_vel=useVelocities, low_res_size=simSizeLow, upres=upRes, keepAll=keepAll)

# ---------------------------------------------
# training
if not outputOnly:
	# manually reshape data to remove third spatial dimension
	dataSize = tiCr.tile_data['outputs_train'].shape[0]
	#print( "Data sizes "+ format(dataSize) +", inputs " + format( tiCr.tile_data['inputs_train'].shape) + ", outputs " + format( tiCr.tile_data['outputs_train'].shape) ) # 16x16, 64x64 
	tiCr.tile_data['inputs_train']  = np.reshape( tiCr.tile_data['inputs_train'],  [dataSize,tileSizeLow,tileSizeLow,  n_inputChannels] )
	tiCr.tile_data['outputs_train'] = np.reshape( tiCr.tile_data['outputs_train'], [dataSize,tileSizeHigh,tileSizeHigh,1 ] )

	startTime = time.time()
	if 1:
		lastSave   = 1
		save_no    = 0
		trainingEpochs = int(trainingEpochs/kerasChunk)
		for epoch in range(trainingEpochs):
			batch_xs, batch_ys = tiCr.selectRandomTiles(batchSize*kerasChunk) 
			batch_xs  = np.reshape( batch_xs,  [batchSize*kerasChunk,tileSizeLow,tileSizeLow,  n_inputChannels] )
			batch_ys  = np.reshape( batch_ys,  [batchSize*kerasChunk,tileSizeHigh,tileSizeHigh,  1] )

			# ...go!
			hist = model.fit( batch_xs, batch_ys, batch_size=batchSize, epochs=1, validation_split=0.05 )

			# save model
			doSave = False
			if (lastSave >= saveInterval) or epoch == (trainingEpochs-1): doSave = True 

			if doSave:
				mpath = test_path + 'model_%04d.kkpt' % save_no
				model.save_weights(mpath)
				print('Saved weights ' + mpath )
				save_no += 1
				lastSave = 1
			else:
				lastSave += 1
			print('\nScript Epoch {:04d}/{:04d}'.format((epoch + 1), trainingEpochs))

	print('\n*****TRAINING %d FINISHED*****' % test_folder_no)
	training_duration = (time.time() - startTime) / 60.0
	#print('Training needed %.02f minutes.' % (training_duration))
	print('To apply the trained model, set "outputOnly" to True, and set id for "loadModelTest". E.g. "out 1 loadModelTest %d".' % test_folder_no)
	
else:
	# ---------------------------------------------
	# outputOnly: apply to a full data set, and re-create full outputs from tiles

	batch_xs, batch_ys = tiCr.tile_inputs_all_complete, tiCr.tile_outputs_all_complete # old, inputs_test, outputs_test
	print('Creating %d tile outputs...' % (len(batch_xs)) )

	tdataSize = len(batch_xs)
	# for prediction, we have to get rid of the third spatial dim
	batch_xs = np.reshape( batch_xs, [tdataSize,tileSizeLow,tileSizeLow,  n_inputChannels] )
	batch_ys = np.reshape( batch_ys, [tdataSize,tileSizeHigh,tileSizeHigh,1 ] )

	resultTiles = model.predict( batch_xs )

	# now restore it...
	batch_xs = np.reshape( batch_xs, [tdataSize,1, tileSizeLow,tileSizeLow,  n_inputChannels] )
	batch_ys = np.reshape( batch_ys, [tdataSize,1, tileSizeHigh,tileSizeHigh,1 ] )
	resultTiles = np.reshape( resultTiles, [tdataSize,1, tileSizeHigh,tileSizeHigh,1 ] )

	# simply concat tiles into images...
	tileSizeHiCrop = upRes * cropTileSizeLow
	tilesPerImg = (simSizeHigh // tileSizeHiCrop) ** 2
	imgCnt = len(tiCr.tile_inputs_all_complete) / tilesPerImg
	tiCr.debugOutputPngsCrop(resultTiles, tileSizeHigh, simSizeHigh, test_path, imageCounter=0, cut_output_to=tileSizeHiCrop, \
		tiles_in_image=tilesPerImg, name='output')

	if outputInputs:
		tiCr.debugOutputPngsSingle(batch_xs,         tileSizeLow, simSizeLow, test_path, imageCounter=0, name='input', channel=0)
		if useVelocities: 
			tiCr.debugOutputPngsSingle(batch_xs,         tileSizeLow, simSizeLow, test_path, imageCounter=0, name='in_vel_x', channel=1)
			tiCr.debugOutputPngsSingle(batch_xs,         tileSizeLow, simSizeLow, test_path, imageCounter=0, name='in_vel_y', channel=2)

	print('Output finished, %d pngs written to %s.' % (imgCnt, test_path) )



