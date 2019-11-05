#******************************************************************************
#
# tempoGAN: A Temporally Coherent, Volumetric GAN for Super-resolution Fluid Flow
# Copyright 2018 You Xie, Erik Franz, Mengyu Chu, Nils Thuerey
#
# This program is free software, distributed under the terms of the
# Apache License, Version 2.0 
# http://www.apache.org/licenses/LICENSE-2.0
#
#******************************************************************************

import time
import os
import shutil
import sys
import math

import tensorflow as tf
from tensorflow.python.client import timeline
import numpy as np

# load manta tools
sys.path.append("../tools")
import tilecreator_t as tc
import uniio
import paramhelpers as ph
from GAN import GAN, lrelu
import fluiddataloader as FDL

# ---------------------------------------------

# initialize parameters / command line params
outputOnly	  = int(ph.getParam( "out",			 False ))>0 		# output/generation mode, main mode switch

basePath		=	 ph.getParam( "basePath",		'../2ddata_gan/' )
randSeed		= int(ph.getParam( "randSeed",		1 )) 				# seed for np and tf initialization
load_model_test = int(ph.getParam( "load_model_test", -1 )) 			# the number of the test to load a model from. can be used in training and output mode. -1 to not load a model
load_model_no   = int(ph.getParam( "load_model_no",   -1 )) 			# nubmber of the model to load

simSizeLow  	= int(ph.getParam( "simSize", 		  64 )) 			# tiles of low res sim
tileSizeLow 	= int(ph.getParam( "tileSize", 		  16 )) 			# size of low res tiles
dt			= float(ph.getParam( "dt", 		  1.0 )) 				# step time of training data
#Data and Output
loadPath		 =	 ph.getParam( "loadPath",		 '../2ddata_sim/' ) 	# path to training data
fromSim		 = int(ph.getParam( "fromSim",		 1000 )) 			# range of sim data to use, start index
toSim		   = int(ph.getParam( "toSim",		   -1   )) 			# end index
dataDimension   = int(ph.getParam( "dataDim",		 2 )) 				# dimension of dataset, can be 2 or 3. in case of 3D any scaling will only be applied to H and W (DHW)
numOut			= int(ph.getParam( "numOut",		  200 )) 			# number ouf images to output (from start of sim)
saveOut	  	= int(ph.getParam( "saveOut",		 False ))>0 		# save output of output mode as .npz in addition to images
loadOut			= int(ph.getParam( "loadOut",		 -1 )) 			# load output from npz to use in output mode instead of tiles. number or output dir, -1 for not use output data
outputImages	=int(ph.getParam( "img",  			  True ))>0			# output images
outputGif		= int(ph.getParam( "gif",  			  False ))>0		# output gif
outputRef		= int(ph.getParam( "ref",			 False ))>0 		# output "real" data for reference in output mode (may not work with 3D)
#models
genModel		 =	 ph.getParam( "genModel",		 'gen_test' ) 	# choose generator model
discModel		=	 ph.getParam( "discModel",		 'disc_test' ) 	# choose discriminator model
#Training
learning_rate   = float(ph.getParam( "learningRate",  0.0002 ))
decayLR		= int(ph.getParam( "decayLR",			 False ))>0 		# decay learning rate?
dropout   		= float(ph.getParam( "dropout",  	  1.0 )) 			# keep prob for all dropout layers during training
dropoutOutput   = float(ph.getParam( "dropoutOutput", dropout )) 		# affects testing, full sim output and progressive output during training
beta			= float(ph.getParam( "adam_beta1",	 0.5 ))			#1. momentum of adam optimizer

weight_dld		= float(ph.getParam( "weight_dld",	1.0)) 			# ? discriminator loss factor ?
k				= float(ph.getParam( "lambda",		  1.0)) 			# influence/weight of l1 term on generator loss
k2				= float(ph.getParam( "lambda2",		  0.0)) 			# influence/weight of d_loss term on generator loss
k_f				= float(ph.getParam( "lambda_f",		  1.0)) 			# changing factor of k
k2_f			= float(ph.getParam( "lambda2_f",		  1.0)) 			# changing factor of k2
k2_l1		   = float(ph.getParam( "lambda2_l1",		   1.0))						 # influence/weight of L1 layer term on discriminator loss
k2_l2		   = float(ph.getParam( "lambda2_l2",		   1.0))						 # influence/weight of L2 layer term on discriminator loss
k2_l3		   = float(ph.getParam( "lambda2_l3",		   1.0))						 # influence/weight of L3 layer term on discriminator loss
k2_l4		   = float(ph.getParam( "lambda2_l4",		   1.0))						 # influence/weight of L4 layer term on discriminator loss
kt			  = float(ph.getParam("lambda_t", 1.0))				    # tempo discriminator loss; 1.0 is good, 0.0 will disable
kt_l		  = float(ph.getParam("lambda_t_l2", 0.0))				# l2 tempo loss (as naive alternative to discriminator); 1.0 is good, 0.0 will disable
batch_size	  = int(ph.getParam( "batchSize",  	  128 ))			# batch size for pretrainig and output, default for batchSizeDisc and batchSizeGen
batch_size_disc = int(ph.getParam( "batchSizeDisc",   batch_size )) 	# batch size for disc runs when training gan
batch_size_gen  = int(ph.getParam( "batchSizeGen",	batch_size )) 	# batch size for gen runs when training gan
trainGAN		= int(ph.getParam( "trainGAN",   	  True ))>0 		# GAN trainng can be switched off to use pretrainig only
trainingIters  = int(ph.getParam( "trainingIters",  100000 )) 		# for GAN training
discRuns 		= int(ph.getParam( "discRuns",  	  1 )) 				# number of discrimiinator optimizer runs per iteration
genRuns  		= int(ph.getParam( "genRuns",  		  1 )) 				# number or generator optimizer runs per iteration
batch_norm		= int(ph.getParam( "batchNorm",	   True ))>0			# apply batch normalization to conv and deconv layers
bn_decay		= float(ph.getParam( "bnDecay",	   0.999 ))			# decay of batch norm EMA
use_spatialdisc = int(ph.getParam( "use_spatialdisc",		   True )) #use spatial discriminator or not

useVelocities   = int(ph.getParam( "useVelocities",   0  )) 			# use velocities or not
useVorticities  = int(ph.getParam( "useVorticities",   0  )) 			# use vorticities or not
premadeTiles	= int(ph.getParam( "premadeTiles",   0  ))		 		# use pre-made tiles?

useDataAugmentation = int(ph.getParam( "dataAugmentation", 0 ))		 # use dataAugmentation or not
minScale = float(ph.getParam( "minScale",	  0.85 ))				 # augmentation params...
maxScale = float(ph.getParam( "maxScale",	  1.15 ))
rot	     = int(ph.getParam( "rot",		  2	 ))		#rot: 1: 90 degree rotations; 2: full rotation; else: nop rotation 
flip	 =   int(ph.getParam( "flip",		  1	 ))

#Test and Save
testPathStartNo = int(ph.getParam( "testPathStartNo", 0  ))
valiInterval	= int(ph.getParam( "valiInterval", 	  20  )) 			# interval in iterations to run validation, should be lower or equal outputInterval
numValis		= int(ph.getParam( "numValis", 		  10  )) 			# number of validation runs to perform from vali data each interval, run as batch
outputInterval	= int(ph.getParam( "outputInterval",  100  ))			# interval in iterations to output statistics
saveInterval	= int(ph.getParam( "saveInterval",	  200  ))	 		# interval in iterations to save model
alwaysSave	    = int(ph.getParam( "alwaysSave",	  True  )) 			#
maxToKeep		= int(ph.getParam( "keepMax",		 3  )) 			# maximum number of model saves to keep in each test-run
genValiImg		= int(ph.getParam( "genValiImg",	  -1 )) 			# if > -1 generate validation image every output interval
note			= ph.getParam( "note",		   "" )					# optional info about the current test run, printed in log and overview
data_fraction	= float(ph.getParam( "data_fraction",		   0.3 ))
frameMax		= int(ph.getParam( "frameMax",		   200 ))
frameMin		= int(ph.getParam( "frameMin",		   0 ))
ADV_flag		= int(ph.getParam( "adv_flag",		   True )) # Tempo parameter, add( or not) advection to pre/back frame to align
saveMD          = int(ph.getParam( "saveMetaData", 0 ))      # profiling, add metadata to summary object? warning - only main training for now
overlap         = int(ph.getParam( "overlap",		   3 )) # parameter for 3d unifile output, overlap of voxels

ph.checkUnusedParams()

useTempoD = False
useTempoL2 = False
if(kt > 1e-6):
	useTempoD = True
if(kt_l > 1e-6):
	useTempoL2 = True
if(kt > 1e-6 and kt_l > 1e-6):
	print("ERROR: temporal loss can only be either discriminator or L2, not both")
	exit(1)

# initialize
upRes	  		= 4 # fixed for now...
simSizeHigh 	= simSizeLow * upRes
tileSizeHigh	= tileSizeLow  * upRes

if not (dataDimension == 2 or dataDimension == 3):
	print('Unsupported data dimension {}. Only 2 and 3 are supported'.format(dataDimension))
	exit(1)

if toSim==-1:
	toSim = fromSim

channelLayout_low = 'd'
lowfilename = "density_low_%04d.uni"
highfilename = "density_high_%04d.uni"
mfl = ["density"]
mfh = ["density"]
if outputOnly: 
	highfilename = None
	mfh = None
if useVelocities:
	channelLayout_low += ',vx,vy,vz'
	mfl= np.append(mfl, "velocity")

dirIDs = np.linspace(fromSim, toSim, (toSim-fromSim+1),dtype='int16')

if (outputOnly): 
	data_fraction = 1.0
	kt = 0.0
	kt_l = 0.0
	useTempoD = False
	useTempoL2 = False
	useDataAugmentation = 0

if ((not useTempoD) and (not useTempoL2)): # should use the full sequence, not use multi_files
	tiCr = tc.TileCreator(tileSizeLow=tileSizeLow, simSizeLow=simSizeLow , dim =dataDimension, dim_t = 1, channelLayout_low = channelLayout_low, upres=upRes, premadeTiles=premadeTiles)
	floader = FDL.FluidDataLoader( print_info=1, base_path=loadPath, filename=lowfilename, oldNamingScheme=False, filename_y=highfilename, filename_index_min=frameMin, filename_index_max=frameMax, indices=dirIDs, data_fraction=data_fraction, multi_file_list=mfl, multi_file_list_y=mfh)
else:
	lowparalen = len(mfl)
	highparalen = len(mfh)
	mfl_tempo= np.append(mfl, mfl)
	mfl= np.append(mfl_tempo, mfl)
	mol = np.append(np.zeros(lowparalen), np.ones(lowparalen))
	mol = np.append(mol, np.ones(lowparalen)*2)
	mfh_tempo = np.append(mfh, mfh)
	mfh= np.append(mfh_tempo, mfh)
	moh = np.append(np.zeros(highparalen), np.ones(highparalen))
	moh = np.append(moh, np.ones(highparalen)*2)
	tiCr = tc.TileCreator(tileSizeLow=tileSizeLow, simSizeLow=simSizeLow , dim =dataDimension, dim_t = 3, channelLayout_low = channelLayout_low, upres=upRes, premadeTiles=premadeTiles)
	floader = FDL.FluidDataLoader( print_info=1, base_path=loadPath, filename=lowfilename, oldNamingScheme=False, filename_y=highfilename, filename_index_max=frameMax, indices=dirIDs, data_fraction=data_fraction, multi_file_list=mfl, multi_file_idxOff=mol, multi_file_list_y=mfh , multi_file_idxOff_y=moh)

if useDataAugmentation:
	tiCr.initDataAugmentation(rot=rot, minScale=minScale, maxScale=maxScale ,flip=flip)
inputx, y, xFilenames  = floader.get()
if (not outputOnly): 
	tiCr.addData(inputx,y)
elif dataDimension == 3:
	simLowLength = inputx.shape[1]
	simLowWidth = inputx.shape[2]
	simLowHeight = inputx.shape[3]

print("Random seed: {}".format(randSeed))
np.random.seed(randSeed)
tf.set_random_seed(randSeed)

# ---------------------------------------------

# 2D: tileSize x tileSize tiles; 3D: tileSize x tileSize x tileSize chunks
n_input  = tileSizeLow  ** 2
n_output = tileSizeHigh ** 2
if dataDimension == 3:
	n_input  *= tileSizeLow
	n_output *= (tileSizeLow*upRes)
n_inputChannels = 1

if useVelocities:
	n_inputChannels += 3
if useVorticities:
	n_inputChannels += 3
n_input *= n_inputChannels

# init paths
if not load_model_test == -1:
	if not os.path.exists(basePath + 'test_%04d/' % load_model_test):
		print('ERROR: Test to load does not exist.')
	load_path = basePath + 'test_%04d/model_%04d.ckpt' % (load_model_test, load_model_no)
	if outputOnly:
		out_path_prefix = 'out_%04d-%04d' % (load_model_test,load_model_no)
		test_path,_ = ph.getNextGenericPath(out_path_prefix, 0, basePath + 'test_%04d/' % load_model_test)

	else:
		test_path,_ = ph.getNextTestPath(testPathStartNo, basePath)

else:
	test_path, load_model_test_new = ph.getNextTestPath(testPathStartNo, basePath)

# logging & info
sys.stdout = ph.Logger(test_path)
print('Note: {}'.format(note))
print("\nCalled on machine '"+ os.uname()[1] +"' with: " + str(" ".join(sys.argv) ) )
print("\nUsing parameters:\n"+ph.paramsToString())
ph.writeParams(test_path+"params.json") # export parameters in human readable format

if outputOnly:
	print('*****OUTPUT ONLY*****')

if not outputOnly:
	os.makedirs(test_path+"/zbu_src")
	uniio.backupFile(__file__, test_path+"/zbu_src/")
	uniio.backupFile("../tools/tilecreator_t.py", test_path+"/zbu_src/")
	uniio.backupFile("../tools/GAN.py", test_path+"/zbu_src/") 
	uniio.backupFile("../tools/fluiddataloader.py", test_path+"/zbu_src/")

#input for gen
x = tf.placeholder(tf.float32, shape=[None, n_input])
#reference input for disc
x_disc = tf.placeholder(tf.float32, shape=[None, n_input])
#real input for disc
y = tf.placeholder(tf.float32, shape=[None, n_output])
kk = tf.placeholder(tf.float32)
kk2 = tf.placeholder(tf.float32)
kkt = tf.placeholder(tf.float32)
kktl = tf.placeholder(tf.float32)
# keep probablity for dropout
keep_prob = tf.placeholder(tf.float32)

print("x: {}".format(x.get_shape()))
# --- main graph setup ---

rbId = 0
def resBlock(gan, inp, s1,s2, reuse, use_batch_norm, filter_size=3):
	global rbId

	# convolutions of resnet block
	if dataDimension == 2:
		filter = [filter_size,filter_size]
		filter1 = [1,1]
	elif dataDimension == 3:
		filter = [filter_size,filter_size,filter_size]
		filter1 = [1,1,1]

	gc1,_ = gan.convolutional_layer(  s1, filter, tf.nn.relu, stride=[1], name="g_cA%d"%rbId, in_layer=inp, reuse=reuse, batch_norm=use_batch_norm, train=train) #->16,64
	gc2,_ = gan.convolutional_layer(  s2, filter, None      , stride=[1], name="g_cB%d"%rbId,               reuse=reuse, batch_norm=use_batch_norm, train=train) #->8,128

	# shortcut connection
	gs1,_ = gan.convolutional_layer(s2, filter1 , None       , stride=[1], name="g_s%d"%rbId, in_layer=inp, reuse=reuse, batch_norm=use_batch_norm, train=train) #->16,64
	resUnit1 = tf.nn.relu( tf.add( gc2, gs1 )  )
	rbId += 1
	return resUnit1

def gen_resnet(_in, reuse=False, use_batch_norm=False, train=None):
	global rbId
	print("\n\tGenerator (resize-resnett3-deep)")
	with tf.variable_scope("generator", reuse=reuse) as scope:

		if dataDimension == 2:
			_in = tf.reshape(_in, shape=[-1, tileSizeLow, tileSizeLow, n_inputChannels]) #NHWC
			patchShape = [2,2]
		elif dataDimension == 3:
			_in = tf.reshape(_in, shape=[-1, tileSizeLow, tileSizeLow, tileSizeLow, n_inputChannels]) #NDHWC
			patchShape = [2,2,2]
		rbId = 0
		gan = GAN(_in)
	
		gan.max_depool()
		inp = gan.max_depool()
		ru1 = resBlock(gan, inp, n_inputChannels*2,n_inputChannels*8,  reuse, use_batch_norm,5)

		ru2 = resBlock(gan, ru1, 128, 128,  reuse, use_batch_norm,5)
		inRu3 = ru2
		ru3 = resBlock(gan, inRu3, 32, 8,  reuse, use_batch_norm,5)
		ru4 = resBlock(gan, ru3, 2, 1,  reuse, False,5)
		resF = tf.reshape( ru4, shape=[-1, n_output] )
		print("\tDOFs: %d , %f m " % ( gan.getDOFs() , gan.getDOFs()/1000000.) ) 
		return resF

############################################discriminator network###############################################################
def disc_binclass(in_low, in_high, reuse=False, use_batch_norm=False, train=None):
	#in_low: low res reference input, same as generator input (condition)
	#in_high: real or generated high res input to classify
	#reuse: variable reuse
	#use_batch_norm: bool, if true batch norm is used in all but the first con layers
	#train: if use_batch_norm, tf bool placeholder
	print("\n\tDiscriminator (conditional binary classifier)")
	with tf.variable_scope("discriminator", reuse=reuse):
		if dataDimension == 2:
			shape = tf.shape(in_low)
			in_low = tf.slice(in_low,[0,0],[shape[0],int(n_input/n_inputChannels)])
			in_low = GAN(tf.reshape(in_low, shape=[-1, tileSizeLow, tileSizeLow, 1])).max_depool(height_factor = upRes,width_factor=upRes) #NHWC
			in_high = tf.reshape(in_high, shape=[-1, tileSizeHigh, tileSizeHigh, 1])
			filter=[4,4]
			stride = [2]
			stride2 = [2]
		elif dataDimension == 3:
			shape = tf.shape(in_low)
			in_low = tf.slice(in_low,[0,0],[shape[0],int(n_input/n_inputChannels)])
			in_low = GAN(tf.reshape(in_low, shape=[-1, tileSizeLow, tileSizeLow, tileSizeLow, 1])).max_depool(depth_factor = upRes,height_factor = upRes,width_factor = upRes) #NDHWC
			in_high = tf.reshape(in_high, shape=[-1, tileSizeHigh, tileSizeHigh, tileSizeHigh, 1]) # dim D is not upscaled
			filter=[4,4,4]
			stride = [2,2]
			stride2 = [2]

		#merge in_low and in_high to [-1, tileSizeHigh, tileSizeHigh, 2]
		gan = GAN(tf.concat([in_low, in_high], axis=-1), bn_decay=bn_decay) #64
		d1,_ = gan.convolutional_layer(32, filter, lrelu, stride=stride2, name="d_c1", reuse=reuse) #32

		d2,_ = gan.convolutional_layer(64, filter, lrelu, stride=stride2, name="d_c2", reuse=reuse, batch_norm=use_batch_norm, train=train) #64

		d3,_ = gan.convolutional_layer(128, filter, lrelu, stride=stride, name="d_c3", reuse=reuse, batch_norm=use_batch_norm, train=train) #128

		d4,_ = gan.convolutional_layer(256, filter, lrelu, stride=[1], name="d_c4", reuse=reuse, batch_norm=use_batch_norm, train=train) #256

		shape=gan.flatten()
		gan.fully_connected_layer(1, None, name="d_l5")

		print("\tDOFs: %d " % gan.getDOFs())
		return gan.y(), d1, d2, d3, d4
		
############################################ Tempo discriminator network ############################################################
def disc_binclass_cond_tempo(in_high, n_t_channels=3, reuse=False, use_batch_norm=False, train=None):
	# NO in_low: low res reference input, same as generator input (no condition)
	# in_high: real or generated high res input to classify, shape should be batch, dim_z, dim_y, dim_x, channels
	# reuse: variable reuse
	# use_batch_norm: bool, if true batch norm is used in all but the first con layers
	# train: if use_batch_norm, tf bool placeholder
	print("\n\tDiscriminator for Tempo (conditional binary classifier)")
	print("\n\tTempo, nearby frames packed as channels, number %d" % n_t_channels)
	with tf.variable_scope("discriminatorTempo", reuse=reuse):
		if dataDimension == 2:
			in_high = tf.reshape(in_high, shape=[-1, tileSizeHigh, tileSizeHigh, n_t_channels])
			filter=[4,4]
			stride = [2]
			stride2 = [2]
		elif dataDimension == 3:
			in_high = tf.reshape(in_high, shape=[-1, tileSizeHigh, tileSizeHigh, tileSizeHigh, n_t_channels]) # dim D is not upscaled
			filter=[4,4,4]
			stride = [2,2]
			stride2 = [2]

		# merge in_low and in_high to [-1, tileSizeHigh, tileSizeHigh, 2]
		gan = GAN(in_high, bn_decay=bn_decay)  # 64
		t1, _ = gan.convolutional_layer(32, filter, lrelu, stride=stride2, name="t_c1", reuse=reuse)  # 32
		t2, _ = gan.convolutional_layer(64, filter, lrelu, stride=stride2, name="t_c2", reuse=reuse,
										batch_norm=use_batch_norm, train=train)  # 64
		t3, _ = gan.convolutional_layer(128, filter, lrelu, stride=stride, name="t_c3", reuse=reuse,
										batch_norm=use_batch_norm, train=train)  # 128
		t4, _ = gan.convolutional_layer(256, filter, lrelu, stride=[1], name="t_c4", reuse=reuse,
										batch_norm=use_batch_norm, train=train)  # 256
		shape = gan.flatten()
		gan.fully_connected_layer(1, None, name="t_l5")

		print("\tDOFs: %d " % gan.getDOFs())
		return gan.y()

############################################gen_test###############################################################
def gen_test(_in, reuse=False, use_batch_norm=False, train=None):
	global rbId
	print("\n\tGenerator-test")
	with tf.variable_scope("generator-test", reuse=reuse) as scope:
		if dataDimension == 2:
			_in = tf.reshape(_in, shape=[-1, tileSizeLow, tileSizeLow, n_inputChannels]) #NHWC
			patchShape = [2,2]
		elif dataDimension == 3:
			_in = tf.reshape(_in, shape=[-1, tileSizeLow, tileSizeLow, tileSizeLow, n_inputChannels]) #NDHWC
			patchShape = [2,2,2]
		rbId = 0
		gan = GAN(_in)

		gan.max_depool()
		i2np,_ = gan.deconvolutional_layer(32, patchShape, None, stride=[1,1], name="g_D1", reuse=reuse, batch_norm=False, train=train, init_mean=0.99)
		gan.max_depool()
		inp,_  = gan.deconvolutional_layer(1                   , patchShape, None, stride=[1,1], name="g_D2", reuse=reuse, batch_norm=False, train=train, init_mean=0.99)
		return 	tf.reshape( inp, shape=[-1, n_output] )

############################################disc_test###############################################################
def disc_test(in_low, in_high, reuse=False, use_batch_norm=False, train=None):
	print("\n\tDiscriminator-test")
	with tf.variable_scope("discriminator_test", reuse=reuse):
		if dataDimension == 2:
			shape = tf.shape(in_low)
			in_low = tf.slice(in_low,[0,0],[shape[0],int(n_input/n_inputChannels)])
			in_low = GAN(tf.reshape(in_low, shape=[-1, tileSizeLow, tileSizeLow, 1])).max_depool(height_factor = upRes,width_factor = upRes) #NHWC
			in_high = tf.reshape(in_high, shape=[-1, tileSizeHigh, tileSizeHigh, 1])
			filter=[4,4]
			stride2 = [2]
		elif dataDimension == 3:
			shape = tf.shape(in_low)
			in_low = tf.slice(in_low,[0,0],[shape[0],int(n_input/n_inputChannels)])
			in_low = GAN(tf.reshape(in_low, shape=[-1, tileSizeLow, tileSizeLow, tileSizeLow, 1])).max_depool(depth_factor = upRes,height_factor = upRes,width_factor = upRes) #NDHWC
			in_high = tf.reshape(in_high, shape=[-1, tileSizeHigh, tileSizeHigh, tileSizeHigh, 1]) # dim D is not upscaled
			filter=[4,4,4]
			stride2 = [2]

		#merge in_low and in_high to [-1, tileSizeHigh, tileSizeHigh, 2]
		gan = GAN(tf.concat([in_low, in_high], axis=-1), bn_decay=bn_decay) #64
		d1,_ = gan.convolutional_layer(32, filter, lrelu, stride=stride2, name="d_c1", reuse=reuse) #32
		shape=gan.flatten()
		gan.fully_connected_layer(1, None, name="d_l5")
		if dataDimension == 2:
			d2 = tf.constant(1., shape = [batch_size, tileSizeLow,tileSizeLow,64])
			d3 = tf.constant(1., shape = [batch_size, int(tileSizeLow/2),int(tileSizeLow/2),128])	
			d4 = tf.constant(1., shape = [batch_size, int(tileSizeLow/2),int(tileSizeLow/2),256])
		elif dataDimension == 3:
			d2 = tf.constant(1., shape = [batch_size, tileSizeLow,tileSizeLow,tileSizeLow,64])
			d3 = tf.constant(1., shape = [batch_size, int(tileSizeLow/2),int(tileSizeLow/2),int(tileSizeLow/2),128])	
			d4 = tf.constant(1., shape = [batch_size, int(tileSizeLow/2),int(tileSizeLow/2),int(tileSizeLow/2),256])
		print("\tDOFs: %d " % gan.getDOFs())
		return gan.y(), d1, d2, d3, d4

#change used models for gen and disc here #other models in NNmodels.py
gen_model = locals()[genModel]
disc_model = locals()[discModel]
disc_time_model = disc_binclass_cond_tempo # tempo dis currently fixed

#set up GAN structure
bn=batch_norm
#training or testing for batch norm
train = tf.placeholder(tf.bool)

if not outputOnly: #setup for training
	gen_part = gen_model(x, use_batch_norm=bn, train=train)
	if use_spatialdisc:
		disc, dy1, dy2, dy3, dy4 = disc_model(x_disc, y, use_batch_norm=bn, train=train)
		gen, gy1, gy2, gy3, gy4 = disc_model(x_disc, gen_part, reuse=True, use_batch_norm=bn, train=train)
	if genValiImg > -1: sampler = gen_part
else: #setup for generating output with trained model
	sampler = gen_model(x, use_batch_norm=bn, train=False)

sys.stdout.flush()

# ---------------------------------------------
# TENSORFLOW SETUP
	
# build the tensorflow graph for tensor(value) re-sampling (at pos)
# value shape (batch, ..., res_x2, res_x1, channels)
# pos shape (batch, ..., res_x2, res_x1, dim)
def tensorResample(value, pos, name='Resample'):
	with tf.name_scope(name) as scope:
		pos_shape = pos.get_shape().as_list()
		dim = len(pos_shape) - 2  # batch and channels are ignored
		assert (dim == pos_shape[-1])
		floors = tf.cast(tf.floor(pos - 0.5), tf.int32)
		ceils = floors + 1

		# clamp min
		floors = tf.maximum(floors, tf.zeros_like(floors))
		ceils = tf.maximum(ceils, tf.zeros_like(ceils))

		# clamp max
		floors = tf.minimum(floors, tf.constant(value.get_shape().as_list()[1:dim + 1], dtype=tf.int32) - 1)

		ceils = tf.minimum(ceils, tf.constant(value.get_shape().as_list()[1:dim + 1], dtype=tf.int32) - 1)

		_broadcaster = tf.ones_like(ceils)
		cell_value_list = []
		cell_weight_list = []
		for axis_x in range(int(pow(2, dim))):  # 3d, 0-7; 2d, 0-3;...
			condition_list = [bool(axis_x & int(pow(2, i))) for i in range(dim)]
			condition_ = (_broadcaster > 0) & condition_list
			axis_idx = tf.cast(
				tf.where(condition_, ceils, floors),
				tf.int32)

			# only support linear interpolation...
			axis_wei = 1.0 - tf.abs((pos - 0.5) - tf.cast(axis_idx, tf.float32))  # shape (..., res_x2, res_x1, dim)
			axis_wei = tf.reduce_prod(axis_wei, axis=-1, keep_dims=True)
			cell_weight_list.append(axis_wei)  # single scalar(..., res_x2, res_x1, 1)
			first_idx = tf.ones_like(axis_wei, dtype=tf.int32)
			first_idx = tf.cumsum(first_idx, axis=0, exclusive=True)
			cell_value_list.append(tf.concat([first_idx, axis_idx], -1))

		values_new = tf.gather_nd(value, cell_value_list[0]) * cell_weight_list[
			0]  # broadcasting used, shape (..., res_x2, res_x1, channels )
		for cell_idx in range(1, len(cell_value_list)):
			values_new = values_new + tf.gather_nd(value, cell_value_list[cell_idx]) * cell_weight_list[cell_idx]
		return values_new  # shape (..., res_x2, res_x1, channels)

if not outputOnly:
	#for discriminator [0,1] output
	if use_spatialdisc:
		disc_sigmoid = tf.reduce_mean(tf.nn.sigmoid(disc))
		gen_sigmoid = tf.reduce_mean(tf.nn.sigmoid(gen))

		# loss of the discriminator with real input 
		disc_loss_disc = tf.reduce_mean(tf.nn.sigmoid_cross_entropy_with_logits(logits=disc, labels=tf.ones_like(disc)))
		#loss of the discriminator with input from generator
		disc_loss_gen = tf.reduce_mean(tf.nn.sigmoid_cross_entropy_with_logits(logits=gen, labels=tf.zeros_like(gen)))
		disc_loss_layer = k2_l1*tf.reduce_mean(tf.nn.l2_loss(dy1 - gy1)) + k2_l2*tf.reduce_mean(tf.nn.l2_loss(dy2 - gy2)) + k2_l3*tf.reduce_mean(tf.nn.l2_loss(dy3 - gy3)) + k2_l4*tf.reduce_mean(tf.nn.l2_loss(dy4 - gy4))
		disc_loss = disc_loss_disc * weight_dld + disc_loss_gen
		#loss of the generator
		gen_loss = tf.reduce_mean(tf.nn.sigmoid_cross_entropy_with_logits(logits=gen, labels=tf.ones_like(gen)))
	
	else:
		gen_loss = tf.zeros([1])
		disc_loss_layer = tf.zeros([1])

	#additional generator losses
	gen_l2_loss = tf.nn.l2_loss(y - gen_part)
	gen_l1_loss = tf.reduce_mean(tf.abs(y - gen_part)) #use mean to normalize w.r.t. output dims. tf.reduce_sum(tf.abs(y - gen_part))

	#uses sigmoid cross entropy and l1 - see cGAN paper
	gen_loss_complete = gen_loss + gen_l1_loss*kk + disc_loss_layer*kk2

	# set up decaying learning rate, if enabled
	lr_global_step = tf.Variable(0, trainable=False)
	learning_rate_scalar = learning_rate
	if decayLR:
		learning_rate = tf.train.polynomial_decay(learning_rate, lr_global_step, trainingIters//2, learning_rate_scalar*0.05, power=1.1)

	update_ops = tf.get_collection(tf.GraphKeys.UPDATE_OPS)
	gen_update_ops = update_ops[:]
	ori_gen_update_ops = update_ops[:]

	#variables to be used in the different otimization steps
	vars = tf.trainable_variables()
	g_var = [var for var in vars if "g_" in var.name]
	if use_spatialdisc:
		dis_update_ops = update_ops[:]
		d_var = [var for var in vars if "d_" in var.name]

	if (useTempoD or useTempoL2):# temporal loss here
		disT_update_ops = []
		ori_gen_loss_complete = gen_loss_complete
		# currently, the update_op gathering is not too nice and very sensitive to the operation order. 
		# TODO: make it flexible!
		n_t = 3
		device_str = '/device:GPU:0'
		if(dataDimension == 3): # have to use a second GPU!
			device_str = '/device:GPU:1'
		with tf.device(device_str): 
			x_t = tf.placeholder(tf.float32, shape=[None, n_input])
			gen_part_t = gen_model(x_t, reuse=True, use_batch_norm=bn, train=train)
		if(ADV_flag):
			y_pos = tf.placeholder(tf.float32, shape=[None, n_output * dataDimension])
			if dataDimension == 2:
				gen_part_t_shape = tf.reshape(gen_part_t, shape=[-1, tileSizeHigh, tileSizeHigh, 1])
				pos_array = tf.reshape(y_pos, shape=[-1, tileSizeHigh, tileSizeHigh, 2])
			elif dataDimension == 3:  # check in 3D
				gen_part_t_shape = tf.reshape(gen_part_t, shape=[-1, tileSizeHigh, tileSizeHigh, tileSizeHigh, 1])
				pos_array = tf.reshape(y_pos, shape=[-1, tileSizeHigh, tileSizeHigh, tileSizeHigh, 3])

			gen_part_t = tensorResample(gen_part_t_shape, pos_array)

		gen_part_t = tf.reshape(gen_part_t, shape = [-1, n_t, n_output])
		gen_part_t = tf.transpose(gen_part_t, perm=[0, 2, 1]) # batch, n_output, channels

		if (useTempoL2): # l2 tempo_loss
			update_ops = tf.get_collection(tf.GraphKeys.UPDATE_OPS)
			for update_op in update_ops:
				if ("/g_" in update_op.name) and ("generator" in update_op.name) and (not ( update_op in gen_update_ops )):
					gen_update_ops.append(update_op)
			
			gen_part_t_list = tf.unstack(gen_part_t, axis = -1) # should have n_t dim
			tl_gen_loss = tf.reduce_mean( tf.square( gen_part_t_list[0] - gen_part_t_list[1] ) )
			for ti in range( 1, n_t-1 ):
				tl_gen_loss = tl_gen_loss + tf.reduce_mean( tf.square( gen_part_t_list[ti] - gen_part_t_list[ti + 1] ) )
			gen_loss_complete = gen_loss_complete + kktl * tl_gen_loss
			
		if (useTempoD):
			# real input for disc
			y_t = tf.placeholder(tf.float32, shape=[None, n_output])
			if (ADV_flag):
				if dataDimension == 2:
					y_t_shape = tf.reshape(y_t, shape=[-1, tileSizeHigh, tileSizeHigh, 1])
				elif dataDimension == 3:  # check in 3D
					y_t_shape = tf.reshape(y_t, shape=[-1, tileSizeHigh, tileSizeHigh, tileSizeHigh, 1])
				y_tR = tensorResample(y_t_shape, pos_array)
			else:
				y_tR = y_t
			y_tR =tf.reshape(y_tR, shape = [-1, n_t, n_output])
			y_tR = tf.transpose(y_tR, perm=[0, 2, 1]) # batch, n_output, channels
			
			gen_t = disc_time_model(gen_part_t, n_t_channels = n_t, use_batch_norm=bn, train=train)
			
			update_ops = tf.get_collection(tf.GraphKeys.UPDATE_OPS)

			for update_op in update_ops:
				if ("/g_" in update_op.name) and ("generator" in update_op.name) and (not ( update_op in gen_update_ops )):
					gen_update_ops.append(update_op)
					disT_update_ops.append(update_op)
					
				if ("/t_" in update_op.name) and ("discriminatorTempo" in update_op.name):
					gen_update_ops.append(update_op)
				
				
			# discrimiinator for tempo only
			disc_t = disc_time_model(y_tR, n_t_channels = n_t, reuse=True, use_batch_norm=bn, train=train)
			
			update_ops = tf.get_collection(tf.GraphKeys.UPDATE_OPS)
			for update_op in update_ops:
				if ("/t_" in update_op.name) and ("discriminatorTempo" in update_op.name):
					disT_update_ops.append(update_op)
					
			t_disc_sigmoid = tf.reduce_mean(tf.nn.sigmoid(disc_t))
			t_gen_sigmoid = tf.reduce_mean(tf.nn.sigmoid(gen_t))
			vars = tf.trainable_variables()
			t_var = [var for var in vars if "t_" in var.name]
			# loss of the discriminator with real input
			t_disc_loss_disc = tf.reduce_mean(
				tf.nn.sigmoid_cross_entropy_with_logits(logits=disc_t, labels=tf.ones_like(disc_t)))
			# loss of the discriminator with input from generator
			t_disc_loss_gen = tf.reduce_mean(
				tf.nn.sigmoid_cross_entropy_with_logits(logits=gen_t, labels=tf.zeros_like(gen_t)))
			t_disc_loss = t_disc_loss_disc * weight_dld + t_disc_loss_gen
			
			t_gen_loss = tf.reduce_mean(tf.nn.sigmoid_cross_entropy_with_logits(logits=gen_t, labels=tf.ones_like(gen_t)))
			
			gen_loss_complete = gen_loss_complete + kkt * t_gen_loss
			
			with tf.control_dependencies(disT_update_ops):
				t_disc_optimizer = tf.train.AdamOptimizer(learning_rate, beta1=beta).minimize(t_disc_loss, var_list=t_var)
				
			with tf.control_dependencies(ori_gen_update_ops):
				# optimizer for generator, can only change variables of the generator,
				ori_gen_optimizer = tf.train.AdamOptimizer(learning_rate, beta1=beta).minimize(ori_gen_loss_complete, var_list=g_var)
	if use_spatialdisc:
		with tf.control_dependencies(dis_update_ops):
			#optimizer for discriminator, uses combined loss, can only change variables of the disriminator
			disc_optimizer_adam = tf.train.AdamOptimizer(learning_rate, beta1=beta)
			disc_optimizer = disc_optimizer_adam.minimize(disc_loss, var_list=d_var)

	with tf.control_dependencies(gen_update_ops):
		# optimizer for generator, can only change variables of the generator,
		gen_optimizer = tf.train.AdamOptimizer(learning_rate, beta1=beta).minimize(gen_loss_complete, var_list=g_var)

# create session and saver
config = tf.ConfigProto(allow_soft_placement=True)
sess = tf.InteractiveSession(config = config)
saver = tf.train.Saver(max_to_keep=maxToKeep)

# init vars or load model
if load_model_test == -1:
	sess.run(tf.global_variables_initializer())
else:
	saver.restore(sess, load_path)
	print("Model restored from %s." % load_path)

if not outputOnly:
	# create a summary to monitor cost tensor
	#training losses
	if use_spatialdisc:
		lossTrain_disc  = tf.summary.scalar("discriminator-loss train",     disc_loss)
		lossTrain_gen  = tf.summary.scalar("generator-loss train",     gen_loss)

	#testing losses
	if use_spatialdisc:
		lossVali_disc_disc   = tf.summary.scalar("discriminator-loss vali real", disc_loss_disc)
		lossVali_disc_gen   = tf.summary.scalar("discriminator-loss vali generated", disc_loss_gen)
		lossVali_disc = tf.summary.scalar("discriminator-loss vali", disc_loss)
		lossVali_gen   = tf.summary.scalar("generator-loss vali", gen_loss)

	#discriminator output [0,1] for real input
	if use_spatialdisc:
		outTrain_disc_real = tf.summary.scalar("discriminator-out train", disc_sigmoid)
		outTrain_disc_gen = tf.summary.scalar("generator-out train", gen_sigmoid)

	#discriminator output [0,1] for generated input
	if use_spatialdisc:
		outVali_disc_real = tf.summary.scalar("discriminator-out vali", disc_sigmoid)
		outVali_disc_gen = tf.summary.scalar("generator-out vali", gen_sigmoid)
	
	if(useTempoD): # all temporal losses
		# training losses, disc & gen
		lossTrain_disc_t = tf.summary.scalar("T discriminator-loss train", t_disc_loss)
		lossTrain_gen_t = tf.summary.scalar("T generator-loss train", t_gen_loss)
		
		# validation losses, discriminator( positive, negative ), generator
		lossVali_disc_disc_t = tf.summary.scalar("T discriminator-loss vali real", t_disc_loss_disc)
		lossVali_disc_gen_t = tf.summary.scalar("T discriminator-loss vali generated", t_disc_loss_gen)
		lossVali_disc_t = tf.summary.scalar("T discriminator-loss vali", t_disc_loss)
		lossVali_gen_t = tf.summary.scalar("T generator-loss vali", t_gen_loss)

		# discriminator output [0,1] for real input, during training
		outTrain_disc_real_t = tf.summary.scalar("T discriminator-out train", t_disc_sigmoid)
		# discriminator output [0,1] for generated input
		outTrain_disc_gen_t = tf.summary.scalar("T generator-out train", t_gen_sigmoid)

		# discriminator output [0,1] for real input, during validation
		outVali_disc_real_t = tf.summary.scalar("T discriminator-out vali", t_disc_sigmoid)
		# discriminator output [0,1] for generated input
		outVali_disc_gen_t = tf.summary.scalar("T generator-out vali", t_gen_sigmoid)
	
	if (useTempoL2):  # all temporal losses
		lossTrain_gen_t_l = tf.summary.scalar("T generator-loss train l2", tl_gen_loss)
		lossVali_gen_t_l = tf.summary.scalar("T generator-loss vali l2", tl_gen_loss)

	merged_summary_op = tf.summary.merge_all()
	summary_writer    = tf.summary.FileWriter(test_path, sess.graph)

save_no = 0
tileSizeHi = upRes * tileSizeLow
if dataDimension == 2:
	tilesPerImg = (simSizeHigh // tileSizeHi) ** 2
else:
	tilesPerImg = (simSizeHigh // tileSizeHi) ** 3
image_no = 0
if not outputOnly:
	os.makedirs(test_path+'test_img/')

def addVorticity(Vel):
	if dataDimension == 2:
		vorout = np.zeros_like(Vel)
		for l in range(vorout.shape[0]):
			for i in range(1, vorout.shape[-3]-1):
				for j in range(1, vorout.shape[-2]-1):
					vorout[l][0][i][j][2] = 0.5 * ((Vel[l][0][i+1][j][1] - Vel[l][0][i-1][j][1]) - (Vel[l][0][i][j+1][0] - Vel[l][0][i][j-1][0]))
	else:
		vorout = np.zeros_like(Vel)
		for l in range(vorout.shape[0]):
			for i in range(1, vorout.shape[-4]-1):
				for j in range(1, vorout.shape[-3]-1):
					for k in range(1, vorout.shape[-2]-1):		
						vorout[l][i][j][k][0] = 0.5 * ((Vel[l][i][j+1][k][2] - Vel[l][i][j-1][k][2]) - (Vel[l][i][j][k+1][1] - Vel[l][i][j][k-1][1]))
						vorout[l][i][j][k][1] = 0.5 * ((Vel[l][i][j][k+1][0] - Vel[l][i][j][k-1][0]) - (Vel[l][i+1][j][k][2] - Vel[l][i-1][j][k][2]))
						vorout[l][i][j][k][2] = 0.5 * ((Vel[l][i+1][j][k][1] - Vel[l][i-1][j][k][1]) - (Vel[l][i][j+1][k][0] - Vel[l][i][j-1][k][0]))
	return vorout

def getInput(index = 1, randomtile = True, isTraining = True, batch_size = 1, useDataAugmentation = False, useVelocities = False, useVorticities = False):
	if randomtile == False:
		batch_xs, batch_ys = tiCr.getFrameTiles(index) 
	else:
		batch_xs, batch_ys = tiCr.selectRandomTiles(selectionSize = batch_size, augment=useDataAugmentation)	

	if useVelocities and useVorticities:
		Velinput = batch_xs[:,:,:,:,tiCr.c_lists[tc.DATA_KEY_LOW][tc.C_KEY_VELOCITY][0]]
		Vorinput = addVorticity(Velinput)
		batch_xs = np.concatenate((batch_xs, Vorinput), axis = 4)
	batch_xs = np.reshape(batch_xs, (-1, n_input))
	batch_ys = np.reshape(batch_ys, (-1, n_output))
	return batch_xs, batch_ys

def getTempoinput(batch_size = 1, isTraining = True, useDataAugmentation = False, useVelocities = False, useVorticities = False, n_t = 3, dt=1.0, adv_flag = 1.0):
	batch_xts, batch_yts, batch_y_pos = tiCr.selectRandomTempoTiles(batch_size, isTraining, useDataAugmentation, n_t, dt, adv_flag)
	if useVelocities and useVorticities:
		real_batch_sz = batch_xts.shape[0]
		if( dataDimension == 2):
			batch_xts = np.reshape(batch_xts,[real_batch_sz,1,tileSizeLow,tileSizeLow,-1])
		else:
			batch_xts = np.reshape(batch_xts,[real_batch_sz,tileSizeLow,tileSizeLow,tileSizeLow,-1])
		Velinput = batch_xts[:,:,:,:,tiCr.c_lists[tc.DATA_KEY_LOW][tc.C_KEY_VELOCITY][0]]
		Vorinput = addVorticity(Velinput)
		batch_xts = np.concatenate((batch_xts, Vorinput), axis = 4)
		batch_xts = np.reshape(batch_xts,[real_batch_sz, -1])
	return batch_xts, batch_yts, batch_y_pos

#evaluate the generator (sampler) on the first step of the first simulation and output result
def generateValiImage(sim_no = fromSim, frame_no = 1, outPath = test_path,imageindex = 0):
	if premadeTiles:
		#todo output for premadetiles
		pass
	else:
		if (not outputOnly):
			batch_xs, _ = getInput(randomtile = False, index = (sim_no-fromSim)*frameMax + frame_no, useVelocities = useVelocities, useVorticities = useVorticities)
		else:
			batch_xs = inputx[frame_no]
		resultTiles = []
		for tileno in range(batch_xs.shape[0]):
			batch_xs_in = np.reshape(batch_xs[tileno],[-1, n_input])
			results = sess.run(sampler, feed_dict={x: batch_xs_in, keep_prob: dropoutOutput, train: False})
			resultTiles.extend(results)
		resultTiles = np.array(resultTiles)
		if dataDimension == 2: # resultTiles may have a different size
			imgSz = int(resultTiles.shape[1]**(1.0/2) + 0.5)
			resultTiles = np.reshape(resultTiles,[resultTiles.shape[0],imgSz,imgSz, 1])
		else:
			imgSz = int(resultTiles.shape[1]**(1.0/3) + 0.5)
			resultTiles = np.reshape(resultTiles,[resultTiles.shape[0],imgSz,imgSz,imgSz])
		tiles_in_image=[int(simSizeHigh/tileSizeHigh),int(simSizeHigh/tileSizeHigh)]
		tc.savePngsGrayscale(resultTiles,outPath, imageCounter=(imageindex+frameMin), tiles_in_image=tiles_in_image)

def generate3DUni(sim_no = fromSim, frame_no = 1, outPath = test_path,imageindex = 0):
	if dataDimension == 2:
		print("ERROR: only for 3D Uni files output!")	
		exit(1)
	if premadeTiles:
		#todo output for premadetiles
		pass
	else:
		if (overlap*2 > tileSizeLow) or (tileSizeLow > simLowLength):
			print("Wrong parameters for 3d output!")	
			exit(1)
		batch_xs = inputx[frame_no]
		if useVelocities and useVorticities:
			batch_xs = np.reshape(batch_xs,[1,simLowLength,simLowWidth,simLowHeight,-1])
			Velinput = batch_xs[:,:,:,:,tiCr.c_lists[tc.DATA_KEY_LOW][tc.C_KEY_VELOCITY][0]]
			Vorinput = addVorticity(Velinput)
			batch_xs = np.concatenate((batch_xs, Vorinput), axis = 4)
			batch_xs = np.reshape(batch_xs,[simLowLength,simLowWidth,simLowHeight,-1])
		tiles = []
		batch_xs=np.reshape(batch_xs,[simLowLength,simLowWidth,simLowHeight,-1])

		lengthnum = ((simLowLength-overlap*2+tileSizeLow-overlap*2-1)//(tileSizeLow-overlap*2))
		widthnum = ((simLowWidth-overlap*2+tileSizeLow-overlap*2-1)//(tileSizeLow-overlap*2))
		heightnum = ((simLowHeight-overlap*2+tileSizeLow-overlap*2-1)//(tileSizeLow-overlap*2))

		for i in range(lengthnum):
			for j in range(widthnum):
				for k in range(heightnum):
					ifrom = (tileSizeLow-overlap*2)*i
					ito = (tileSizeLow-overlap*2)*i+tileSizeLow
					jfrom = (tileSizeLow-overlap*2)*j
					jto = (tileSizeLow-overlap*2)*j+tileSizeLow
					kfrom = (tileSizeLow-overlap*2)*k
					kto = (tileSizeLow-overlap*2)*k+tileSizeLow
					if ito >simLowLength:
						ifrom = simLowLength-tileSizeLow
						ito = simLowLength
					if jto >simLowWidth:
						jfrom = simLowWidth-tileSizeLow
						jto = simLowWidth
					if kto >simLowHeight:
						kfrom = simLowHeight-tileSizeLow
						kto = simLowHeight
					low = batch_xs[ifrom:ito, jfrom:jto, kfrom:kto, :]
					tiles.append(low)
		batch_xs = np.array(tiles)
		resultTiles = []
		for tileno in range(batch_xs.shape[0]):
			batch_xs_in = np.reshape(batch_xs[tileno],[-1, n_input])
			results = sess.run(sampler, feed_dict={x: batch_xs_in, keep_prob: dropoutOutput, train : False})
			results = np.array(results)
			resultTiles.extend(results)
		resultTiles = np.array(resultTiles)
		resulttiles = np.reshape(resultTiles,[resultTiles.shape[0],tileSizeHigh,tileSizeHigh,tileSizeHigh])
		high = np.zeros([simLowLength*upRes,simLowWidth*upRes,simLowHeight*upRes])
		for i in range(lengthnum):
			for j in range(widthnum):
				for k in range(heightnum):
					ihighfrom = (tileSizeLow-overlap*2)*upRes*(i-1)+(tileSizeLow-overlap)*upRes
					ihighto = ihighfrom + (tileSizeLow-overlap*2)*upRes
					jhighfrom = (tileSizeLow-overlap*2)*upRes*(j-1)+(tileSizeLow-overlap)*upRes
					jhighto = jhighfrom+(tileSizeLow-overlap*2)*upRes
					khighfrom = (tileSizeLow-overlap*2)*upRes*(k-1)+(tileSizeLow-overlap)*upRes
					khighto = khighfrom+(tileSizeLow-overlap*2)*upRes
					ifrom = overlap*upRes
					ito = (tileSizeLow-overlap)*upRes
					jfrom = overlap*upRes
					jto = (tileSizeLow-overlap)*upRes
					kfrom = overlap*upRes
					kto = (tileSizeLow-overlap)*upRes
					if i == 0:
						ifrom = 0
						ito = (tileSizeLow-overlap)*upRes
						ihighfrom = 0
						ihighto = (tileSizeLow-overlap)*upRes
					if j == 0:
						jfrom = 0
						jto = (tileSizeLow-overlap)*upRes
						jhighfrom = 0
						jhighto = (tileSizeLow-overlap)*upRes
					if k == 0:
						kfrom = 0
						kto = (tileSizeLow-overlap)*upRes
						khighfrom = 0
						khighto = (tileSizeLow-overlap)*upRes
					if i == lengthnum-1:
						ifrom = overlap*upRes
						ito = tileSizeLow*upRes
						ihighfrom = simLowLength*upRes-tileSizeLow*upRes+overlap*upRes
						ihighto = simLowLength*upRes
					if j == widthnum-1:
						jfrom = overlap*upRes
						jto = tileSizeLow*upRes
						jhighfrom = simLowWidth*upRes-tileSizeLow*upRes+overlap*upRes
						jhighto = simLowWidth*upRes
					if k == heightnum-1:
						kfrom = overlap*upRes
						kto = tileSizeLow*upRes
						khighfrom = simLowHeight*upRes-tileSizeLow*upRes+overlap*upRes
						khighto = simLowHeight*upRes
					high[ihighfrom: ihighto, jhighfrom:jhighto, khighfrom:khighto] = resulttiles[i*widthnum*heightnum+j*heightnum+k][ifrom:ito,jfrom:jto,kfrom:kto]

		high = np.reshape(high,[simLowLength*upRes,simLowWidth*upRes,simLowHeight*upRes])
		
		head, _ = uniio.readUni(loadPath + "sim_%04d/density_high_%04d.uni"%(sim_no,frame_no+frameMin))
		head['dimX'] = simLowHeight*upRes
		head['dimY'] = simLowWidth*upRes
		head['dimZ'] = simLowLength*upRes
		uniio.writeUni(outPath+'source_%04d.uni'%(frame_no+frameMin), head, high)

def saveModel(cost, exampleOut=-1, imgPath = test_path):
	global save_no
	saver.save(sess, test_path + 'model_%04d.ckpt' % save_no)
	msg = 'Saved Model %04d with cost %f.' % (save_no, cost)
	if exampleOut > -1:
		generateValiImage(imageindex = save_no, outPath = imgPath)
	save_no += 1
	return msg

# write summary to test overview
loaded_model = ''
if not load_model_test == -1:
	loaded_model = ', Loaded %04d, %04d' % (load_model_test , load_model_no)
with open(basePath + 'test_overview.log', "a") as text_file:
	if not outputOnly:
		text_file.write(test_path[-10:-1] + ': {}D, \"{}\"\n'.format(dataDimension, note))
		text_file.write('\t{} Iters, gen: {}, disc: {}'.format(trainingIters, gen_model.__name__, disc_model.__name__) + loaded_model + '\n')
		text_file.write('\tgen-runs: {}, disc-runs: {}, lambda: {}, dropout: {:.4f}({:.4f})'.format(genRuns, discRuns, k, dropout, dropoutOutput) + '\n')
	else:
		text_file.write('Output:' + loaded_model + ' (' + test_path[-28:-1] + ')\n')
		text_file.write('\ttile size: {}, seed: {}, dropout-out: {:.4f}'.format(tileSizeLow, randSeed, dropoutOutput) + '\n')

	
# ---------------------------------------------
# ---------------------------------------------
# START TRAINING
training_duration = 0.0
cost = 0.0

if not outputOnly and trainGAN:
	try:
		print('\n*****TRAINING STARTED*****\n')
		print('(stop with ctrl-c)')
		avgCost_disc = 0
		avgCost_gen = 0
		avgL1Cost_gen = 0
		avgOut_disc = 0
		avgOut_gen = 0

		avgValiCost_disc_real = 0
		avgValiCost_disc_gen = 0
		avgValiCost_gen = 0
		avgValiOut_disc_real = 0
		avgValiOut_disc_gen = 0
		validations = 0
		startTime = time.time()
		intervalTime = startTime
		lastOut = 1
		lastSave = 1
		lastCost = 1e10
		saved = False
		saveMsg = ''
		kkin = k
		kk2in = k2

		disc_cost = 0
		gen_cost = 0
		
		avgTemCost_gen = 0
		avgTemCost_gen_l = 0
		avgTemCost_disc = 0
		kktin = kt
		kktin_l = kt_l

		avgOut_disc_t = 0
		avgOut_gen_t = 0
		avgValiCost_disc_real_t = 0
		avgValiOut_disc_real_t = 0
		avgValiCost_disc_gen_t = 0
		avgValiOut_disc_gen_t = 0
		avgValiCost_gen_t = 0
		avgValiCost_gen_t_l = 0
		
		for iteration in range(trainingIters):
			lrgs = max(0, iteration-(trainingIters//2)) # LR counter, start decay at half time... (if enabled) 
			run_options = None; run_metadata = None
			if saveMD:
				run_options = tf.RunOptions(trace_level=tf.RunOptions.FULL_TRACE)
				run_metadata = tf.RunMetadata()


			# TRAIN MODEL
			# discriminator variables; with real and generated input
			if use_spatialdisc:
				for runs in range(discRuns):
					batch_xs, batch_ys = getInput(batch_size = batch_size_disc, useDataAugmentation = useDataAugmentation, useVelocities = useVelocities, useVorticities = useVorticities)
					_, disc_cost, summary,disc_sig,gen_sig = sess.run([disc_optimizer, disc_loss, lossTrain_disc,disc_sigmoid,gen_sigmoid], feed_dict={x: batch_xs, x_disc: batch_xs, y: batch_ys, keep_prob: dropout, train: True, lr_global_step: lrgs}     , options=run_options, run_metadata=run_metadata )
					avgCost_disc += disc_cost
					summary_writer.add_summary(summary, iteration)
					if saveMD: summary_writer.add_run_metadata(run_metadata, 'dstep%d' % iteration)

			# temporal discriminator
			if(useTempoD):
				for runs in range(discRuns):
					batch_xts, batch_yts, batch_y_pos = getTempoinput(batch_size_disc, n_t = 3, dt=dt, useVelocities = useVelocities, useVorticities = useVorticities, useDataAugmentation = useDataAugmentation, adv_flag = ADV_flag)
					dict_train = {x_t:batch_xts, y_t:batch_yts, keep_prob: dropout, train: True}
					if(ADV_flag): dict_train[y_pos] = batch_y_pos
					_, t_disc_cost, summary, t_disc_sig, t_gen_sig = sess.run(
						[t_disc_optimizer, t_disc_loss, lossTrain_disc_t, t_disc_sigmoid, t_gen_sigmoid], feed_dict=dict_train)
					avgTemCost_disc += t_disc_cost
					summary_writer.add_summary(summary, iteration)
					
			# generator variables
			for runs in range(genRuns):
				batch_xs, batch_ys = getInput(batch_size = batch_size_disc, useDataAugmentation = useDataAugmentation, useVelocities = useVelocities, useVorticities = useVorticities)
				kkin = k_f*kkin
				kk2in = k2_f*kk2in
				# TODO a decay for weights, kktin = kt_f * kktin (kt_f<1.0)
				
				train_dict = {x: batch_xs, x_disc: batch_xs, y: batch_ys, keep_prob: dropout, train: True, kk: kkin,
							  kk2: kk2in, lr_global_step: lrgs}
				if use_spatialdisc:
					getlist = [gen_optimizer, gen_loss, disc_loss_layer, gen_l1_loss, lossTrain_gen, gen_l2_loss]
				else:
					getlist = [gen_optimizer, gen_l1_loss, gen_l2_loss]
				if(useTempoD or useTempoL2):
					batch_xts, batch_yts, batch_y_pos = getTempoinput(batch_size_disc, n_t = 3, dt=dt, useVelocities = useVelocities, useVorticities = useVorticities, useDataAugmentation=useDataAugmentation, adv_flag = ADV_flag)
					train_dict[x_t] = batch_xts
					if(ADV_flag):
						train_dict[y_pos] = batch_y_pos
					if(useTempoD): 
						train_dict[kkt] = kktin
						getlist.append(t_gen_loss)
					if(useTempoL2): 
						train_dict[kktl] = kktin_l
						getlist.append(tl_gen_loss)

				result_list = sess.run(getlist, feed_dict=train_dict, options=run_options, run_metadata=run_metadata)
				if (useTempoD and (not useTempoL2)):
					if use_spatialdisc:
						_, gen_cost, layer_cost, gen_l1_cost, summary, gen_l2_cost, gen_tem_cost = result_list
					else:
						_, gen_l1_cost, gen_l2_cost, gen_tem_cost = result_list
					gen_tem_cost_l = 0
				elif ((not useTempoD) and useTempoL2):
					if use_spatialdisc:
						_, gen_cost, layer_cost, gen_l1_cost, summary, gen_l2_cost, gen_tem_cost_l = result_list
					else:
						_, gen_l1_cost, gen_l2_cost, gen_tem_cost_l = result_list
					gen_tem_cost = 0
				elif (useTempoD and useTempoL2):
					if use_spatialdisc:
						_, gen_cost, layer_cost, gen_l1_cost, summary, gen_l2_cost, gen_tem_cost, gen_tem_cost_l = result_list
					else:
						_, gen_l1_cost, gen_l2_cost, gen_tem_cost, gen_tem_cost_l = result_list
					
				else:
					if use_spatialdisc:
						_, gen_cost, layer_cost, gen_l1_cost, summary, gen_l2_cost = result_list
					else:
						_, gen_l1_cost, gen_l2_cost = result_list
					gen_tem_cost = 0
					gen_tem_cost_l = 0
				avgL1Cost_gen += gen_l1_cost
				avgTemCost_gen += gen_tem_cost
				avgTemCost_gen_l += gen_tem_cost_l
				if use_spatialdisc:
					avgCost_gen += gen_cost
					summary_writer.add_summary(summary, iteration)
				if saveMD: summary_writer.add_run_metadata(run_metadata, 'gstep%d' % iteration)


			# save model
			if ((disc_cost+gen_cost < lastCost) or alwaysSave) and (lastSave >= saveInterval):
				lastSave = 1
				lastCost = disc_cost+gen_cost
				saveMsg = saveModel(lastCost)
				saved = True
			else:
				lastSave += 1
				saved = False

			# validate model
			if (iteration + 1) % valiInterval == 0:
				if use_spatialdisc:
					# gather statistics from training
					batch_xs, batch_ys = getInput(batch_size = numValis, useVelocities = useVelocities, useVorticities = useVorticities)
					disc_out, summary_disc_out, gen_out, summary_gen_out = sess.run([disc_sigmoid, outTrain_disc_real, gen_sigmoid, outTrain_disc_gen], feed_dict={x: batch_xs, x_disc: batch_xs, y: batch_ys, keep_prob: dropout, train: False})
					summary_writer.add_summary(summary_disc_out, iteration)
					summary_writer.add_summary(summary_gen_out, iteration)
					avgOut_disc += disc_out
					avgOut_gen += gen_out

					# validation starts here...
					# get vali data
					batch_xs, batch_ys = getInput(batch_size = numValis, isTraining=False, useVelocities = useVelocities, useVorticities = useVorticities)
					#disc with real imput
					disc_out_real, summary_vali_out, disc_vali_cost_real, summary_vali = sess.run([disc_sigmoid, outVali_disc_real, disc_loss_disc, lossVali_disc_disc], feed_dict={x: batch_xs, x_disc: batch_xs, y: batch_ys, keep_prob: dropoutOutput, train: False})
					summary_writer.add_summary(summary_vali, iteration)
					summary_writer.add_summary(summary_vali_out, iteration)
					avgValiCost_disc_real += disc_vali_cost_real
					avgValiOut_disc_real += disc_out_real
					#disc with generated input
					disc_out_gen, summary_vali_out, disc_vali_cost_gen, summary_vali = sess.run([gen_sigmoid, outVali_disc_gen, disc_loss_gen, lossVali_disc_gen], feed_dict={x: batch_xs, x_disc: batch_xs, keep_prob: dropoutOutput, train: False})
					summary_writer.add_summary(summary_vali, iteration)
					summary_writer.add_summary(summary_vali_out, iteration)
					avgValiCost_disc_gen += disc_vali_cost_gen
					avgValiOut_disc_gen += disc_out_gen
				
				if(useTempoD): # temporal logs
					# T disc output with training data
					batch_xts, batch_yts, batch_y_pos = getTempoinput(numValis, useVelocities = useVelocities, useVorticities = useVorticities, n_t = 3, dt=dt, adv_flag = ADV_flag)
					vali_dict = {x_t: batch_xts, y_t: batch_yts, keep_prob: dropout, train: False}
					if(ADV_flag):
						vali_dict[y_pos] = batch_y_pos
					t_disc_out, summary_disc_out_t, t_gen_out, summary_gen_out_t = sess.run(
						[t_disc_sigmoid, outTrain_disc_real_t, t_gen_sigmoid, outTrain_disc_gen_t],
						feed_dict=vali_dict)
					summary_writer.add_summary(summary_disc_out_t, iteration)
					summary_writer.add_summary(summary_gen_out_t, iteration)
					avgOut_disc_t += t_disc_out
					avgOut_gen_t += t_gen_out

					# validation data
					batch_xts, batch_yts, batch_y_pos = getTempoinput(numValis, isTraining=False, useVelocities = useVelocities, useVorticities = useVorticities, n_t = 3, dt=dt, adv_flag = ADV_flag)
					# disc with real input
					vali_dict = {x_t: batch_xts, y_t: batch_yts, keep_prob: dropout, train: False}
					if(ADV_flag):
						vali_dict[y_pos] = batch_y_pos
					t_disc_out_real, summary_vali_out_t, t_disc_vali_cost_real, summary_vali_t = sess.run(
						[t_disc_sigmoid, outVali_disc_real_t, t_disc_loss_disc, lossVali_disc_disc_t],
						feed_dict=vali_dict)
					summary_writer.add_summary(summary_vali_t, iteration)
					summary_writer.add_summary(summary_vali_out_t, iteration)
					avgValiCost_disc_real_t += t_disc_vali_cost_real
					avgValiOut_disc_real_t += t_disc_out_real
					# disc with generated input
					vali_dict = {x_t: batch_xts, y_t: batch_yts, keep_prob: dropout, train: False}
					if(ADV_flag):
						vali_dict[y_pos] = batch_y_pos
					t_disc_out_gen, summary_vali_out_t, t_disc_vali_cost_gen, summary_vali_t = sess.run(
						[t_gen_sigmoid, outVali_disc_gen_t, t_disc_loss_gen, lossVali_disc_gen_t],
						feed_dict=vali_dict)
					summary_writer.add_summary(summary_vali_t, iteration)
					summary_writer.add_summary(summary_vali_out_t, iteration)
					avgValiCost_disc_gen_t += t_disc_vali_cost_gen
					avgValiOut_disc_gen_t += t_disc_out_gen
					
				#gen
				train_dict = {x: batch_xs, x_disc: batch_xs, keep_prob: dropoutOutput, train: False}
				if (useTempoD or useTempoL2):  # add tempo logs
					train_dict[x_t] = batch_xts
					if(ADV_flag):
						train_dict[y_pos] = batch_y_pos
					if (useTempoD):
						train_dict[kkt] = kktin
						if use_spatialdisc:
							gen_vali_cost, summary_vali, gen_tem_cost, summary_vali_gen \
								= sess.run([gen_loss, lossVali_gen, t_gen_loss, lossVali_gen_t], feed_dict=train_dict)
						else:
							gen_tem_cost, summary_vali_gen \
								= sess.run([t_gen_loss, lossVali_gen_t], feed_dict=train_dict)
						avgValiCost_gen_t += gen_tem_cost
					if (useTempoL2):
						train_dict[kktl] = kktin_l
						if use_spatialdisc:
							gen_vali_cost, summary_vali, gen_tem_cost, summary_vali_gen \
								= sess.run([gen_loss, lossVali_gen, tl_gen_loss, lossVali_gen_t_l], feed_dict=train_dict)
						else:
							gen_tem_cost, summary_vali_gen \
								= sess.run([tl_gen_loss, lossVali_gen_t_l], feed_dict=train_dict)
						avgValiCost_gen_t_l += gen_tem_cost
					summary_writer.add_summary(summary_vali_gen, iteration)

				else:
					if use_spatialdisc:
						gen_vali_cost, summary_vali = sess.run([gen_loss, lossVali_gen], feed_dict=train_dict)
				if use_spatialdisc:	
					summary_writer.add_summary(summary_vali, iteration)
					avgValiCost_gen += gen_vali_cost

				validations += 1

			# output statistics
			if (iteration + 1) % outputInterval == 0:
				# training average costs
				avgCost_disc /= (outputInterval * discRuns)
				avgCost_gen /= (outputInterval * genRuns)
				avgL1Cost_gen /= (outputInterval * genRuns)
				# validation average costs
				if not (validations == 0):
					avgOut_disc /= validations
					avgOut_gen /= validations
					avgValiCost_disc_real /= validations
					avgValiCost_disc_gen /= validations
					avgValiCost_gen /= validations
					avgValiOut_disc_real /= validations
					avgValiOut_disc_gen /= validations
					
				if(useTempoD):
					avgTemCost_gen /= (outputInterval * genRuns)
					avgTemCost_disc /= (outputInterval * discRuns)
					if( not validations == 0):
						avgOut_disc_t /= validations
						avgOut_gen_t /= validations
						avgValiCost_disc_real_t /= validations
						avgValiOut_disc_real_t /= validations
						avgValiCost_disc_gen_t /= validations
						avgValiOut_disc_gen_t /= validations
						avgValiCost_gen_t /= validations
						
				if (useTempoL2):
					avgTemCost_gen_l /= (outputInterval * genRuns)
					if (not validations == 0):
						avgValiCost_gen_t_l /= validations
						
				print('\nIter {:05d}/{}, Cost:'.format((iteration + 1), trainingIters))
				print('\tdisc: loss: train_loss={:.6f} - vali-real={:.6f} - vali-generated={:.6f}, out: train={:.6f} - vali={:.6f}'.
					format(avgCost_disc, avgValiCost_disc_real, avgValiCost_disc_gen, avgOut_disc, avgValiOut_disc_real))
				print('\tT D : loss[ -train (total={:.6f}), -vali (real&1={:.6f}) (generated&0={:.6f})]'.
					format(avgTemCost_disc, avgValiCost_disc_real_t, avgValiCost_disc_gen_t))
				print('\t	sigmoidout[ -vali (real&1={:.6f}) (generated&0={:.6f})'.
					format(avgValiOut_disc_real_t, avgValiOut_disc_gen_t))
				print('\t gen: loss: train={:.6f} - L1(*k)={:.3f} - vali={:.6f}, DS out: train={:.6f} - vali={:.6f}'
					.format(avgCost_gen, avgL1Cost_gen * k, avgValiCost_gen, avgOut_gen, avgValiOut_disc_gen))
				print('\t gen: loss[ -train (total Temp(*k)={:.6f}) -vali (total Temp(*k)={:.6f})], DT out: real={:.6f} - gen={:.6f}'
					.format(avgTemCost_gen * kt, avgValiCost_gen_t * kt, avgOut_disc_t, avgOut_gen_t))
				if use_spatialdisc:
					print('\tdisc: loss: disc=%f'%(disc_sig))
					print('\tgen: loss: gen=%f'%(gen_sig))
					print('\t layer_cost: %f'%(layer_cost))
				if(useTempoD):
					print('\tTdisc: loss: disc=%f' % (t_disc_sig))
					print('\tTgen: loss: gen=%f' % (t_gen_sig))
				
				
				if(useTempoD): print('\t tempo_cost: %f' % (gen_tem_cost))
				print('\t l1_cost: %f'%(gen_l1_cost))
				print('\t l2 tempo loss[ -train (total Temp(*k)={:.6f}) -vali (total Temp(*k)={:.6f})]'
					.format(avgTemCost_gen_l * kt_l, avgValiCost_gen_t_l * kt_l))
				
				iterTime = (time.time() - startTime) / (iteration + 1)
				print('\t{} Iterations took {:.2f} seconds. (Est. next: {})'.format(outputInterval, (time.time() - intervalTime), time.ctime(time.time() + outputInterval * iterTime)))
				remainingTime = (trainingIters - iteration) * iterTime
				print('\tEstimated remaining time: {:.2f} minutes. (Est. end: {})'.format(remainingTime / 60.0, time.ctime(time.time() + remainingTime)))
				if saved:
					print('\t' + saveMsg) # print save massage here for clarity
				if genValiImg > -1:
					generateValiImage(outPath = test_path+'test_img/', imageindex = image_no)
					image_no +=1
				sys.stdout.flush()
				intervalTime = time.time()
				avgCost_disc = 0
				avgCost_gen = 0
				avgL1Cost_gen = 0
				avgOut_disc = 0
				avgOut_gen = 0
				avgValiCost_disc_real = 0
				avgValiCost_disc_gen = 0
				avgValiCost_gen = 0
				avgValiOut_disc_real = 0
				avgValiOut_disc_gen = 0
				validations = 0
				lastOut = 0
				
				if(useTempoD):
					avgTemCost_gen = 0
					avgTemCost_disc = 0
					avgOut_disc_t = 0
					avgOut_gen_t = 0
					avgValiCost_disc_real_t = 0
					avgValiOut_disc_real_t = 0
					avgValiCost_disc_gen_t = 0
					avgValiOut_disc_gen_t = 0
					avgValiCost_gen_t = 0
					
				if (useTempoL2):
					avgTemCost_gen_l = 0
					avgValiCost_gen_t_l = 0

			lastOut +=1

	except KeyboardInterrupt:
		print("training interrupted")
		sys.stdout.flush()
		with open(basePath + 'test_overview.log', "a") as text_file:
			text_file.write('\ttraining interrupted after %d iters' % (iteration + 1) + '\n')

	print('\n*****TRAINING FINISHED*****')
	training_duration = (time.time() - startTime) / 60.0
	print('Training needed %.02f minutes.' % (training_duration))
	print('To apply the trained model, call the script with command line params: "out 1  load_model_test %d  load_model_no %d" ' % (load_model_test_new, (save_no-1)) )
	sys.stdout.flush()
	with open(basePath + 'test_overview.log', "a") as text_file:
		text_file.write('\ttraining duration: %.02f minutes' % training_duration + '\n')


### OUTPUT MODE ###

elif outputOnly: 
	print('*****OUTPUT ONLY*****')

	for layerno in range(0,frameMax-frameMin):
		print('Generating %d' % (layerno))
		if dataDimension == 2:
			generateValiImage(fromSim,layerno,outPath = test_path, imageindex = layerno)
		else:
			generate3DUni(fromSim,layerno,outPath = test_path, imageindex = layerno)

	if outputGif and dataDimension == 2: #write gif
		print("Writing gif")
		tc.pngs_to_gif(test_path, start_idx=frameMin, end_idx=frameMax)

	print('Test finished, %d outputs written to %s.' % (frameMax-frameMin, test_path) )

