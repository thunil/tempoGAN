#******************************************************************************
#
# tempoGAN: A Temporally Coherent, Volumetric GAN for Super-resolution Fluid Flow
# Copyright 2018 You Xie, Erik Franz, Mengyu Chu, Nils Thuerey
#
# This program is free software, distributed under the terms of the
# Apache License, Version 2.0 
# http://www.apache.org/licenses/LICENSE-2.0
#
# script to test data loading, tilecreator or output reference data
#
#******************************************************************************

import os,sys,time
import traceback
from itertools import repeat
import tilecreator_t as tc
import fluiddataloader as fdl
import numpy as np
import paramhelpers as ph


out_path		=	 ph.getParam( "basePath",		'../test_out/' )
sim_path		=	 ph.getParam( "simPath",		'../data_sim/' )
randSeed		= int(ph.getParam( "randSeed",		1 )) 				# seed for np and tf initialization

simSizeHigh  	= int(ph.getParam( "simSizeHigh", 		256 )) 			# size of high res sim
tileSizeHigh 	= int(ph.getParam( "tileSizeHigh", 		64 )) 			# size of high res tiles
simSizeLow  	= int(ph.getParam( "simSizeLow", 		64 )) 			# size of low res sim
tileSizeLow 	= int(ph.getParam( "tileSizeLow", 		16 )) 			# size of low res tiles
upRes	  		= float(ph.getParam( "upRes", 		4 )) 				# single generator scaling factor
dim   = int(ph.getParam( "dim",		 2 )) 								# dimension of dataset

useDummyData	= int(ph.getParam( "dummyData",		0 ))				# create dummy arrays instead of loading real data
dummySizeLow	= int(ph.getParam( "dummyLow",		simSizeLow ))				# sim size of dummy data if used. use to test size mismatches
dummySizeHigh	= int(ph.getParam( "dummyHigh",		simSizeHigh ))				#

useVel = int(ph.getParam( "vel", 1 ))									# currently not in use
augment = int(ph.getParam( "aug", 1 ))									# use dataAugmentation or not

# no high res data in TC, using high data in TC's low res
useScaledData = int(ph.getParam( "scaled", 0 ))
mainIsLow = int(ph.getParam( "mainIsLow", 1 if useScaledData else 0 ))							# use high or low data as main data. TileCreator requires the main data to be the smaller one.

useLabelData = int(ph.getParam( "label", 0 ))
useDataBlocks = int(ph.getParam( "block", 0 ))
blockSize = int(ph.getParam( "blockSize", 1 ))

batchCount = int(ph.getParam( "batchCount", 1 ))						# number of batches to create
saveImages = int(ph.getParam( "img", 1 ))
saveRef = int(ph.getParam( "ref", 0 ))

fail = int(ph.getParam( "fail", 0 ))									# wether the test is supposed to fail. used for output to indicate a "successful fail".

ph.checkUnusedParams()
np.random.seed(randSeed)
#tf.set_random_seed(randSeed)


#if not os.path.exists(out_path):
#	os.makedirs(out_path)
test_path,_ = ph.getNextTestPath(0, out_path)
sys.stdout = ph.Logger(test_path)
sys.stderr = ph.ErrorLogger(test_path)

def testFailed(msg):
	print('\n{}:\n{}'.format(msg, traceback.format_exc()))
	if fail:
		print('')
		print('--- TEST FAILED AS INDICATED ---')
		print('')
	else:
		print('')
		print('--- TEST FAILED ---')
		print('')
	exit()

print('')
print('--- TEST STARTED ---')
print('')

print("\nUsing parameters:\n"+ph.paramsToString())

recursionDepth = 0

if dim==2:
	fromSim = 1018
	index_min = 0
	index_max = 200
	fileType = '.npz'
	rgb_channels = [[1,2]]
	rgb_range = [-2,2]
elif dim==3:
	fromSim = 3006
	if useScaledData:
		index_min = 30
		index_max = 60
	else:
		index_min = 30
		index_max = 60
	fileType = '.npz'
	rgb_channels = [[1,2,3]]
	rgb_range = [-2,2]#[-0.1,0.1]
else:
	print('dim must be 2 or 3.')
	exit()
toSim = fromSim

dirIDs = np.linspace(fromSim, toSim, (toSim-fromSim+1),dtype='int16')
lowfilename = "density_low_%04d" + fileType
highfilename = "density_high_%04d" + fileType

# this is no longer supported by the TileCreator. main must always be the smaller one.
if mainIsLow:
	if not useScaledData:
		highfilename = None
		upRes = 1
	simSize = simSizeLow
	tileSize = tileSizeLow
else:
	lowfn = lowfilename
	lowfilename = highfilename
	if useScaledData:
		highfilename = lowfn
		upRes = 1/upRes
	else:
		highfilename = None
		upRes = 1
	simSize = simSizeHigh
	tileSize = tileSizeHigh
#if useVel:
#	fl = ["density", "velocity"]
#	ol = [0,0]
#else:

#load data
mfl = ["density", "velocity"] if not useDataBlocks else ["density", "velocity", "density", "velocity", "density", "velocity" ]
mol  = [0,0] if not useDataBlocks else [0,0,1,1,2,2]
mfh = None
moh = None
if useScaledData:
	mfh = ["density", "velocity"] if not useDataBlocks else ["density", "velocity", "density", "velocity", "density", "velocity" ]
	moh  = [0,0] if not useDataBlocks else [0,0,1,1,2,2]

if not useDummyData:
	print('\n  - LOADING DATA -\n')
	pt1_start = time.perf_counter()
	pt2_start = time.process_time()
	floader = fdl.FluidDataLoader( print_info=1, base_path=sim_path, filename=lowfilename, oldNamingScheme=False, filename_y=highfilename, filename_index_min=index_min, filename_index_max=index_max, indices=dirIDs, data_fraction=0.5, multi_file_list=mfl, multi_file_idxOff=mol, multi_file_list_y=mfh , multi_file_idxOff_y=moh)
	x, y, xFilenames  = floader.get()
	pt2_end = time.process_time()
	pt1_end = time.perf_counter()
	pt1 = pt1_end - pt1_start
	pt2 = pt2_end - pt2_start
	print('Loading Process Time: Total {:.04f}s; avg/frame {:.04f}s'.format(pt2,pt2/len(x)))
	print('Loading Time: Total {:.04f}s; avg/frame {:.04f}s'.format(pt1,pt1/len(x)))
	
else:
	print('\n  - CREATING DATA -\n')
	shapeLow = (40 ,(dummySizeLow if dim==3 else 1),dummySizeLow,dummySizeLow, dim+1)
	shapeHigh = (40 ,(dummySizeHigh if dim==3 else 1),dummySizeHigh,dummySizeHigh, dim+1)
	
	x = np.ones(shapeLow if mainIsLow else shapeHigh)
	if useScaledData:
		y = np.ones(shapeHigh if mainIsLow else shapeLow)


tile_format='NYXC'
z_axis = 1
if useDataBlocks or dim==3:
	tile_format='NBYXC'
	z_axis = 2

print('Loaded x shape: {}'.format(x.shape))
if useScaledData: print('Loaded y shape: {}'.format(y.shape))
l=None
if useLabelData:
	l = list(range(len(x)))
	print('\tUsing label data:\n{}'.format(l))
b=None
if useDataBlocks:
	b = [i for item in range(len(x)) for i in repeat(item, 3)]
	x = tc.blockFromChannelsToSequence(x, 3)
	print('\tUsing block data:\n{}'.format(b))
	print('Extracted blocks x shape: {}'.format(x.shape))
	if useScaledData:
		y = tc.blockFromChannelsToSequence(y, 3)
		print('Extracted blocks y shape: {}'.format(y.shape))
	if useLabelData:
		l = b[:]

#save ref:
if saveRef:
	print('Output reference')
	#tileShape = (x.shape[0],simSize,simSize,x.shape[-1])
	#tiles = np.reshape(x, tileShape) ,27:34
	tc.savePngs(x[11:12], test_path + 'ref_main_', tile_format=tile_format,imageCounter=0, tiles_in_image=[1,1], plot_vel_x_y=False, channels=[0], save_rgb = rgb_channels, rgb_interval=rgb_range)
	if useScaledData:
		#tileShape = (y.shape[0],simSize*upRes,simSize*upRes,y.shape[-1])
		#tiles = np.reshape(y, tileShape)
		tc.savePngs(y[:1], test_path + 'ref_scaled_',imageCounter=0, tiles_in_image=[1,1], plot_vel_x_y=False, channels=[0], save_rgb = [[2,3]], rgb_interval=[-2,2])

channel_layout = 'd,vx,vy'
if dim==3:
	channel_layout += ',vz'
# tilecreator
print('\n  - INIT TILECREATOR -\n')
try:
	TC = tc.TileCreator(tileSize=tileSize, simSize=simSize , dim=dim, densityMinimum=0.1, scaleFactor=upRes, channelLayout_main=channel_layout, channelLayout_scaled=channel_layout, useScaledData=useScaledData, useLabels=useLabelData, useDataBlocks=useDataBlocks, logLevel=10)
except tc.TilecreatorError as e:
	testFailed('TileCreator Error on construction')
if augment:
	print('\n  - INIT DATA AUGMENTATION -\n')
	try:
		TC.initDataAugmentation(2)
	except tc.TilecreatorError as e:
		testFailed('TileCreator Error on augmentation init')

# strip zero z vel of 2D data
if dim==2:
	x,_ = np.split(x, [3], axis=-1)
	if useScaledData:
		y,_ = np.split(y, [3], axis=-1)

# add low data with dummy labels
print('\n  - ADDING DATA -\n')
try:
	TC.addData(x, y if useScaledData else None, l, b)
except tc.TilecreatorError as e:
	testFailed('TileCreator Error when adding data')
#bx,by = TC.selectRandomTiles(64, True, augment=True)



#test batch:
if True:
	print('\n  - CREATING BATCH -\n')
	imageCounter=0
	pt1_start = time.perf_counter()
	pt2_start = time.process_time()
	for batch_number in range(batchCount):
		i=0
		print('\nOutput batch #{}'.format(batch_number))
		try:
			batch = TC.selectRandomTiles(selectionSize = 8, augment=augment, isTraining=True, blockSize=blockSize, squeezeZ=True)
		except tc.TilecreatorError as e:
			testFailed('TileCreator Error when creating batch')
		print('batch_x shape: {}'.format(batch[0].shape))
		
		#tileShape = (batch[0].shape[0],tileSize,tileSize,batch[0].shape[-1])
		tiles = batch[i]
		i+=1
		#print('tiles_x shape: {}'.format(tiles.shape))
		if dim==3: tiles = (tiles[2:6,6:8] if useScaledData else tiles[2:6,24:32])
		if saveImages: ic=tc.savePngs(tiles, test_path, tile_format=tile_format,imageCounter=imageCounter, tiles_in_image=[1,1], plot_vel_x_y=False, channels=[0], save_rgb = rgb_channels, rgb_interval=rgb_range)
		if useScaledData:
			print('batch_y shape: {}'.format(batch[1].shape))
			#tileShape = (batch[1].shape[0],tileSize*upRes,tileSize*upRes,batch[1].shape[-1])
			tiles = batch[i]
			i+=1
			#print('tiles_y shape: {}'.format(tiles.shape))
			if dim==3: tiles = tiles[2:6,24:32]
			if saveImages: tc.savePngs(tiles, test_path + 'high_', tile_format=tile_format, imageCounter=imageCounter, tiles_in_image=[1,1], plot_vel_x_y=False, channels=[0], save_rgb = rgb_channels, rgb_interval=rgb_range)
		if useLabelData:
			tiles = batch[i]
			i+=1
			print('labels: {}'.format(tiles))
		if saveImages: print('-> images {} to {}'.format(imageCounter, ic-1))
		if saveImages: imageCounter=ic
	
	pt2_end = time.process_time()
	pt1_end = time.perf_counter()
	pt1 = pt1_end - pt1_start
	pt2 = pt2_end - pt2_start
	print('Process Time: Total {:.04f}s; avg/batch {:.04f}s; avg/tile {:.04f}s'.format(pt2,pt2/batchCount,pt2/(batchCount*8)))
	print('Time: Total {:.04f}s; avg/batch {:.04f}s; avg/tile {:.04f}s'.format(pt1,pt1/batchCount,pt1/(batchCount*8)))

#test online scaled batch
# NOT YET IMPLEMENTED
if False:
	factors=[upRes**r for r in range(1, recursionDepth+1)] # can be arbitrary, 1 is always included
	print('Output online scaled batch. factors: {}'.format(factors))
	batch_scaled = TC.selectRandomTilesRecScale(selectionSize = 4, factors=[upRes**r for r in range(1, recursionDepth+1)], augment=augment, isTraining=True)
	for r in range(recursionDepth +1):
		print('batch {}: {}'.format(r, batch_scaled[r].shape))
		imgSz = batch_scaled[r].shape[2]#int((batch_scaled[r].shape[1]//4)**(1.0/2) + 0.5)
		tileShape = (batch_scaled[r].shape[0], imgSz,imgSz,3)
		tile = np.reshape(batch_scaled[r], tileShape)
		print('rec {}: tile shape: {}'.format(r, tile.shape))
		tc.savePngsGrayscale(tile, test_path + 'rec_{}_'.format(r),imageCounter=0, tiles_in_image=[1,1], channels=[0], save_rgb = [[1,2]], rgb_interval=[-2,2], plot_vel_x_y=False)

print('')
if fail:
	print('--- TEST FINISHED DESPITE INDICATED FAILURE---')
else:
	print('--- TEST FINISHED ---')
print('')
