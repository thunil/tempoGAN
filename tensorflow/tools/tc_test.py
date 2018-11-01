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

import os,sys
from itertools import repeat
import tilecreator_t as tc
import fluiddataloader as fdl
import numpy as np
import paramhelpers as ph


out_path		=	 ph.getParam( "basePath",		'../test_out/' )
sim_path		=	 ph.getParam( "basePath",		'../data_sim/' )
randSeed		= int(ph.getParam( "randSeed",		1 )) 				# seed for np and tf initialization

simSize  	= int(ph.getParam( "simSize", 		256 )) 			# tiles of low res sim
tileSize 	= int(ph.getParam( "tileSize", 		64 )) 			# size of low res tiles
upRes	  		= int(ph.getParam( "upRes", 		4 )) 				# single generator scaling factor
dim   = int(ph.getParam( "dim",		 2 )) 				# dimension of dataset

augment = int(ph.getParam( "aug", 1 ))		 # use dataAugmentation or not

# no high res data in TC, using high data in TC's low res
useScaledData = int(ph.getParam( "scaled", 0 ))
useLabelData = int(ph.getParam( "label", 0 ))
useDataBlocks = int(ph.getParam( "block", 0 ))
blockSize = int(ph.getParam( "blockSize", 1 ))

batchCount = int(ph.getParam( "batchCount", 1 ))


ph.checkUnusedParams()
np.random.seed(randSeed)
#tf.set_random_seed(randSeed)


#if not os.path.exists(out_path):
#	os.makedirs(out_path)
test_path,_ = ph.getNextTestPath(0, out_path)
sys.stdout = ph.Logger(test_path)
sys.stderr = ph.ErrorLogger(test_path)

print('')
print('--- TEST STARTED ---')
print('')

print("\nUsing parameters:\n"+ph.paramsToString())

recursionDepth = 0

fromSim = 1018
toSim = fromSim

dirIDs = np.linspace(fromSim, toSim, (toSim-fromSim+1),dtype='int16')
lowfilename = "density_low_%04d.npz"
highfilename = "density_high_%04d.npz"


if not useScaledData:
	lowfilename = highfilename
	highfilename = None
else:
	simSize = simSize//upRes
	tileSize = tileSize//upRes

#load data
mfl = ["density", "velocity"] if not useDataBlocks else ["density", "velocity", "density", "velocity", "density", "velocity" ]
mol  = [0,0] if not useDataBlocks else [0,0,1,1,2,2]
mfh = None
moh = None
if useScaledData:
	mfh = ["density", "velocity"] if not useDataBlocks else ["density", "velocity", "density", "velocity", "density", "velocity" ]
	moh  = [0,0] if not useDataBlocks else [0,0,1,1,2,2]

floader = fdl.FluidDataLoader( print_info=1, base_path=sim_path, filename=lowfilename, oldNamingScheme=False, filename_y=highfilename, filename_index_max=200, indices=dirIDs, data_fraction=0.5, multi_file_list=mfl, multi_file_idxOff=mol, multi_file_list_y=mfh , multi_file_idxOff_y=moh)
x, y, xFilenames  = floader.get()

tile_format='NYXC'
z_axis = 1
if useDataBlocks:
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
if False:
	print('Output reference')
	tileShape = (x.shape[0],simSize,simSize,x.shape[-1])
	tiles = np.reshape(x, tileShape)
	tc.savePngs(tiles[:1], test_path + 'ref_low_',imageCounter=0, tiles_in_image=[1,1], plot_vel_x_y=False, channels=[0], save_rgb = [[2,3]], rgb_interval=[-2,2])
	if useScaledData:
		tileShape = (y.shape[0],simSize*upRes,simSize*upRes,y.shape[-1])
		tiles = np.reshape(y, tileShape)
		tc.savePngs(tiles[:1], test_path + 'ref_high_',imageCounter=0, tiles_in_image=[1,1], plot_vel_x_y=False, channels=[0], save_rgb = [[2,3]], rgb_interval=[-2,2])

# tilecreator
TC = tc.TileCreator(tileSize=tileSize, simSize=simSize , dim=dim, densityMinimum=0.1, scaleFactor=upRes, channelLayout_main='d,vx,vy', channelLayout_scaled='d,vx,vy', useScaledData=useScaledData, useLabels=useLabelData, useDataBlocks=useDataBlocks, logLevel=10)
if augment:
	TC.initDataAugmentation(2)

# strip zero z vel of 2D data
if dim==2:
	x,_ = np.split(x, [3], axis=-1)
	if useScaledData:
		y,_ = np.split(y, [3], axis=-1)

# add low data with dummy labels
TC.addData(x, y if useScaledData else None, l, b)

#bx,by = TC.selectRandomTiles(64, True, augment=True)



#test batch:
if True:
	imageCounter=0
	for batch_number in range(batchCount):
		i=0
		print('\nOutput batch {}'.format(batch_number))
		batch = TC.selectRandomTiles(selectionSize = 8, augment=augment, isTraining=True, blockSize=blockSize, squeezeZ=True)
		print('batch_x shape: {}'.format(batch[0].shape))
		
		#tileShape = (batch[0].shape[0],tileSize,tileSize,batch[0].shape[-1])
		tiles = batch[i]
		i+=1
		#print('tiles_x shape: {}'.format(tiles.shape))
		ic=tc.savePngs(tiles, test_path, tile_format=tile_format,imageCounter=imageCounter, tiles_in_image=[1,1], plot_vel_x_y=False, channels=[0], save_rgb = [[1,2]], rgb_interval=[-2,2])
		if useScaledData:
			print('batch_y shape: {}'.format(batch[1].shape))
			#tileShape = (batch[1].shape[0],tileSize*upRes,tileSize*upRes,batch[1].shape[-1])
			tiles = batch[i]
			i+=1
			#print('tiles_y shape: {}'.format(tiles.shape))
			tc.savePngs(tiles, test_path + 'high_', tile_format=tile_format, imageCounter=imageCounter, tiles_in_image=[1,1], plot_vel_x_y=False, channels=[0], save_rgb = [[1,2]], rgb_interval=[-2,2])
		if useLabelData:
			tiles = batch[i]
			i+=1
			print('labels: {}'.format(tiles))
		print('-> images {} to {}'.format(imageCounter, ic-1))
		imageCounter=ic

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
print('--- TEST FINISHED ---')
print('')
