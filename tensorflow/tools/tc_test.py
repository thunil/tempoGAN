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
import tilecreator_t as tc
import fluiddataloader as fdl
import numpy as np

upRes=4
simSize=256
tileSize=64
recursionDepth = 0

sim_path = '../data_sim/'
out_path = '../test_out/'

fromSim = 1018
toSim = fromSim

dirIDs = np.linspace(fromSim, toSim, (toSim-fromSim+1),dtype='int16')
lowfilename = "density_low_%04d.npz"
highfilename = "density_high_%04d.npz"

# no high res data in TC, using high data in TC's low res
highIsLabel = True
if highIsLabel:
	lowfilename = highfilename
	highfilename = None
else:
	simSize = simSize//upRes
	tileSize = tileSize//upRes

#load data
mfl = ["density", "velocity"]#, "density", "velocity", "density", "velocity" ]
mol  = [0,0]
mfh = None if highIsLabel else ["density", "velocity"]
moh  = None if highIsLabel else [0,0]

floader = fdl.FluidDataLoader( print_info=1, base_path=sim_path, filename=lowfilename, oldNamingScheme=False, filename_y=highfilename, filename_index_max=200, indices=dirIDs, data_fraction=0.5, multi_file_list=mfl, multi_file_idxOff=mol, multi_file_list_y=mfh , multi_file_idxOff_y=moh)
x, y, xFilenames  = floader.get()

print(x.shape)
if not highIsLabel: print(y.shape)
#save ref:
if False:
	print('Output reference')
	tileShape = (x.shape[0],simSize,simSize,x.shape[-1])
	tiles = np.reshape(x, tileShape)
	tc.savePngsGrayscale(tiles[:1], out_path + 'ref_low_',imageCounter=0, tiles_in_image=[1,1], plot_vel_x_y=False, channels=[0], save_rgb = [[2,3]], rgb_interval=[-2,2])
	if not highIsLabel:
		tileShape = (y.shape[0],simSize*upRes,simSize*upRes,y.shape[-1])
		tiles = np.reshape(y, tileShape)
		tc.savePngsGrayscale(tiles[:1], out_path + 'ref_high_',imageCounter=0, tiles_in_image=[1,1], plot_vel_x_y=False, channels=[0], save_rgb = [[2,3]], rgb_interval=[-2,2])

# tilecreator
TC = tc.TileCreator(tileSizeLow=tileSize, simSizeLow=simSize , dim =2, dim_t = 1,densityMinimum=0.1, upres=upRes, channelLayout_low='d,vx,vy', channelLayout_high='d,vx,vy', premadeTiles=False, highIsLabel=highIsLabel)
TC.initDataAugmentation(2)

# strip zero z vel of 2D data
x,_ = np.split(x, [3], axis=-1)
if not highIsLabel: y,_ = np.split(y, [3], axis=-1)

# add low data with dummy labels
TC.addData(x,np.zeros(x.shape[0]) if highIsLabel else y)

#bx,by = TC.selectRandomTiles(64, True, augment=True)

if not os.path.exists(out_path):
	os.makedirs(out_path)


#test batch:
if True:
	print('Output normal batch')
	batch_x, batch_y = TC.selectRandomTiles(selectionSize = 8, augment=True, isTraining=True)
	print('batch_x shape: {}'.format(batch_x.shape))
	print('batch_y shape: {}'.format(batch_y.shape))
	
	tileShape = (batch_x.shape[0],tileSize,tileSize,batch_x.shape[-1])
	tiles = np.reshape(batch_x, tileShape)
	tc.savePngsGrayscale(tiles, out_path,imageCounter=0, tiles_in_image=[1,1], plot_vel_x_y=False, channels=[0], save_rgb = [[1,2]], rgb_interval=[-2,2])
	if not highIsLabel:
		tileShape = (batch_y.shape[0],tileSize*upRes,tileSize*upRes,batch_y.shape[-1])
		tiles = np.reshape(batch_y, tileShape)
		tc.savePngsGrayscale(tiles, out_path + 'high_',imageCounter=0, tiles_in_image=[1,1], plot_vel_x_y=False, channels=[0], save_rgb = [[1,2]], rgb_interval=[-2,2])

#test online scaled batch
# NOT YET IMPLEMENTED
if False:
	factors=[upRes**r for r in range(1, recursionDepth+1)] # can be arbitrary, 1 is always included
	print('Output online scaled batch. factors: {}'.format(factors))
	batch_scaled = TC.selectRandomTilesRecScale(selectionSize = 4, factors=[upRes**r for r in range(1, recursionDepth+1)], augment=True, isTraining=True)
	for r in range(recursionDepth +1):
		print('batch {}: {}'.format(r, batch_scaled[r].shape))
		imgSz = batch_scaled[r].shape[2]#int((batch_scaled[r].shape[1]//4)**(1.0/2) + 0.5)
		tileShape = (batch_scaled[r].shape[0], imgSz,imgSz,3)
		tile = np.reshape(batch_scaled[r], tileShape)
		print('rec {}: tile shape: {}'.format(r, tile.shape))
		tc.savePngsGrayscale(tile, out_path + 'rec_{}_'.format(r),imageCounter=0, tiles_in_image=[1,1], channels=[0], save_rgb = [[1,2]], rgb_interval=[-2,2], plot_vel_x_y=False)


