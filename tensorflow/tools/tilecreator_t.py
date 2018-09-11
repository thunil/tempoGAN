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

import os, math
import shutil, sys
from random import seed, random, randrange
import uniio
import numpy as np
import scipy.misc
import scipy.ndimage 
import imageio

# check whether matplotlib is available to generate vector/quiver plots
import imp
try:
	imp.find_module('matplotlib')
	import matplotlib.pyplot
	found_matplotlib = True
except ImportError:
    found_matplotlib = False
#import matplotlib.pyplot as plt

# global channel keys, have to be one char
C_KEY_DEFAULT = 'd'
C_KEY_VELOCITY = 'v'
C_KEY_VORTICITY = 'x'
C_KEY_POSITION = 'p'

DATA_KEY_LOW = 0
DATA_KEY_HIGH= 1

#keys for augmentation operations
AOPS_KEY_ROTATE = 'rot'
AOPS_KEY_SCALE = 'scale'
AOPS_KEY_ROT90 = 'rot90'
AOPS_KEY_FLIP = 'flip'

seed( 42 )

# default channel layouts
C_LAYOUT = {
	'dens':C_KEY_DEFAULT,
	'dens_vel':'d,vx,vy,vz'
	}

class TileCreator(object):

	def __init__(self, tileSizeLow, simSizeLow=64, upres=2, dim=2, dim_t=1, overlapping=0, densityMinimum=0.02, premadeTiles=False, partTrain=0.8, partTest=0.2, partVal=0, channelLayout_low=C_LAYOUT['dens_vel'], channelLayout_high=C_LAYOUT['dens'], highIsLabel=False, loadPN=False, padding=0):
		'''
			tileSizeLow, simSizeLow: int, [int,int] if 2D, [int,int,int]
			channelLayout: 'key,key,...'
				the keys are NOT case sensitive and leading and trailing whitespace characters are REMOVED.
				key:
					default:	d
					velocity:	v[label](x|y|z)
						label can be arbitrary or empty,
						key must be unique and x,y must exist while z is optional in 2D, x,y,z must exist in 3D.
						if x does not exist y,z will be ignored (treaded as 'd').
					rest is not yet supported
				
			premadeTiles: cut regular tiles when loading data, can't use data augmentation
			part(Train|Test|Val): relative size of the different data sets
			highIsLabel: high data is not augmented
			loadHigh: 
			simPath: path to the uni simulation files
			loadPath: packed simulations are stored here
		'''
		
		# DATA DIMENSION
		self.dim_t = dim_t # same for hi_res or low_res
		if dim!=2 and dim!=3:
			self.TCError('Data dimension must be 2 or 3.')
		self.dim = dim
		# TILE SIZE
		if np.isscalar(tileSizeLow):
			self.tileSizeLow = [tileSizeLow, tileSizeLow, tileSizeLow]
		elif len(tileSizeLow)==2 and self.dim==2:
			self.tileSizeLow = [1]+tileSizeLow
		elif len(tileSizeLow)==3:
			self.tileSizeLow = tileSizeLow
		else:
			self.TCError('Tile size mismatch.')
		self.tileSizeLow = np.asarray(self.tileSizeLow)
		#SIM SIZE
		if np.isscalar(simSizeLow):
			self.simSizeLow = [simSizeLow, simSizeLow, simSizeLow]
		elif len(simSizeLow)==2 and self.dim==2:
			self.simSizeLow = [1]+simSizeLow
		elif len(simSizeLow)==3:
			self.simSizeLow = simSizeLow
		else:
			self.TCError('Simulation size mismatch.')
		self.simSizeLow = np.asarray(self.simSizeLow)
		
		if upres < 1:
			self.TCError('Upres must be at least 1.')
		self.upres = upres
		if not highIsLabel:
			self.tileSizeHigh = self.tileSizeLow*upres
			self.simSizeHigh = self.simSizeLow*upres
		else:
			self.tileSizeHigh = np.asarray([1])
			self.simSizeHigh = np.asarray([1])
		
		if self.dim==2:
			self.tileSizeLow[0]=1
			self.tileSizeHigh[0]=1
			self.simSizeLow[0]=1
			self.simSizeHigh[0]=1
		if np.less(self.simSizeLow, self.tileSizeLow).any():
			self.TCError('Tile size {} can not be larger than sim size {}, {}.'.format(self.tileSizeLow,self.simSizeLow))
		
		
		if densityMinimum<0.:
			self.TCError('densityMinimum can not be negative.')
		self.densityMinimum = densityMinimum
		self.premadeTiles = premadeTiles
		self.useDataAug = False
		
		#CHANNELS
		self.c_lists = {}
		self.c_low, self.c_lists[DATA_KEY_LOW] = self.parseChannels(channelLayout_low)
		self.c_high, self.c_lists[DATA_KEY_HIGH] = self.parseChannels(channelLayout_high)

		# print info
		print('\n')
		print('Dimension: {}, time dimension: {}'.format(self.dim,self.dim_t))
		print('Low-res data:')
		print('  channel layout: {}'.format(self.c_low))
		print('  default channels: {}'.format(self.c_lists[DATA_KEY_LOW][C_KEY_DEFAULT]))
		if len(self.c_lists[DATA_KEY_LOW][C_KEY_VELOCITY])>0: 
			print('  velocity channels: {}'.format(self.c_lists[DATA_KEY_LOW][C_KEY_VELOCITY]))
		if len(self.c_lists[DATA_KEY_LOW][C_KEY_VORTICITY])>0: 
			print('  vorticity channels: {}'.format(self.c_lists[DATA_KEY_LOW][C_KEY_VORTICITY]))
		print('High-res data:')
		if highIsLabel:
			print('  is Label')
		print('  channel layout: {}'.format(self.c_high))
		print('  default channels: {}'.format(self.c_lists[DATA_KEY_HIGH][C_KEY_DEFAULT]))
		if len(self.c_lists[DATA_KEY_HIGH][C_KEY_VELOCITY])>0: 
			print('  velocity channels: {}'.format(self.c_lists[DATA_KEY_HIGH][C_KEY_VELOCITY]))
		if len(self.c_lists[DATA_KEY_HIGH][C_KEY_VORTICITY])>0: 
			print('  vorticity channels: {}'.format(self.c_lists[DATA_KEY_HIGH][C_KEY_VORTICITY]))
		#self.channels=len(self.c)
		
		self.data_flags = {
			DATA_KEY_LOW:{
				'isLabel':False,
				'channels':len(self.c_low),
				C_KEY_VELOCITY:len(self.c_lists[DATA_KEY_LOW][C_KEY_VELOCITY])>0,
				C_KEY_VORTICITY:len(self.c_lists[DATA_KEY_LOW][C_KEY_VORTICITY])>0,
				C_KEY_POSITION:False
			},
			DATA_KEY_HIGH:{
				'isLabel':highIsLabel,
				'channels':len(self.c_high),
				C_KEY_VELOCITY:len(self.c_lists[DATA_KEY_HIGH][C_KEY_VELOCITY])>0,
				C_KEY_VORTICITY:len(self.c_lists[DATA_KEY_HIGH][C_KEY_VORTICITY])>0,
				C_KEY_POSITION:False
			}
		}
		if loadPN:
			self.TCError('prev and next tiles not supported.')
		self.hasPN = loadPN
		self.padding=padding
		
		#if self.hasPN:
		#[z,y,x, velocities an/or position if enabled (density,vel,vel,vel, pos, pos [,pos])]
		
		#DATA SHAPES
		self.tile_shape_low = np.append(self.tileSizeLow, [self.data_flags[DATA_KEY_LOW]['channels']])
		self.frame_shape_low = np.append(self.simSizeLow, [self.data_flags[DATA_KEY_LOW]['channels']])
		if not self.data_flags[DATA_KEY_HIGH]['isLabel']:
			self.tile_shape_high = np.append(self.tileSizeHigh, [self.data_flags[DATA_KEY_HIGH]['channels']])
			self.frame_shape_high = np.append(self.simSizeHigh, [self.data_flags[DATA_KEY_HIGH]['channels']])
		else:
			self.tile_shape_high = self.tileSizeHigh[:]
			self.frame_shape_high = self.simSizeHigh[:]
		
		self.densityThreshold = (self.densityMinimum * self.tile_shape_low[0] * self.tile_shape_low[1] * self.tile_shape_low[2])

		self.data = {
			DATA_KEY_LOW:[],
			DATA_KEY_HIGH:[]
		}
		
		all=partTrain+partTest+partVal
		self.part_train=partTrain/all
		self.part_test=partTest/all
		self.part_validation=partVal/all
	
	def initDataAugmentation(self, rot=2, minScale=0.85, maxScale=1.15 ,flip=True):
		'''
			set up data augmentation
			rot: 1: 90 degree rotations; 2: full rotation; else: nop rotation
			Scale: if both 1 disable scaling
		'''
		self.useDataAug = True
		
		"""
			specify the special augmentation operation needed for some channel types here
			will only be applyed if the specified channel type is in the data
			** Tempo Datum may have multiple channels as coherent frames, [batch, z, y, x, t*channels]
			** They are reshaped first before these aops, [batch, z, y, x, t, channels], and then reshape back
			** Because of this extra time dimention, all aops can only do isolate calculations, for e.g., value scaling,
			** Any calculation relay on neighborhood will be wrong, for e.g., spacial scaling (zoom).
		"""
		self.aops = {
			DATA_KEY_LOW:{
				AOPS_KEY_ROTATE:{
					C_KEY_VELOCITY:self.rotateVelocities,
					C_KEY_VORTICITY:self.rotateVelocities
				},
				AOPS_KEY_SCALE:{
					C_KEY_VELOCITY:self.scaleVelocities,
					C_KEY_VORTICITY:self.scaleVelocities
				},
				AOPS_KEY_ROT90:{
					C_KEY_VELOCITY:self.rotate90Velocities,
					C_KEY_VORTICITY:self.rotate90Velocities
				},
				AOPS_KEY_FLIP:{
					C_KEY_VELOCITY:self.flipVelocities,
					C_KEY_VORTICITY:self.flipVelocities
				}
				
			},
			DATA_KEY_HIGH:{
				AOPS_KEY_ROTATE:{
					C_KEY_VELOCITY:self.rotateVelocities,
					C_KEY_VORTICITY:self.rotateVelocities
				},
				AOPS_KEY_SCALE:{
					C_KEY_VELOCITY:self.scaleVelocities,
					C_KEY_VORTICITY:self.scaleVelocities
				},
				AOPS_KEY_ROT90:{
					C_KEY_VELOCITY:self.rotate90Velocities,
					C_KEY_VORTICITY:self.rotate90Velocities
				},
				AOPS_KEY_FLIP:{
					C_KEY_VELOCITY:self.flipVelocities,
					C_KEY_VORTICITY:self.flipVelocities
				}
			}
		}

		msg = 'data augmentation: '

		if rot==2:
			self.do_rotation = True
			self.do_rot90 = False
			msg += 'rotation, '
		elif rot==1:
			self.do_rotation = False
			self.do_rot90 = True
			msg += 'rot90, '
			z=(2,1)
			nz=(1,2)
			x=(1,0)
			y=(0,2)
			nx=(0,1)
			ny=(2,0)
			# thanks to http://www.euclideanspace.com/maths/discrete/groups/categorise/finite/cube/
			self.cube_rot = {2: [[],[z],[z,z],[nz]], 3: [[],[x],[y],[x,x],[x,y],[y,x],[y,y],[nx],[x,x,y],[x,y,x],[x,y,y],[y,x,x],[y,y,x],[ny],[nx,y],[x,x,y,x],[x,x,y,y],[x,y,x,x],[x,ny],[y,nx],[ny,x],[nx,y,x],[x,y,nx],[x,ny,x]]}
		else:
			self.do_rotation = False
			self.do_rot90 = False
			
		self.scaleFactor = [minScale, maxScale]
		if (self.scaleFactor[0]==1 and self.scaleFactor[1]==1):
			self.do_scaling = False
		else:
			self.do_scaling = True
			msg += 'scaling, '
			
		self.do_flip = flip
		if self.do_flip:
			msg += 'flip'
		msg += '.'
		print(msg)
		self.interpolation_order = 1
		self.fill_mode = 'constant'
		
	
	
	def addData(self, low, high):
		'''
			add data, tiles if premadeTiles, frames otherwise.
			low, high: list of or single 3D data np arrays
		'''
		# check data shape
		low = np.asarray(low)
		high = np.asarray(high)
		
		if not self.data_flags[DATA_KEY_HIGH]['isLabel']:
			if len(low.shape)!=len(high.shape): #high-low mismatch
				self.TCError('Data shape mismatch. Dimensions: {} low vs {} high. Dimensions must match or use highIsLabel.'.format(len(low.shape),len(high.shape)) )
		if not (len(low.shape)==4 or len(low.shape)==5): #not single frame or sequence of frames
			self.TCError('Input must be single 3D data or sequence of 3D data. Format: ([batch,] z, y, x, channels). For 2D use z=1.')

		if (low.shape[-1]!=(self.dim_t * self.data_flags[DATA_KEY_LOW]['channels'])):
			self.TCError('Dim_t ({}) * Channels ({}, {}) configured for LOW-res data don\'t match channels ({}) of input data.'.format(self.dim_t, self.data_flags[DATA_KEY_LOW]['channels'], self.c_low,  low.shape[-1]) )
		if not self.data_flags[DATA_KEY_HIGH]['isLabel']:
			if (high.shape[-1]!=(self.dim_t * self.data_flags[DATA_KEY_HIGH]['channels'])):
				self.TCError('Dim_t ({}) * Channels ({}, {}) configured for HIGH-res data don\'t match channels ({}) of input data.'.format(self.dim_t, self.data_flags[DATA_KEY_HIGH]['channels'], self.c_high, high.shape[-1]) )
		
		low_shape = low.shape
		high_shape = high.shape
		if len(low.shape)==5: #sequence
			if low.shape[0]!=high.shape[0]: #check amount
				self.TCError('Unequal amount of low ({}) and high ({}) data.'.format(low.shape[1], high.shape[1]))
			# get single data shape
			low_shape = low_shape[1:]
			if not self.data_flags[DATA_KEY_HIGH]['isLabel']:
				high_shape = high_shape[1:]
			else: high_shape = [1]
		else: #single
			low = [low]
			high = [high]
		
		if self.premadeTiles:
			if not (self.dim_t == 1):
				self.TCError('Currently, Dim_t = {} > 1 is not supported by premade tiles'.format(self.dim_t))
			if not np.array_equal(low_shape, self.tile_shape_low) or not np.array_equal(high_shape,self.tile_shape_high):
				self.TCError('Tile shape mismatch: is - specified\n\tlow: {} - {}\n\thigh {} - {}'.format(low_shape, self.tile_shape_low, high_shape,self.tile_shape_high))
		else:
			single_frame_low_shape = list(low_shape)
			single_frame_high_shape = list(high_shape)
			single_frame_low_shape[-1] = low_shape[-1] // self.dim_t
			if not self.data_flags[DATA_KEY_HIGH]['isLabel']:
				single_frame_high_shape[-1] = high_shape[-1] // self.dim_t
			
			if not np.array_equal(single_frame_low_shape, self.frame_shape_low) or not np.array_equal(single_frame_high_shape,self.frame_shape_high):
				self.TCError('Frame shape mismatch: is - specified\n\tlow: {} - {}\n\thigh: {} - {}, given dim_t as {}'.format(single_frame_low_shape, self.frame_shape_low, single_frame_high_shape,self.frame_shape_high, self.dim_t))

		self.data[DATA_KEY_LOW].extend(low)
		self.data[DATA_KEY_HIGH].extend(high)
		
		print('\n')
		print('Added {} datasets. Total: {}'.format(low.shape[0], len(self.data[DATA_KEY_LOW])))
		self.splitSets()
	
	def splitSets(self):
		'''
			calculate the set borders for training, testing and validation set
		'''
		length = len(self.data[DATA_KEY_LOW])
		
		end_train = int( length * self.part_train )
		end_test = end_train + int( length * self.part_test )
		#just store the borders of the different sets to avoid data duplication
		self.setBorders = [end_train, end_test, length]
		
		print('Training set: {}'.format(self.setBorders[0]))
		print('Testing set:  {}'.format(self.setBorders[1]-self.setBorders[0]))
		print('Validation set:  {}'.format(self.setBorders[2]-self.setBorders[1]))
		
	def clearData(self):
		'''
			clears the data buffer
		'''
		self.data = {
			DATA_KEY_LOW:[],
			DATA_KEY_HIGH:[]
		}
	
	def createTiles(self, data, tileShape, strides=-1): 
		'''
			create tiles from a single frame. fixed, regular pattern
			strides: <=0 or tileShape is normal, otherwise create overlapping tiles
		'''
		dataShape = data.shape #2D sim: [1,res,res,channels]
		pad = [self.padding,self.padding,self.padding,0]
		if np.isscalar(strides):
			if strides <= 0:
				strides = tileShape
			else:
				strides = [strides,strides,strides]
		if dataShape[0]<=1:
			pad[0] = 0
			strides[0] = 1
		channels = dataShape[3]
		noTiles = [ (dataShape[0]-tileShape[0])//strides[0]+1, (dataShape[1]-tileShape[1])//strides[1]+1, (dataShape[2]-tileShape[2])//strides[2]+1 ]
		tiles = []

		for tileZ in range(0, noTiles[0]):
			for tileY in range(0, noTiles[1]):
				for tileX in range(0, noTiles[2]):
					idx_from=[tileZ*strides[0], tileY*strides[1], tileX*strides[2]]
					idx_to=[idx_from[0]+tileShape[0], idx_from[1]+tileShape[1], idx_from[2]+tileShape[2]]
					currTile=data[ idx_from[0]:idx_to[0], idx_from[1]:idx_to[1], idx_from[2]:idx_to[2], :]
					if self.padding > 0:
						currTile = np.pad(currTile, pad, 'edge')

					tiles.append(currTile)
		return np.array(tiles)
		
	def cutTile(self, data, tileShape, offset=[0,0,0]): 
		'''
			cut a tile of with shape and offset 
		'''
		offset = np.asarray(offset)
		tileShape = np.asarray(tileShape)
		tileShape[-1] = data.shape[-1]
		if np.less(data.shape[:3], tileShape[:3]+offset[:3]).any():
			self.TCError('Can\'t cut tile with shape {} and offset{} from data with shape {}.'.format(tileShape, offset, data.shape))
		
		tile = data[offset[0]:offset[0]+tileShape[0], offset[1]:offset[1]+tileShape[1], offset[2]:offset[2]+tileShape[2], :]
		
		if not np.array_equal(tile.shape,tileShape):
			self.TCError('Wrong tile shape after cutting. is: {}. goal: {}.'.format(tile.shape,tileShape))
		return tile
		
#####################################################################################
# batch creation
#####################################################################################
	
	
	def selectRandomTiles(self, selectionSize, isTraining=True, augment=False, tile_t = 1):
		'''
			main method to create baches
			Return:
				shape: [selectionSize, z, y, x, channels * tile_t]
				if 2D z = 1
				channels: density, [vel x, vel y, vel z], [pos x, pox y, pos z]
		'''
		if isTraining:
			if self.setBorders[0]<1:
				self.TCError('no training data.')
		else:
			if (self.setBorders[1] - self.setBorders[0])<1:
				self.TCError('no test data.')
		if(tile_t > self.dim_t):
			self.TCError('not enough coherent frames. Requested {}, available {}'.format(tile_t, self.dim_t))
		batch_low = []
		batch_high = []
		for i in range(selectionSize):
			if augment and self.useDataAug: #data augmentation
				low, high = self.generateTile(isTraining, tile_t)
			else: #cut random tile without augmentation
				low, high = self.getRandomDatum(isTraining, tile_t)
				if not self.premadeTiles:
					low, high = self.getRandomTile(low, high)
			batch_low.append(low)
			batch_high.append(high)
			
		return np.asarray(batch_low), np.asarray(batch_high)
		
	def generateTile(self, isTraining=True, tile_t = 1):
		'''
			generates a random low-high pair of tiles (data augmentation)
		'''
		# get a frame, is a copy to avoid transormations affecting the raw dataset
		data = {}
		data[DATA_KEY_LOW], data[DATA_KEY_HIGH] = self.getRandomDatum(isTraining, tile_t)
		
		if not self.premadeTiles:
			#cut a tile for faster transformation
			if self.do_scaling or self.do_rotation:
				factor = 1
				if self.do_rotation: # or self.do_scaling:
					factor*=1.5 # scaling: to avoid size errors caused by rounding
				if self.do_scaling:
					scaleFactor = np.random.uniform(self.scaleFactor[0], self.scaleFactor[1])
					factor/= scaleFactor 
				tileShapeLow = np.ceil(self.tile_shape_low*factor)
				if self.dim==2:
					tileShapeLow[0] = 1
				data[DATA_KEY_LOW], data[DATA_KEY_HIGH] = self.getRandomTile(data[DATA_KEY_LOW], data[DATA_KEY_HIGH], tileShapeLow.astype(int))
		
		
			#random scaling, changes resolution
			if self.do_scaling:
				data = self.scale(data, scaleFactor)
		
		
			bounds = np.zeros(4)
		
			#rotate
			if self.do_rotation:
				bounds = np.array(data[DATA_KEY_LOW].shape)*0.16 #bounds applied on all sides, 1.5*(1-2*0.16)~1
				data = self.rotate(data)
		
			#get a tile
			data[DATA_KEY_LOW], data[DATA_KEY_HIGH] = self.getRandomTile(data[DATA_KEY_LOW], data[DATA_KEY_HIGH], bounds=bounds) #includes "shifting"
		
		if self.do_rot90:
			rot = np.random.choice(self.cube_rot[self.dim])
			for axis in rot:
				data = self.rotate90(data, axis)
			
		#flip once
		if self.do_flip:
			axis = np.random.choice(4)
			if axis < 3: # axis < self.dim
				data = self.flip(data, [axis])
		
		# check tile size
		target_shape_low = np.copy(self.tile_shape_low)
		target_shape_high = np.copy(self.tile_shape_high)
		target_shape_low[-1] *= tile_t
		target_shape_high[-1] *= tile_t
		
		if not np.array_equal(data[DATA_KEY_LOW].shape,target_shape_low) or (not np.array_equal(data[DATA_KEY_HIGH].shape,target_shape_high) and not self.data_flags[DATA_KEY_HIGH]['isLabel']):
			self.TCError('Wrong tile shape after data augmentation. is: {},{}. goal: {},{}.'.format(data[DATA_KEY_LOW].shape, data[DATA_KEY_HIGH].shape, target_shape_low, target_shape_high))
		
		return data[DATA_KEY_LOW], data[DATA_KEY_HIGH]
	
	def getRandomDatum(self, isTraining=True, tile_t = 1):
		'''returns a copy of a random frame'''
		if isTraining:
			randNo = randrange(0, self.setBorders[0])
		else:
			randNo = randrange(self.setBorders[0], self.setBorders[1])
		randFrame = 0
		if tile_t<self.dim_t:
			randFrame = randrange(0, self.dim_t - tile_t)
		else:
			tile_t = self.dim_t

		return self.getDatum(randNo*self.dim_t+randFrame, tile_t)

	def getDatum(self, index, tile_t = 1):
		'''returns a copy of the indicated frame or tile'''
		begin_ch = 0
		if(self.dim_t > 1):
			begin_ch = (index % self.dim_t) * self.tile_shape_low[-1]
		end_ch = begin_ch + tile_t * self.tile_shape_low[-1]
		begin_ch_y = 0
		if(self.dim_t > 1):
			begin_ch_y = (index % self.dim_t) * self.tile_shape_high[-1]
		end_c_h_y = begin_ch_y + tile_t * self.tile_shape_high[-1]
		
		if not self.data_flags[DATA_KEY_HIGH]['isLabel']:
			return np.copy(self.data[DATA_KEY_LOW][index//self.dim_t][:,:,:,begin_ch:end_ch]), np.copy(self.data[DATA_KEY_HIGH][index//self.dim_t][:,:,:,begin_ch_y:end_c_h_y])
		else:
			return np.copy(self.data[DATA_KEY_LOW][index//self.dim_t][:,:,:,begin_ch:end_ch]), np.copy(self.data[DATA_KEY_HIGH][index//self.dim_t])


	def getRandomTile(self, low, high, tileShapeLow=None, bounds=[0,0,0,0]): #bounds to avoid mirrored parts
		'''
			cut a random tile (low and high) from a given frame, considers densityMinimum
			bounds: ignore edges of frames, used to discard mirrored parts after rotation
		'''
		
		if tileShapeLow is None:
			tileShapeLow = np.copy(self.tile_shape_low) # use copy is very important!!!
		tileShapeHigh = tileShapeLow*self.upres
		
		frameShapeLow = np.asarray(low.shape)
		if len(low.shape)!=4 or len(tileShapeLow)!=4:
			self.TCError('Data shape mismatch.')
		if  len(high.shape)!=4 and not self.data_flags[DATA_KEY_HIGH]['isLabel']:
			self.TCError('Data shape mismatch.')

		start = np.ceil(bounds)
		end = frameShapeLow - tileShapeLow + np.ones(4) - start
		
		offset_up = np.array([self.upres, self.upres, self.upres])
		
		if self.dim==2:
			start[0] = 0
			end[0] = 1
			offset_up[0] = 1
			tileShapeHigh[0] = 1
			
		# check if possible to cut tile
		if np.amin((end-start)[:3]) < 0:
			self.TCError('Can\'t cut tile {} from frame {} with bounds {}.'.format(tileShapeLow, frameShapeLow, start))
	
		# cut tile
		hasMinDensity = False
		i = 1
		while (not hasMinDensity) and i<20:
			offset = np.asarray([randrange(start[0], end[0]), randrange(start[1], end[1]), randrange(start[2], end[2])])
			lowTile = self.cutTile(low, tileShapeLow, offset)
			offset *= offset_up
			if not self.data_flags[DATA_KEY_HIGH]['isLabel']:
				highTile = self.cutTile(high, tileShapeHigh, offset)
			else:
				highTile = high
			hasMinDensity = self.hasMinDensity(lowTile)
			i+=1
		return lowTile, highTile

#####################################################################################
# AUGMENTATION
#####################################################################################

	def special_aug(self, data, ops_key, param):
		"""
			wrapper to call the augmentation operations specified in self.aops in initAugmentation
		"""
		for data_key in data:
			if self.data_flags[data_key]['isLabel']: continue
			orig_shape = data[data_key].shape
			tile_t = orig_shape[-1] // self.data_flags[data_key]['channels']
			data_array = data[data_key]
			if(tile_t > 1): data_array = data[data_key].reshape( (-1, tile_t, self.data_flags[data_key]['channels']) )
			for c_key, op in self.aops[data_key][ops_key].items():
				if self.data_flags[data_key][c_key]:
					data_array = op(data_array, self.c_lists[data_key][c_key], param)
			if (tile_t > 1): data[data_key] = data_array.reshape(orig_shape)
		return data
		
	def rotate(self, data):
		'''
			random uniform rotation of low and high data of a given frame
		'''
		#check if single frame
		
		#2D:
		if self.dim==2:
			theta = np.pi * np.random.uniform(0, 2)
			rotation_matrix = np.array([[1, 0, 0, 0 ],
										[0, np.cos(theta), -np.sin(theta), 0], 
										[0, np.sin(theta), np.cos(theta) , 0],
										[0, 0, 0, 1]  ])
			
		#3D:
		elif self.dim==3:
			# random uniform rotation in 3D
			quat = np.random.normal(size=4)
			quat/= np.linalg.norm(quat)
			
			q = np.outer(quat, quat)*2
			rotation_matrix = np.array([[1-q[2, 2]-q[3, 3],   q[1, 2]-q[3, 0],   q[1, 3]+q[2, 0], 0],
										[  q[1, 2]+q[3, 0], 1-q[1, 1]-q[3, 3],   q[2, 3]-q[1, 0], 0],
										[  q[1, 3]-q[2, 0],   q[2, 3]+q[1, 0], 1-q[1, 1]-q[2, 2], 0],
										[                0,                 0,                 0, 1]])
		
		data = self.special_aug(data, AOPS_KEY_ROTATE, rotation_matrix)
		
		for data_key in data:
			if not self.data_flags[data_key]['isLabel']:
				data[data_key] = self.applyTransform(data[data_key], rotation_matrix.T)
			
			
		return data
		
	def rotate_simple(self, low, high, angle):
		'''
			use a different method for rotation. about 30-40% faster than with rotation matrix, but only one axis.
		'''
		if len(low.shape)!=4 or len(high.shape)!=4:
			self.TCError('Data shape mismatch.')
		#test rot around z (axis order z,y,x,c)
		low = scipy.ndimage.rotate(low, angle, [1,2] , reshape=False, order=self.interpolation_order, mode=self.fill_mode, cval=1.0)
		high = scipy.ndimage.rotate(high, angle, [1,2] , reshape=False, order=self.interpolation_order, mode=self.fill_mode, cval=1.0)
		return low, high
		
	def rotateVelocities(self, datum, c_list, rotationMatrix):
		'''
			rotate vel vectors (channel 1-3)
		'''
		
		rotation3 = rotationMatrix[:3, :3]
		rotation2 = rotationMatrix[1:3, 1:3]
		channels = np.split(datum, datum.shape[-1], -1)
		for v in c_list:
			if len(v) == 3: # currently always ends here!! even for 2D, #z,y,x to match rotation matrix
				vel = np.stack([channels[v[2]].flatten(),channels[v[1]].flatten(),channels[v[0]].flatten()]) 
				vel = rotation3.dot(vel)
				channels[v[2]] = np.reshape(vel[0], channels[v[2]].shape)
				channels[v[1]] = np.reshape(vel[1], channels[v[1]].shape)
				channels[v[0]] = np.reshape(vel[2], channels[v[0]].shape)
			if len(v) == 2:
				vel = np.concatenate([channels[v[1]],channels[v[0]]], -1) #y,x to match rotation matrix
				shape = vel.shape
				vel = np.reshape(vel, (-1, 2))
				vel = np.reshape(rotation2.dot(vel.T).T, shape)
				vel = np.split(vel, 2, -1)
				channels[v[1]] = vel[0]
				channels[v[0]] = vel[1]
				
		return np.concatenate(channels, -1)
		
	def rotate90(self, data, axes):
		'''
			rotate the frame by 90 degrees from the first axis counterclockwise to the second
			axes: 2 int, from axis to axis; see np.rot90 
				0,1,2 -> z,y,x
		'''
		if len(axes)!=2:
			self.TCError('need 2 axes for rotate90.')
		
		for data_key in data:
			if not self.data_flags[data_key]['isLabel']:
				data[data_key] = np.rot90(data[data_key], axes=axes)
		
		data = self.special_aug(data, AOPS_KEY_ROT90, axes)
		
		return data
		
	def rotate90Velocities(self, datum, c_list, axes):
		if len(axes)!=2:
			self.TCError('need 2 axes for rotate90.')
		
		channels = np.split(datum, datum.shape[-1], -1)
		for v in c_list: #axes z,y,x -> vel x,y,z: 0,1,2 -> 2,1,0
			channels[v[-axes[0]+2]], channels[v[-axes[1]+2]] = -channels[v[-axes[1]+2]], channels[v[-axes[0]+2]]
			
		return np.concatenate(channels, -1)
	
	def flip(self, data, axes, isFrame=True): #axes=list, flip multiple at once
		'''
			flip low and high data (single frame/tile) along the specified axes
			low, high: data format: (z,x,y,c)
			axes: list of axis indices 0,1,2-> z,y,x
		'''
		# axis: 0,1,2 -> z,y,x
		
		if not isFrame:
			axes = np.asarray(axes) + np.ones(axes.shape)
		
		#flip tiles/frames
		for axis in axes:
			for data_key in data:
				if not self.data_flags[data_key]['isLabel']:
					data[data_key] = np.flip(data[data_key], axis)
			
		
		data = self.special_aug(data, AOPS_KEY_FLIP, axes)
		
		return data
	
	def flipVelocities(self, datum, c_list, axes):
		'''
			flip velocity vectors along the specified axes
			low: data with velocity to flip (4 channels: d,vx,vy,vz)
			axes: list of axis indices 0,1,2-> z,y,x
		'''
		
		# !axis order: data z,y,x
		channels = np.split(datum, datum.shape[-1], -1)
		for v in c_list: # x,y,[z], 2,1,0
			if 2 in axes: # flip vel x
				channels[v[0]] *= (-1)
			if 1 in axes:
				channels[v[1]] *= (-1)
			if 0 in axes and len(v)==3:
				channels[v[2]] *= (-1)
			
		return np.concatenate(channels, -1)
		
	
	def scale(self, data, factor):
		'''
			 changes frame resolution to "round((factor) * (original resolution))"
		'''
		# only same factor in every dim for now. how would it affect vel scaling?
		# check for 2D
		
		scale = [factor, factor, factor, 1] #single frame
		if self.dim==2:
			scale[0] = 1
		
		# to ensure high/low ration stays the same
		scale = np.round(np.array(data[DATA_KEY_LOW].shape) * scale )/np.array(data[DATA_KEY_LOW].shape)
		if len(data[DATA_KEY_LOW].shape)==5: #frame sequence
			scale = np.append([1],scale)
			
		#apply transform
		#low = self.applyTransform(low, zoom_matrix)
		#high = self.applyTransform(high, zoom_matrix)
		
		#changes the size of the frame. should work well with getRandomTile(), no bounds needed
		for data_key in data:
			if not self.data_flags[data_key]['isLabel']:
				data[data_key]  = scipy.ndimage.zoom( data[data_key], scale, order=self.interpolation_order, mode=self.fill_mode, cval=0.0)
		
		#necessary?
		data = self.special_aug(data, AOPS_KEY_SCALE, factor)
		
		return data
	
	def scaleVelocities(self, datum, c_list, factor):
		#scale vel? vel*=factor
		channels = np.split(datum, datum.shape[-1], -1)
		for v in c_list: # x,y,[z]; 2,1,0
			channels[v[0]] *= factor
			channels[v[1]] *= factor
			if len(v)==3:
				channels[v[2]] *= factor
			
		return np.concatenate(channels, -1)
		
	def applyTransform(self, data, transform_matrix, data_dim=3):
		# change axis order from z,y,x to x,y,z? (invert axis order +channel)
		if len(data.shape)!=4:
			self.TCError('Data shape mismatch.')
			
		#set transform to center; from fluiddatagenerator.py
		offset = np.array(data.shape) / 2 - np.array([0.5, 0.5, 0.5, 0])
		offset_matrix = np.array([[1, 0, 0, offset[0]], [0, 1, 0, offset[1]], [0, 0, 1, offset[2]], [0, 0, 0, 1]])
		reset_matrix  = np.array([[1, 0, 0,-offset[0]], [0, 1, 0,-offset[1]], [0, 0, 1,-offset[2]], [0, 0, 0, 1]])
		transform_matrix = np.dot(np.dot(offset_matrix, transform_matrix), reset_matrix)
		
		data = np.rollaxis(data, 3, 0) #channel to front
		channel_data = [scipy.ndimage.interpolation.affine_transform(
			channel,
			transform_matrix[:data_dim,:data_dim],
			transform_matrix[:data_dim, data_dim],
			order=self.interpolation_order,
			mode=self.fill_mode,
			cval=0.) for channel in data]
		data = np.stack(channel_data, axis=-1) # stack axis=-1 ?\
		return data

#####################################################################################
# HELPER METHODS
#####################################################################################

	def concatTiles(self, tiles, frameShape ,tileBorder=[0,0,0,0]):
		'''
			build a frame by concatenation of the given tiles.
			tiles: numpy array of same shaped tiles [batch,z,y,x,c]
			frameShape: the shape of the frame in tiles [z,y,x]
			tileBorder: cut off borders of the tiles. [z,y,x,c]
		'''
		if len(tiles.shape)!=5 or len(frameShape)!=3 or len(tileBorder)!=4:
			self.TCError('Data shape mismatch.')
		tiles_in_frame = frameShape[0]*frameShape[1]*frameShape[2]
		if tiles_in_frame != len(tiles):
			self.TCError('given tiles do not match required tiles.')
			
		# cut borders
		tileBorder = np.asarray(tileBorder)
		if np.less(np.zeros(4),tileBorder).any():
			tileShape = tiles.shape[1:] - 2*tileBorder
			tiles_cut = []
			for tile in tiles:
				tiles_cut.append(self.cutTile(tile, tileShape, tileBorder))
			tiles = tiles_cut
		
		#combine tiles to image
		frame = []
		for z in range(frameShape[0]):
			frame_slices = []
			for y in range(frameShape[1]):
				offset=z*frameShape[1]*frameShape[2] + y*frameShape[2]
				frame_slices.append(np.concatenate(tiles[offset:offset+frameShape[2]],axis=2)) #combine x
			frame.append(np.concatenate(frame_slices, axis=1)) #combine y
		frame = np.concatenate(frame, axis=0) #combine z
		
		return frame
		
	def hasMinDensity(self, tile):
		return self.getTileDensity(tile) >= (self.densityMinimum * tile.shape[0] * tile.shape[1] * tile.shape[2])
		
	def getTileDensity(self, tile):
		if self.data_flags[DATA_KEY_LOW]['channels'] > 1:
			tile = np.split(tile, [1], axis=-1)[0]
		return tile.sum( dtype=np.float64 )
	
	def getFrameTiles(self, index):
		''' returns the frame as tiles'''
		low, high = self.getDatum(index)
		return self.createTiles(low, self.tile_shape_low), self.createTiles(high, self.tile_shape_high)

#####################################################################################
# CHANNEL PARSING
#####################################################################################

	def parseChannels(self, channelString):
		''' arbitrary channel structure from string, expand if necessary. USE GLOBAL KEYS ^
			'd': default/ density; data that needs no special operations during augmentation
			'v[label](x|y|z)': vector/velocity; is transformed according to the augmentation
		'''
		#need this for low and high, +high only labels
		c = channelString.lower().split(',')
		for i in range(len(c)):
			c[i] = c[i].strip()
		
		c_types = {
			C_KEY_DEFAULT:[], #list of indices of default channels. e.g. for normal sim data [0]
			C_KEY_VELOCITY:[], #list of ordered (x,y,z; chnage to z,y,x to match data?) triples of velocity sets. e.g. for normal sim data [[1,2,3]]
			C_KEY_VORTICITY:[]
			}
		
		self.parse = {
			C_KEY_DEFAULT:self.parseCDefault,
			C_KEY_VELOCITY:self.parseCVelocity,
			C_KEY_VORTICITY:self.parseCVorticity
			}
			
		for i in range(len(c)):
			if len(c[i])==0: # check empty key
				self.TCError('empty channel key.'.format(i))
			try:
				self.parse[c[i][0]](c, i, c_types)
			except KeyError:
				self.TCError('channel {}: unknown channel key \"{}\".'.format(i, c[i]))
			
		# TODO check unused channels here 
		return c, c_types
	
	def parseCDefault(self, c, i, c_types):
		# check length
		if c[i]=='d':
			c_types[C_KEY_DEFAULT].append(i)
		else:
			self.TCError('channel {}: unknown channel key \"{}\".'.format(i, c[i]))
			
	def parseCVector(self, c, i, c_types, c_key, c_name='vector'):
		# c_key[label](x|y|z)
		if c[i][-1] == 'x' or c[i][-1] == 'y' or c[i][-1] == 'z':
			label = c[i][1:-1] #can be empty
			
			#get matching keys
			v_x = c_key+label+'x'
			v_y = c_key+label+'y'
			v_z = c_key+label+'z'
			
			#check for duplicates
			if c.count(v_x)>1:
				self.TCError('Duplicate {} ({}) x-channel with label \"{}\": {}. Vector keys must be unique.'.format(c_name, c_key, label, v_x))
			if c.count(v_y)>1:
				self.TCError('Duplicate {} ({}) y-channel with label \"{}\": {}. Vector keys must be unique.'.format(c_name, c_key, label, v_y))
			if c.count(v_z)>1:
				self.TCError('Duplicate {} ({}) z-channel with label \"{}\": {}. Vector keys must be unique.'.format(c_name, c_key, label, v_z))
			
			#check missing
			if c.count(v_x)==0:
				self.TCError('Missing {} ({}) x-channel with label \"{}\": {}'.format(c_name, c_key, label, v_x))
			if c.count(v_y)==0:
				self.TCError('Missing {} ({}) y-channel with label \"{}\": {}'.format(c_name, c_key, label, v_y))
			if self.dim==3 and c.count(v_z)==0:
				self.TCError('Missing {} ({}) z-channel with label \"{}\": {}'.format(c_name, c_key, label, v_z))
			
			if c[i][-1] == 'x':
				if(c.count(v_z)==0 and self.dim==2):
					c_types[C_KEY_VELOCITY].append([c.index(v_x),c.index(v_y)])
				else:
					c_types[C_KEY_VELOCITY].append([c.index(v_x),c.index(v_y),c.index(v_z)])
		
		# check wrong suffix
		else:
			self.TCError('Channel {}, \"{}\": unknown {} ({}) channel suffix \"{}\". Valid suffixes are \"x\", \"y\", \"z\".'.format(i, c[i], c_name, c_key, c[i][-1]))
			
		
	
	def parseCVelocity(self, c, i, c_types):
		# C_KEY_VELOCITY[label](x|y|z)
		self.parseCVector(c, i, c_types, C_KEY_VELOCITY, 'velociy')
	
	def parseCVorticity(self, c, i, c_types):
		# C_KEY_VELOCITY[label](x|y|z)
		self.parseCVector(c, i, c_types, C_KEY_VELOCITY, 'vorticity')

#####################################################################################
# ERROR HANDLING
#####################################################################################
	
	def TCError(self, msg):
		raise TilecreatorError(msg)
		
class TilecreatorError(Exception):
	''' Tilecreator errors '''
	

#####################################################################################
# IMAGE OUTPUT
#####################################################################################


# save summary images of all channels in a batch generated by the tile creator
# projects 3D data onto different axes, data has to be B-ZYX-C format
batchCounterGlob=0
def savePngsBatch(low,high, TC, path, batchCounter=-1, save_vels=False, dscale=1., vscale=1.):
	global batchCounterGlob
	if(low.shape[1]==1):
		dim=2
	else:
		dim=3

	# figure out good tile size, and project all axes for 3D
	batch = low.shape[0]
	tileX = 4 if batch>=4 else batch
	if batch%tileX != 0: tileX = batch
	# file names
	if batchCounter < 0:
		batchCounter = batchCounterGlob
		batchCounterGlob += 1
	path = path+"batch{:04d}_".format(batchCounter)

	# show scalar channels, in tiled images 
	aNames = ["xy_","xz_","yz_"]
	for axis in range(1 if dim==2 else 3):
		suff = aNames[axis]

		if dim == 3:
			highD = np.average(high, axis=axis+1) * dscale
			lowD  = np.average(low, axis=axis+1) * dscale
		if dim == 2:
			highD = high*brightness ; lowD = low * dscale
			lowD.shape = (batch, tll, tll, cl)
			highD.shape = (batch, tlh, tlh, ch)

		# note - outputs all channels as images, also vel channels...
		clout = np.arange(low.shape[4])
		savePngsGrayscale(lowD, path+'low_'+suff, tiles_in_image=[batch//tileX,tileX], channels=clout )
		chout = np.arange(high.shape[4])
		savePngsGrayscale(tiles=highD, path=path+'high_'+suff, imageCounter=0, tiles_in_image=[batch//tileX,tileX], channels=chout )

	# plot velocities , for individual samples
	if save_vels:
		for i in range(low.shape[0]):
			saveVelChannels(low[i], TC.c_lists[DATA_KEY_LOW][C_KEY_VELOCITY], path=path+'low_vel_i{:02d}_'.format(i), name="", scale=vscale )
		for i in range(high.shape[0]):
			saveVelChannels(high[i], TC.c_lists[DATA_KEY_HIGH][C_KEY_VELOCITY], path=path+'high_vel_i{:02d}_'.format(i), name="", scale=vscale )


# simpler function to output multiple tiles into grayscale pngs
def savePngsGrayscale(tiles, path, imageCounter=0, tiles_in_image=[1,1], channels=[0], save_gif=False, plot_vel_x_y=False, save_rgb=None, rgb_interval=[-1,1]):
	'''
		tiles_in_image: (y,x)
		tiles: shape: (tile,y,x,c)
	'''
	tilesInImage = tiles_in_image[0]*tiles_in_image[1]
	if len(tiles)%tilesInImage!=0: 
		print('ERROR: number of tiles does not match tiles per image')
		return
	tiles = np.asarray(tiles)
	noImages = len(tiles)//tilesInImage
	if save_gif:
		gif=[]
		
	for image in range(noImages):
		img = []
		#combine tiles to image
		for y in range(tiles_in_image[0]):
			offset=image*tilesInImage + y*tiles_in_image[1]
			img.append(np.concatenate(tiles[offset:offset+tiles_in_image[1]],axis=1)) #combine x
		img = np.concatenate(img, axis=0) #combine y
		# move channels to first dim.
		img_c = np.rollaxis(img, -1, 0)
		if len(img_c)>1 and (plot_vel_x_y or save_rgb!=None):
			if plot_vel_x_y: saveVel(img, path, imageCounter+image)
			if save_rgb!=None: saveRGBChannels(img,path, save_rgb,value_interval=rgb_interval, imageCounter=imageCounter+image)
		if len(channels) == 1:
			scipy.misc.toimage(img_c[channels[0]], cmin=0.0, cmax=1.0).save(path + 'img_{:04d}.png'.format(imageCounter*noImages+image))
		else:
			for i in channels:
				scipy.misc.toimage(img_c[i], cmin=0.0, cmax=1.0).save(path + 'img_{:04d}_c{:04d}.png'.format(imageCounter*noImages+image, i))

	
# store velocity as quiver plot
def saveVel(tile, path, imageCounter=0, name='vel-x-y'):
	# origin is in upper left corner, transform acordingly
	y, x = np.mgrid[-tile.shape[0]:0, 0:tile.shape[1]]
	vx = None; vy = None
	if tile.shape[-1]==4:
		d, vx, vy, vz = np.split(tile, 4, -1)
	elif tile.shape[-1]==2:
		vx, vy  = np.split(tile, 2, -1)
	else:
		print('ERROR: unknown nr of channels for vel input '+format(tile.shape))
	vx = vx[::-1, ...]
	vy = vy[::-1, ...]
	if found_matplotlib:
		matplotlib.pyplot.quiver(x,y,vx.flatten(),vy.flatten(), units = 'xy', scale = 1)
		matplotlib.pyplot.axis('equal')
		matplotlib.pyplot.savefig(path + '{}_{:04d}.png'.format(name,imageCounter))
		matplotlib.pyplot.clf()
	
# save velocity channels from the tilecreator with multiple axis projections (uses saveVel)
def saveVelChannels(data, c_idx, path, average=False, scale=1.0, normalize=True, name=''):
	channels = np.split(data, data.shape[-1], -1)

	vpath = path
	vcnt = 0
	for v in c_idx:
		if(len(c_idx)>1):
			vpath = path + "vc{}".format(vcnt)
		vcnt += 1
		# compute scale factor
		vscale = scale
		if normalize:
			vavg = np.concatenate( [ channels[v[0]],channels[v[1]] ], -1)
			if(len(v))>2: vavg = np.concatenate( [ vavg, channels[v[2]] ],-1)
			vscale *= (1./(np.max( vavg )+1e-10)) # normalize

		vavg = np.concatenate( [ channels[v[0]],channels[v[1]] ] , -1)
		vavg = np.average(vavg, axis=0)
		vavg *= vscale
		saveVel(vavg, path=vpath, name='_xy' )

		if(len(v))>2: # also output xz,yz
			vavg = np.concatenate( [ channels[v[0]],channels[v[2]] ] , -1)
			vavg = np.average(vavg, axis=1)
			vavg *= vscale
			saveVel(vavg, path=vpath, name='_xz' )

			vavg = np.concatenate( [ channels[v[1]],channels[v[2]] ] , -1)
			vavg = np.average(vavg, axis=2)
			vavg *= vscale
			saveVel(vavg, path=vpath, name='_yz' )


def saveRGBChannels(data, path, channel_list, imageCounter=0, value_interval=[-1,1]):
	"""
		data: shape[y,x,c]
		channels: list of triples of channel ids saved as RGB image
	"""

	cmin = value_interval[0]
	cmax = value_interval[1]
	num_channels = data.shape[-1]
	
	channels = np.split(data, num_channels, -1)
	for i in channel_list:
		if len(i)==2:
			img = np.concatenate([channels[i[0]], channels[i[1]], np.ones_like(channels[i[0]])*cmin], -1)
		else:
			img = np.concatenate([channels[i[0]], channels[i[1]], channels[i[2]]], -1)
		scipy.misc.toimage(img, cmin=-1.0, cmax=1.0).save(path + 'img_rgb_{:04d}.png'.format(imageCounter))

def save3DasUni(tiles, path, motherUniPath, imageCounter=0, tiles_in_image=[1,1]):
	'''
		tiles_in_image: (y,x)
		tiles: shape: (image,y,x,c)
	'''

	tilesInImage = tiles_in_image[0]*tiles_in_image[1]
	if len(tiles)%tilesInImage!=0: 
		print('ERROR: number of tiles does not match tiles per image')
		return
	tiles = np.asarray(tiles)
	noImages = len(tiles)//tilesInImage

	tiles = np.reshape(tiles,(len(tiles), tiles.shape[1], tiles.shape[2], tiles.shape[-1]))

	#only save the average
	img_all = []
	for image in range(noImages):
		img = []
		#combine tiles to image
		for y in range(tiles_in_image[0]):
			offset=( image) * tilesInImage + (y)*tiles_in_image[1]
			img.append(np.concatenate(tiles[offset:offset+tiles_in_image[1]],axis=2)) #combine y
		img = np.array(img)
		img = np.concatenate(img, axis=1) #combine x
		img = np.array(img)
		# move channels to first dim.
		img_c = np.rollaxis(img, 0, 0)
		img_all.append(img_c)
	img_all = np.array(img_all)
	img_all = np.concatenate(img_all, axis=0)
	img_all = np.array(img_all)
	TDarrayToUni(img_all, path + 'source_{:04d}.uni'.format(imageCounter), motherUniPath, img_all.shape[0], img_all.shape[1], img_all.shape[2])


def TDarrayToUni(input, savePath, motherUniPath, imageHeight, imageWidth, imageDepth, is_vel=False):
	head, _ = uniio.readUni(motherUniPath)
	head['dimX'] = imageWidth
	head['dimY'] = imageHeight
	head['dimZ'] = imageDepth

	if not is_vel:
		fixedArray = np.zeros((imageHeight, imageWidth, imageDepth), dtype='f')
		for x in range(0, imageHeight):
			for y in range(0, imageWidth):
				for z in range(0, imageDepth):
					fixedArray[x][y][z] = input[imageDepth - 1 - z][y][(imageHeight - 1) - x]
	else:
		fixedArray = np.zeros((imageHeight, imageWidth, imageDepth, 3), dtype='f')
		for x in range(0, imageHeight):
			for y in range(0, imageWidth):
				for z in range(0, imageDepth):
					fixedArray[x][y][z] = input[imageDepth - 1 - z][y][(imageHeight - 1) - x]

	uniio.writeUni(savePath, head, fixedArray)


# ******************************************************************************
# faster functions, batch operations
#

# grid interpolation method, order: only linear tested
# for velocity, macgridbatch.shape should be [b,z,y,x,3] (in 2D, should be [b,1,ny,nx,3])
# for density , macgridsource.shape should be [z,y,x,1] (in 2D, should be [1,ny,nx,1])
def gridInterpolBatch(macgridbatch, targetshape, order=1):
	assert (targetshape[-1] == macgridbatch.shape[-1])  # no interpolation between channels
	assert (len(targetshape) == 5 and len(macgridbatch.shape) == 5)
	dim = 3
	if (macgridbatch.shape[1] == 1 and targetshape[1] == 1):  dim = 2
	
	x_ = np.linspace(0.5, targetshape[3] - 0.5, targetshape[3])
	y_ = np.linspace(0.5, targetshape[2] - 0.5, targetshape[2])
	z_ = np.linspace(0.5, targetshape[1] - 0.5, targetshape[1])
	c_ = np.linspace(0, targetshape[4] - 1, targetshape[4])  # no interpolation between channels
	b_ = np.linspace(0, targetshape[0] - 1, targetshape[0])  # no interpolation between batches
	
	b, z, y, x, c = np.meshgrid(b_, z_, y_, x_, c_, indexing='ij')
	
	# scale
	fx = float(macgridbatch.shape[3]) / targetshape[3]
	fy = float(macgridbatch.shape[2]) / targetshape[2]
	fz = float(macgridbatch.shape[1]) / targetshape[1]
	
	mactargetbatch = scipy.ndimage.map_coordinates(macgridbatch, [b, z * fz, y * fy, x * fx, c], order=order, mode='reflect')
	
	return mactargetbatch;


# macgrid_batch shape b, z, y, x, 3
# return a matrix in size [b,z,y,x,3] ( 2D: [b,y,x,2]), last channel in z-y-x order!(y-x order for 2D)
def getMACGridCenteredBatch(macgrid_batch, is3D):
	_bn, _zn, _yn, _xn, _cn = macgrid_batch.shape
	
	valid_idx = list(range(1, _xn))
	valid_idx.append(_xn - 1)
	add_x = macgrid_batch.take(valid_idx, axis=3)[:, :, :, :, 0]  # shape, [b,z,y,x]
	valid_idx = list(range(1, _yn))
	valid_idx.append(_yn - 1)
	add_y = macgrid_batch.take(valid_idx, axis=2)[:, :, :, :, 1]  # shape, [b,z,y,x]
	
	add_y = add_y.reshape([_bn, _zn, _yn, _xn, 1])
	add_x = add_x.reshape([_bn, _zn, _yn, _xn, 1])
	if (is3D):
		valid_idx = list(range(1, _zn))
		valid_idx.append(_zn - 1)
		add_z = macgrid_batch.take(valid_idx, axis=1)[:, :, :, :, 2]  # shape, [b,z,y,x]
		add_z = add_z.reshape([_bn, _zn, _yn, _xn, 1])
		resultgrid = 0.5 * (macgrid_batch[:, :, :, :, ::-1] + np.concatenate((add_z, add_y, add_x), axis=-1))
		return resultgrid.reshape([_bn, _zn, _yn, _xn, 3])
	
	resultgrid = 0.5 * (macgrid_batch[:, :, :, :, -2::-1] + np.concatenate((add_y, add_x), axis=4))
	return resultgrid.reshape([_bn, _yn, _xn, 2])


# macgrid_batch shape b, z, y, x, 3 ( b,1,y,x,3 for 2D )
# return the re-sampling positions as a matrix, in size of [b,z,y,x,3] ( 2D: [b,y,x,2])
def getSemiLagrPosBatch(macgrid_batch, dt, cube_len_output=-1):  # check interpolation later
	assert (len(macgrid_batch.shape) == 5)
	_bn, _zn, _yn, _xn, _cn = macgrid_batch.shape
	assert (_cn == 3)
	is3D = (_zn > 1)
	if (cube_len_output == -1): cube_len_output = _xn
	factor = float(_xn) / cube_len_output
	x_ = np.linspace(0.5, int(_xn / factor + 0.5) - 0.5, int(_xn / factor + 0.5))
	y_ = np.linspace(0.5, int(_yn / factor + 0.5) - 0.5, int(_yn / factor + 0.5))
	interp_shape = [_bn, int(_zn / factor + 0.5), int(_yn / factor + 0.5), int(_xn / factor + 0.5), 3]
	if (not is3D): interp_shape[1] = 1
	
	if (is3D):
		z_ = np.linspace(0.5, int(_zn / factor + 0.5) - 0.5, int(_zn / factor + 0.5))
		z, y, x = np.meshgrid(z_, y_, x_, indexing='ij')
		posArray = np.stack((z, y, x), axis=-1)  # shape, z,y,x,3
		tarshape = [1, int(_zn / factor + 0.5), int(_yn / factor + 0.5), int(_xn / factor + 0.5), 3]
	else:
		y, x = np.meshgrid(y_, x_, indexing='ij')
		posArray = np.stack((y, x), axis=-1)  # shape, y,x,2
		tarshape = [1, int(_yn / factor + 0.5), int(_xn / factor + 0.5), 2]
	posArray = posArray.reshape(tarshape)
	if (cube_len_output == _xn):
		return (posArray - getMACGridCenteredBatch(macgrid_batch, is3D) * dt)
	# interpolate first
	
	inter_mac_batch = gridInterpolBatch(macgrid_batch, interp_shape, 1)
	inter_mac_batch = getMACGridCenteredBatch(inter_mac_batch, is3D) / factor
	return (posArray - (inter_mac_batch) * dt)


# error! muss not shuffle for fluid data loader!
def selectRandomTempoTiles(self, selectionSize, isTraining=True, augment=False, n_t=3, dt=0.5, adv_flag = 1.0):
	'''
		main method to create coherent baches
		Return:
			shape: [n_t * selectionSize//n_t, z, y, x, channels]
			if 2D z = 1
			channels: density, [vel x, vel y, vel z], [pos x, pox y, pos z]
	'''
	batch_sz = int( max( 1, selectionSize // n_t) )
	batch_low, batch_high = self.selectRandomTiles(batch_sz, isTraining, augment, tile_t = n_t)
	real_batch_sz = batch_sz * n_t
	ori_input_shape = batch_low.reshape((batch_sz, self.tileSizeLow[0], self.tileSizeLow[1], self.tileSizeLow[2], n_t, -1))
	ori_input_shape = np.transpose(ori_input_shape, (0,4,1,2,3,5))
	ori_input_shape = ori_input_shape.reshape((real_batch_sz, self.tileSizeLow[0], self.tileSizeLow[1], self.tileSizeLow[2], -1))
	vel_pos_high_inter = None
	if adv_flag:	
		# TODO check velocity channels and 3D 
		macgrid_input = ori_input_shape[:, :, :, :, self.c_lists[DATA_KEY_LOW][C_KEY_VELOCITY][0]]
		macgrid_input = macgrid_input.reshape( (real_batch_sz, self.tileSizeLow[0], self.tileSizeLow[1], self.tileSizeLow[2], 3))
		dtArray = np.array([i * dt for i in range(n_t // 2, -n_t // 2, -1)] * batch_sz, dtype=np.float32)
		if (self.dim == 2):
			dtArray = dtArray.reshape((-1, 1, 1, 1))
		else:
			dtArray = dtArray.reshape((-1, 1, 1, 1, 1))  

		vel_pos_high_inter = getSemiLagrPosBatch(macgrid_input, dtArray, self.tileSizeHigh[1]).reshape((real_batch_sz, -1))
	
	# return reshape_input, selectedOutputs, vel_pos_high_inter
	batch_high = batch_high.reshape((batch_sz, self.tileSizeHigh[0], self.tileSizeHigh[1], self.tileSizeHigh[2], n_t, -1))
	batch_high = np.transpose(batch_high, (0,4,1,2,3,5))
	return ori_input_shape.reshape((real_batch_sz, -1)), batch_high.reshape((real_batch_sz, -1)), vel_pos_high_inter

TileCreator.selectRandomTempoTiles = selectRandomTempoTiles

def pngs_to_gif(path, start_idx=0, end_idx=199, step=1, fps=20, mask="img_%04d.png"):
	print("creating gif from {} to {} with {} fps".format(mask % start_idx, mask % end_idx, fps))
	with imageio.get_writer(path + 'step%d.gif'%step, mode='I', fps=fps) as writer:
		for i in range(start_idx - 1, end_idx, step):
			image = imageio.imread(path+(mask % (i+1)))
			writer.append_data(image)

