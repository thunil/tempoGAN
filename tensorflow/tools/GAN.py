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

import numpy as np
import tensorflow as tf
import math
from keras import backend as kb

class GAN(object):
	#---------------------------------------------------------------------------------
	def __init__(self, _image, bn_decay=0.999):
		self.layer = _image 
		self.batch_size = tf.shape(_image)[0]
		self.DOFs = 0
		# stack 
		self.preFlatShapes = []
		self.weight_stack = []
		self.layer_num = 0
		self.layer_num_gen = 0
		self.layer_num_disc = 0
		
		self.bn_decay=bn_decay
		
		print("Input: {}".format(self.layer.get_shape()))
	
	#---------------------------------------------------------------------------------
	# thanks to http://robromijnders.github.io/tensorflow_basic/
	def weight_image(self):
		W = self.weight_stack[-1]
		# compute size of the image
		s = W.get_shape()
		out_channels = 1
		if int(s[3]) % 3 == 0:
			out_channels = 3
		print("Shape {}".format(s))
		weight_patches = int(s[2]) * int(s[3]) / out_channels # e.g. the number of [3,3] patches in a CNN
		side_length = int(math.ceil(math.sqrt(weight_patches))) # image side length (in patches)
		image_patches = side_length * side_length # max number of patches that fit in the image
		# split into per filter weights
		ws = []
		ws_dim3 = tf.split(3, s[3] / out_channels, W) # e.g. [ [3,3,3,1], [3,3,3,1], ... ]
		for w in ws_dim3:
			# split these further
			ws.extend(tf.split(2, s[2], w))  # e.g. [ [3,3,1,1], [3,3,1,1], ... ]
		# pad image
		padding = image_patches - weight_patches
		for i in range(padding):
			ws.append(tf.zeros([s[0], s[1], 1, out_channels]))
		# build rows of image
		rows = []
		for i in range(side_length):
			start = i * side_length
			end = start + side_length
			rows.append(tf.concat(axis=0, values=ws[start:end]))
		# combine rows to image
		image = tf.concat(axis=1, values=rows) # [sidelength * ]
		s = [int(image.get_shape()[0]), int(image.get_shape()[1])]
		image = tf.reshape(image, [1, s[0], s[1], out_channels])
		image = tf.image.resize_images(image, [int(s[1] * 50), int(s[0] * 50)], 1)
		image_tag = "l" + str(self.layer_num) + "_weight_image"
		tf.image_summary(image_tag, image)
		print("Image Summary: save weights as image")
		
	#---------------------------------------------------------------------------------
	# outChannels: int
	# _patchShape: 2D: [H,W]; 3D: [D,H,W]
	# stride 3D: if 1D: [DHW], if 2D:[D,HW], if 3D:[D,H,W]
	# returns both normalized and linearized versions
	def convolutional_layer(self, outChannels, _patchShape, activation_function=tf.nn.tanh, stride=[1], name="conv",reuse=False, batch_norm=False, train=None, in_layer=None):
		if in_layer==None:
			in_layer = self.layer
		with tf.variable_scope(name):
			self.layer_num += 1
			# set the input and output dimension
			inChannels = int(in_layer.get_shape()[-1])
			#outChannels = int(inChannels * _filterSpread)
			# create a weight matrix
			if len(_patchShape) == 2:
				W = self.weight_variable([_patchShape[0], _patchShape[1], inChannels, outChannels], name=name)
				self.layer = self.conv2d(in_layer, W, stride)
				self.DOFs += _patchShape[0]* _patchShape[1]* inChannels* outChannels
			elif len(_patchShape) == 3:
				W = self.weight_variable([_patchShape[0], _patchShape[1], _patchShape[2], inChannels, outChannels], name=name)
				self.layer = self.conv3d(in_layer, W, stride)
				self.DOFs += _patchShape[0]* _patchShape[1]* _patchShape[2]* inChannels* outChannels
				#batch_norm = False
				
			self.weight_stack.append(W)
			# create a bias vector
			b = self.bias_variable([outChannels], name=name)
			self.layer = self.layer + b
			self.DOFs += outChannels
			
			if batch_norm:
				#self.layer = self.conv_batch_norm(self.layer, train=train)
				self.layer = tf.contrib.layers.batch_norm(self.layer, decay=self.bn_decay, scale=True, scope=tf.get_variable_scope(), reuse=reuse, fused=False, is_training=train)
			layer_lin = self.layer
			if activation_function:
				self.layer = activation_function(self.layer)
			# user output
			if activation_function:
				print("Convolutional Layer \'{}\' {} ({}) : {}, BN:{}".format(name, W.get_shape(), activation_function.__name__,self.layer.get_shape(),batch_norm))
			else:
				print("Convolutional Layer \'{}\' {} ({}) : {}, BN:{}".format(name, W.get_shape(), 'None',self.layer.get_shape(),batch_norm))
			return self.layer, layer_lin
	
	#---------------------------------------------------------------------------------
	# s1: outChannels of intermediate conv layer
	# s2: outChannels of final and skip conv layer
	# filter: 2D: [H,W]; 3D: [D,H,W]
	# returns both normalized and linearized versions
	def residual_block(self, s1,s2, filter, activation_function=tf.nn.tanh, name="RB", reuse=False, batch_norm=False, train=None, in_layer=None):
		# note - leaky relu (lrelu) not too useful here
		if in_layer==None:
			in_layer = self.layer
		# convolutions of resnet block
		if len(filter) == 2:
			filter1 = [1,1]
		elif len(filter) == 3:
			filter1 = [1,1,1]
			
		print("Residual Block:")
		A,_ = self.convolutional_layer(s1, filter, activation_function, stride=[1], name=name+"_A", in_layer=in_layer, reuse=reuse, batch_norm=batch_norm, train=train)
		B,_ = self.convolutional_layer(s2, filter, None               , stride=[1], name=name+"_B",                    reuse=reuse, batch_norm=batch_norm, train=train)
		# shortcut connection
		s,_ = self.convolutional_layer(s2, filter1, None              , stride=[1], name=name+"_s", in_layer=in_layer, reuse=reuse, batch_norm=batch_norm, train=train)
		
		self.layer = tf.add( B, s)
		layer_lin = self.layer
		if activation_function:
			self.layer = activation_function(self.layer )
			
		return self.layer, layer_lin
	
	
	#---------------------------------------------------------------------------------
	# 2 x 2 max pool operation
	def max_pool(self, window_size=[2], window_stride=[2]):
		if len(self.layer.get_shape()) == 4:
			self.layer = tf.nn.max_pool(self.layer, ksize=[1, window_size[0], window_size[0], 1], strides=[1, window_stride[0], window_stride[0], 1], padding="VALID")
		elif len(self.layer.get_shape()) == 5:
			self.layer = tf.nn.max_pool3d(self.layer, ksize=[1, window_size[0], window_size[0], window_size[0], 1], strides=[1, window_stride[0], window_stride[0], window_stride[0], 1], padding="VALID")
		# user output
		print("Max Pool {}: {}".format(window_size, self.layer.get_shape()))
		return self.layer
	
	#---------------------------------------------------------------------------------
	def avg_pool(self, window_size=[2], window_stride=[2]):
		if len(self.layer.get_shape()) == 4:
			self.layer = tf.nn.avg_pool(self.layer, ksize=[1, window_size[0], window_size[0], 1], strides=[1, window_stride[0], window_stride[0], 1], padding="VALID")
		elif len(self.layer.get_shape()) == 5:
			self.layer = tf.nn.avg_pool3d(self.layer, ksize=[1, window_size[0], window_size[0], window_size[0], 1], strides=[1, window_stride[0], window_stride[0], window_stride[0], 1], padding="VALID")
		# user output
		print("Avg Pool {}: {}".format(window_size, self.layer.get_shape()))
		return self.layer
	
	#---------------------------------------------------------------------------------
	# make layer flat
	# e.G. [1, 4, 4, 2] -> [1, 32]
	def flatten(self):
		# get unflat shape
		layerShape = self.layer.get_shape()
		self.preFlatShapes.append(layerShape)
		# compute flat size
		flatSize = int(layerShape[1]) * int(layerShape[2]) * int(layerShape[3])
		if len(layerShape) == 5:
			flatSize *= int(layerShape[4])
		# make flat
		self.layer = tf.reshape(self.layer, [-1, flatSize])
		# user output
		print("Flatten: {}".format(self.layer.get_shape()))
		return flatSize
	
	#---------------------------------------------------------------------------------
	def fully_connected_layer(self, _numHidden, _act, name="full"):
		with tf.variable_scope(name):
			self.layer_num += 1
			# get previous layer size
			numInput = int(self.layer.get_shape()[1])
			# build layer variables
			W = self.weight_variable([numInput, _numHidden], name=name)
			b = self.bias_variable([_numHidden], name=name)
			self.DOFs += numInput*_numHidden + _numHidden
			# activate
			self.layer = tf.matmul(self.layer, W) + b
			if _act:
				self.layer = _act(self.layer)  # ??
			# user output
			if _act:
				print("Fully Connected Layer \'{}\': {}".format(name, self.layer.get_shape()))
			else:
				print("Linear Layer \'{}\': {}".format(name, self.layer.get_shape()))
			return self.layer
	
	#---------------------------------------------------------------------------------
	# make layer 3D (from previously stored)
	# e.G. [1, 32] -> [1, 4, 4, 2]
	def unflatten(self):
		unflatShape = self.preFlatShapes.pop()
		if len(unflatShape) == 4:
			unflatShape = [-1, int(unflatShape[1]), int(unflatShape[2]), int(unflatShape[3])]
		elif len(unflatShape) == 5:
			unflatShape = [-1, int(unflatShape[1]), int(unflatShape[2]), int(unflatShape[3]), int(unflatShape[4])]
		self.layer = tf.reshape(self.layer, unflatShape)
		print("Unflatten: {}".format(self.layer.get_shape()))
		return self.layer


	#---------------------------------------------------------------------------------
	# inverse of 2 x 2 max pool , note window size&stride given as [x,y] pair
	# does not support 3D
	def max_depool(self, in_layer=None, depth_factor =2, height_factor=2, width_factor=2):
		if in_layer==None:
			in_layer = self.layer

		if len(self.layer.get_shape()) == 4:
			self.layer = kb.resize_images(self.layer, height_factor, width_factor, 'channels_last')
			print("Max Depool : {}".format(self.layer.get_shape()))	
		if len(self.layer.get_shape()) == 5:
			self.layer = kb.resize_volumes(self.layer, depth_factor, height_factor, width_factor, 'channels_last')
			print("Max Depool : {}".format(self.layer.get_shape()))
		return self.layer

	#---------------------------------------------------------------------------------
	# resizes H and W dimensions of NHWC or NDHWC (scale only H and W for 3D DHW data)
	# https://stackoverflow.com/questions/43814367/resize-3d-data-in-tensorflow-like-tf-image-resize-images
	def avg_depool(self, window_size=[2, 2], window_stride=[2,2]):
		is3D = False
		if len(self.layer.get_shape()) == 5: # 3D data, merge D into C to have a 4D tensor (like 2D data)
			is3D = True
			self.layer = tf.transpose(self.layer, [0,2,3,1,4]) # NDHWC -> NHWDC
			s=self.layer.get_shape() # NHWDC
			self.layer = tf.reshape(self.layer, [-1, int(s[1]), int(s[2]), int(s[3])*int(s[4])]) # NHWDC -> NHW(D*C)
			
		outWidth = self.layer.get_shape()[2] * window_stride[0] + window_size[0] - window_stride[0]
		outHeight = self.layer.get_shape()[1] * window_stride[1] + window_size[1] -  window_stride[1]
		self.layer = tf.image.resize_images(self.layer, [int(outHeight), int(outWidth)], 0) #0 = ResizeMethod.BILINEAR
		
		if is3D: # recover D dimension
			self.layer = tf.reshape(self.layer, [-1, int(outHeight), int(outWidth), int(s[3]), int(s[4])])
			self.layer = tf.transpose(self.layer, [0,3,1,2,4]) # -> NDHWC
			s=self.layer.get_shape() # NDHWC
			self.layer = tf.reshape(self.layer, [-1, int(s[1]), int(s[2]), int(s[3])*int(s[4])])# NDHWC ->	 NDH(W*C)
			self.layer = tf.image.resize_images(self.layer, [int(s[1]*window_size[0]), int(s[2])], 0) #0 = ResizeMethod.BILINEAR		
			self.layer = tf.reshape(self.layer, [-1, int(s[1]*window_size[0]), int(s[2]), int(s[3]), int(s[4])])# NDHWC

		print("Avg Depool {}: {}".format(window_size, self.layer.get_shape()))
		return self.layer

	#---------------------------------------------------------------------------------
	# outChannels: int
	# _patchShape: 2D: [H,W]; 3D: [D,H,W]
	# stride 3D: if 1D: [DHW], if 2D:[D,HW], if 3D:[D,H,W]
	def deconvolutional_layer(self, outChannels, _patchShape, activation_function=tf.nn.tanh, stride=[1], name="deconv",reuse=False, batch_norm=False, train=None, init_mean=0., strideOverride=None):
		if init_mean==1.:
			name = name+"_EXCLUDE_ME_"
		with tf.variable_scope(name):
			self.layer_num += 1
			shape = self.layer.get_shape()
			# spread channels
			inChannels = int(self.layer.get_shape()[-1])
			#outChannels = int(inChannels / _filterSpread) # must always come out even

			dcStride = stride
			if strideOverride is not None:
				dcStride = strideOverride

			if len(_patchShape) == 2:
				if len(stride) == 1:
					stride = [stride[0],stride[0]]
				# create a weight matrix
				W = self.weight_variable([_patchShape[0], _patchShape[1], outChannels, inChannels], name=name, init_mean=init_mean)
				self.layer = self.deconv2d(self.layer, W, [self.batch_size, int(shape[1]*stride[0]), int(shape[2]*stride[1]), outChannels], dcStride)
				self.DOFs += _patchShape[0]* _patchShape[1]* outChannels* inChannels
			if len(_patchShape) == 3:
				if len(stride) == 1:
					stride = [stride[0],stride[0],stride[0]]
				elif len(stride) == 2:
					stride = [stride[0],stride[1],stride[1]]
				# create a weight matrix
				W = self.weight_variable([_patchShape[0], _patchShape[1], _patchShape[2], outChannels, inChannels], name=name, init_mean=init_mean)
				self.layer = self.deconv3d(self.layer, W, [self.batch_size, int(shape[1]*stride[0]), int(shape[2]*stride[1]), int(shape[3]*stride[2]), outChannels], dcStride)
				self.DOFs += _patchShape[0]* _patchShape[1]* _patchShape[2]* outChannels* inChannels
				#batch_norm = False
				
			# create a bias vector
			b = self.bias_variable([outChannels], name=name)
			self.layer = self.layer + b
			self.DOFs += outChannels
			
			if len(_patchShape) == 2:
				self.layer = tf.reshape(self.layer, [-1, int(shape[1]*stride[0]), int(shape[2]*stride[1]), outChannels])
			if len(_patchShape) == 3:
				self.layer = tf.reshape(self.layer, [-1, int(shape[1]*stride[0]), int(shape[2]*stride[1]), int(shape[3]*stride[2]), outChannels])
				
			if batch_norm:
				self.layer = tf.contrib.layers.batch_norm(self.layer, decay=self.bn_decay, scale=True, scope=tf.get_variable_scope(), reuse=reuse, fused=False, is_training=train)
			layer_lin = self.layer
			if activation_function:
				self.layer = activation_function(self.layer)
			# user output
			if activation_function:
				print("Deconvolutional Layer \'{}\' {} ({}): {}, BN:{}".format(name, W.get_shape(), activation_function.__name__, self.layer.get_shape(),batch_norm))
			else:
				print("Deconvolutional Layer \'{}\' {} ({}): {}, BN:{}".format(name, W.get_shape(), 'None', self.layer.get_shape(),batch_norm))
			return self.layer, layer_lin
	
	#---------------------------------------------------------------------------------
	#adds noise to the current layer
	#channels: number of noise channels to add, uses channels of current layer if < 1
	def noise(self, channels=-1):
		shape=tf.shape(self.layer)
		if channels > 0:
			shape[-1] = channels
		noise = tf.random_normal(shape=shape, mean=0.0, stddev=0.04, dtype=tf.float32)
		self.layer = tf.concat([self.layer, noise], axis=-1)
		print("Noise {}: {}".format(noise.get_shape(), self.layer.get_shape()))
		return self.layer

	#---------------------------------------------------------------------------------
	#adds the given tensor to self.layer on axis -1(channels)
	def concat(self, layer):
		self.layer = tf.concat(values=[self.layer, layer], axis=-1)
		print("Concat {}: {}".format(layer.get_shape(), self.layer.get_shape()))
		return self.layer
	#---------------------------------------------------------------------------------
	#applys the given operation to self.layer
	def apply(self, op):
		self.layer = op(self.layer)
		print("Apply \'{}\': {}".format(op.__name__, self.layer.get_shape()))
		return self.layer
	#---------------------------------------------------------------------------------
	def dropout(self, keep_prob):
		self.layer = tf.nn.dropout(self.layer, keep_prob)
		print("Dropout: {}".format(self.layer.get_shape()))
		return self.layer
	
	#---------------------------------------------------------------------------------
	def y(self):
		return self.layer

	#---------------------------------------------------------------------------------
	def getDOFs(self):
		return self.DOFs

	#---------------------------------------------------------------------------------
	# generate random valued weight field
	def weight_variable(self, shape, name="w", init_mean=0.):
		# use tf.get_variable() instead of tf.Variable() to be able to reuse variables
		s = 0.04
		if init_mean==1.:
			s = 0. # enforce value
		v = tf.get_variable("weight", shape, initializer=tf.random_normal_initializer(stddev=s, mean=init_mean))
		return v

	#---------------------------------------------------------------------------------
	# gemerate biases for the nodes
	def bias_variable(self, shape, name="b"):
		return tf.get_variable("bias", shape, initializer=tf.constant_initializer(0.1))

	#---------------------------------------------------------------------------------
	def conv2d(self, x, W, stride=[1]):
		if len(stride) == 1: #[HW]
			strides = [1, stride[0], stride[0], 1]
		elif len(stride) == 2: #[H,W]
			strides = [1, stride[0], stride[1], 1]
		return tf.nn.conv2d(x, W, strides=strides, padding="SAME")
		
	def conv3d(self, x, W, stride=[1]):
		if len(stride) == 1: #[DHW]
			strides = [1, stride[0], stride[0], stride[0], 1]
		elif len(stride) == 2: #[D,HW] for use when striding time and space separately
			strides = [1, stride[0], stride[1], stride[1], 1]
		elif len(stride) == 3: #[D,H,W]
			strides = [1, stride[0], stride[1], stride[2], 1]
		return tf.nn.conv3d(x, W, strides=strides, padding="SAME")

	#---------------------------------------------------------------------------------
	def deconv2d(self, x, W, output_shape, stride=[1]):
		if len(stride) == 1:
			strides = [1, stride[0], stride[0], 1]
		elif len(stride) == 2:
			strides = [1, stride[0], stride[1], 1]
		return tf.nn.conv2d_transpose(x, W, output_shape=output_shape, strides=strides, padding="SAME")
		
	def deconv3d(self, x, W, output_shape, stride=[1]):
		if len(stride) == 1:
			strides = [1, stride[0], stride[0], stride[0], 1]
		elif len(stride) == 2: # for use when striding time and space separately
			strides = [1, stride[0], stride[1], stride[1], 1]
		elif len(stride) == 3:
			strides = [1, stride[0], stride[1], stride[2], 1]
		return tf.nn.conv3d_transpose(x, W, output_shape=output_shape, strides=strides, padding="SAME")

	def variable_summaries(self, var, name):
		"""Attach a lot of summaries to a Tensor."""
		with tf.name_scope('summaries'):
			mean = tf.reduce_mean(var)
			tf.summary.scalar('mean/' + name, mean)
			with tf.name_scope('stddev'):
				stddev = tf.sqrt(tf.reduce_sum(tf.square(var - mean)))
			tf.summary.scalar('sttdev/' + name, stddev)
			tf.summary.scalar('max/' + name, tf.reduce_max(var))
			tf.summary.scalar('min/' + name, tf.reduce_min(var))
			tf.summary.histogram(name, var)
			
	
# thanks to https://github.com/bamos/dcgan-completion.tensorflow
def lrelu(x, leak=0.2, name="lrelu"):
    with tf.variable_scope(name):
        f1 = 0.5 * (1 + leak)
        f2 = 0.5 * (1 - leak)
        return f1 * x + f2 * abs(x)	

