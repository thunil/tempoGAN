# ----------------------------------------------------------------------------
#
# MantaFlow fluid solver framework
# Copyright 2017 Kiwon Um, Nils Thuerey
#
# This program is free software, distributed under the terms of the
# Apache License, Version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Utilities for creating the neural newtork of TensorFlow
#
# ----------------------------------------------------------------------------

import tensorflow as tf

tf.set_random_seed(1)

def parse_act(act):
	if act == 'tanh':    return tf.nn.tanh
	if act == 'sigmoid': return tf.nn.sigmoid
	if act == 'relu':    return tf.nn.relu
	return tf.nn.tanh

# init_weights = { 'w': {'mean': float, 'stddev': float}, 'b': {'value': float} }
def weight_variable(shape, init_weights=None):
	params = dict(shape=shape)
	if init_weights and 'w' in init_weights: params.update(init_weights['w'])
	initial = tf.truncated_normal(**params)
	return tf.Variable(initial)

def bias_variable(shape, init_weights=None):
	params = dict(value=0.1, shape=shape)
	if init_weights and 'b' in init_weights: params.update(init_weights['b'])
	initial = tf.constant(**params)
	return tf.Variable(initial)

def variable_summaries(var, name):
	"""Attach a lot of summaries to a Tensor."""
	with tf.name_scope('summaries'):
		mean = tf.reduce_mean(var)
		tf.summary.scalar(name='mean/' + name, tensor=mean)
		with tf.name_scope('stddev'):
			stddev = tf.sqrt(tf.reduce_sum(tf.square(var - mean)))

		tf.summary.scalar(name='stddev/' + name, tensor=stddev)
		tf.summary.scalar(name='max/' + name, tensor=tf.reduce_max(var))
		tf.summary.scalar(name='min/' + name, tensor=tf.reduce_min(var))
		tf.summary.histogram(name=name, values=var)

def nn_layer(input_tensor, input_dim, output_dim, layer_name,
			 init_weights=None, keep_prob=None, bn=False, is_training=True, act=tf.nn.tanh):
	with tf.name_scope(layer_name):
		with tf.name_scope('weights'):
			weights = weight_variable([input_dim, output_dim], init_weights)
			variable_summaries(weights, layer_name + '/weights')
		with tf.name_scope('biases'):
			biases = bias_variable([output_dim], init_weights)
			variable_summaries(biases, layer_name + '/biases')
		with tf.name_scope('Wx_plus_b' if not bn else 'BN_Wx_plus_b'):
			if bn:
				# with batch normalization
				epsilon =  1e-4
				z_bn        = tf.matmul(input_tensor, weights if keep_prob is None else tf.nn.dropout(weights, keep_prob))
				scale, beta = tf.Variable(tf.ones([output_dim])), tf.Variable(tf.zeros([output_dim]))
				p_mean      = tf.Variable(tf.zeros([output_dim]), trainable=False)
				p_var       = tf.Variable(tf.ones([output_dim]),  trainable=False)
				if is_training:
					decay         = 0.999
					b_mean, b_var = tf.nn.moments(z_bn, [0])
					t_mean        = tf.assign(p_mean, p_mean*decay + b_mean*(1.0-decay))
					t_var         = tf.assign(p_var, p_var*decay + b_var*(1.0-decay))
					with tf.control_dependencies([t_mean, t_var]):
						preactivate = tf.nn.batch_normalization(
							x=z_bn, mean=b_mean, variance=b_var, offset=beta, scale=scale, variance_epsilon=epsilon)

				else:
					preactivate = tf.nn.batch_normalization(
						x=z_bn, mean=p_mean, variance=p_var, offset=beta, scale=scale, variance_epsilon=epsilon)

			else:
				# without batch normalization
				preactivate = tf.matmul(input_tensor, weights if keep_prob is None else tf.nn.dropout(weights, keep_prob)) + biases

			tf.summary.histogram(name=layer_name + '/pre_activations', values=preactivate)

		activations = act(preactivate, 'activation')
		tf.summary.histogram(name=layer_name + '/activations', values=activations)
		return activations

def build_network(layers, layers_act=None, input_x_holder=None, input_y_holder=None, dropout_holder=None,
				  init_weights=None, bn=False, is_training=True, scope=''):
	with tf.name_scope('{}input'.format(scope)):
		x  = tf.placeholder(tf.float32, shape=[None, layers[ 0]], name='x-input') if input_x_holder is None else input_x_holder
		y_ = tf.placeholder(tf.float32, shape=[None, layers[-1]], name='y-input') if input_y_holder is None else input_y_holder

	lp = x

	for i in range(1, len(layers)-1):
		lp = nn_layer(input_tensor=lp, input_dim=layers[i-1], output_dim=layers[i], layer_name='{}layer{}'.format(scope, i),
					  init_weights=init_weights, keep_prob=dropout_holder, bn=bn, is_training=is_training,
					  act=layers_act[i] if layers_act else tf.nn.tanh)

	y = nn_layer(input_tensor=lp, input_dim=layers[i], output_dim=layers[i+1], layer_name='{}layer_full'.format(scope),
				 init_weights=init_weights, keep_prob=dropout_holder, bn=bn, is_training=is_training,
				 act=layers_act[i+1] if layers_act else tf.nn.tanh)

	return x, y_, y
