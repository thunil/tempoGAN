#!/usr/bin/python2

# ----------------------------------------------------------------------------
#
# MantaFlow fluid solver framework
# Copyright 2017 Kiwon Um, Nils Thuerey
#
# This program is free software, distributed under the terms of the
# Apache License, Version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Training code for ML-FLIP
#
# ----------------------------------------------------------------------------

import sys, os, argparse
if(len(sys.argv)<2): print("Call with list of data directories, at least one. E.g., 'python tf_train.py ../data/manta-flip/training_data ' \n"); exit(1); 

import numpy as np
import tensorflow as tf

import tf_datasets, tf_network

import logging, pickle

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

parser = argparse.ArgumentParser(description='Generate Training Data', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-o', '--output', default='../data/mlflip-tf/',        help='output directory')
parser.add_argument('-s', '--steps',  default=10000, type=int,   help='maximum training steps')
parser.add_argument('-b', '--batch',  default=100,   type=int,   help='batch size for one step training')
parser.add_argument('-d', '--dnet',   default='27-34-2',         help='detection networks int-int-...')
parser.add_argument('-m', '--mnet',   default='27-34-2',         help='modification networks int-int-...')
parser.add_argument(      '--dact',   default='none-tanh-tanh',  help='activation function for detection networks')
parser.add_argument(      '--mact',   default='none-tanh-tanh',  help='activation function for modification networks')
parser.add_argument('-v', '--mve',    action="store_true",       help='turn on mean-variance learning')
parser.add_argument(      '--bn',     action="store_true",       help='turn on batch normalization')
parser.add_argument(      '--nosmax', action="store_true",       help='do not use the softmax model')
parser.add_argument('-r', '--decay',  default=0.1,   type=float, help='regularization coefficient')
parser.add_argument(      '--ddrop',  default=0.1,   type=float, help='dropout rate (detection)')
parser.add_argument(      '--mdrop',  default=0.1,   type=float, help='dropout rate (modification)')
parser.add_argument('datadirs', action="store", nargs="+",       help='path(s) to the training data (e.g., ../data/manta-flip/training_data)')
pargs = parser.parse_args()

pargs.output = os.path.normpath(pargs.output)
os.path.isdir(pargs.output) or os.makedirs(pargs.output)
logging.basicConfig(filename='{}/training-info.log'.format(pargs.output), level=logging.INFO)
logging.getLogger().addHandler(logging.StreamHandler())

with open(pargs.output+'/run_args.pickle', 'wb') as f: pickle.dump(vars(pargs), f)
with open(pargs.output+'/run_cmd.txt', 'w') as f: f.write(' '.join(os.uname()) + '\n' + ' '.join(sys.argv))

data_sets, N_tuple = tf_datasets.read_data_sets(dirs=sorted(pargs.datadirs), use_softmax=(not pargs.nosmax))
scale = {}
for i in data_sets.train.get_data(): scale[i] = max(abs(data_sets.train.get_data()[i].min()), abs(data_sets.train.get_data()[i].max()))
with open(pargs.output+'/scale.pickle', 'wb') as f: pickle.dump(scale, f)
logging.info(pargs)
logging.info('{} tuples have been loaded; randomly selected {} for the training set and {} for the test set'.format(
	N_tuple, data_sets.train._num_examples, data_sets.test._num_examples))
logging.info(scale)

# statistics
with PdfPages(pargs.output+'/histogram.pdf') as pdf:
	l = data_sets.train.get_data()['labels']
	if not pargs.nosmax: l = l[:,0]
	for i in data_sets.train.get_data():
		d = data_sets.train.get_data()[i][(l==1).reshape(-1)]
		for j in range(d.shape[1]):
			plt.figure()
			plt.hist(d[:,j].reshape(-1), bins='auto')
			plt.title('Histogram of {}[{}]'.format(i, j))
			plt.savefig(pdf, format='pdf')
			plt.close()

sess = tf.InteractiveSession()

################################################################################
# neural networks
logging.info('Neural network structure: detection {} and modification {}'.format(pargs.dnet, pargs.mnet))
dlayers    = list(map(int, pargs.dnet.split('-')))
mlayers    = list(map(int, pargs.mnet.split('-')))
dact       = list(map(tf_network.parse_act, pargs.dact.split('-')))
mact       = list(map(tf_network.parse_act, pargs.mact.split('-')))
init_w     = {'w': {'stddev': 0.1}, 'b': {'value': 0.5}}
x          = tf.placeholder(tf.float32, shape=[None, dlayers[0]], name='x-input')
keep_prob  = tf.placeholder(tf.float32, name='keep_prob_detector') if pargs.ddrop>0.0 else None
keep_prob2 = tf.placeholder(tf.float32, name='keep_prob_modifier') if pargs.mdrop>0.0 else None
y_,  y     = tf_network.build_network(dlayers, dact, init_weights=init_w, input_x_holder=x, dropout_holder=keep_prob,  bn=pargs.bn, scope='detector/')[1:]
y2_, y2    = tf_network.build_network(mlayers, mact, init_weights=init_w, input_x_holder=x, dropout_holder=keep_prob2, bn=pargs.bn, scope='modifier/')[1:]
if pargs.mve:
	s      = tf_network.build_network(mlayers, mact, init_weights=init_w, input_x_holder=x, input_y_holder=y2_, dropout_holder=keep_prob2, bn=pargs.bn, scope='modifier_var/')[2]

################################################################################
# evaluation functions
log_dict = {}
with tf.name_scope('accuracy'):
	with tf.name_scope('correct_prediction'):
		if pargs.nosmax: corr, appx = tf.cast(tf.less(y_, 0.5), tf.int64), tf.cast(tf.less(y, 0.5), tf.int64) # f: splashing, t: non-splashing
		else:            corr, appx = tf.argmax(y_, 1), tf.argmax(y, 1)                                       # 0: splashing, 1: non-splashing

		correct_prediction = tf.equal(corr, appx)
		accuracy           = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))
		y_b_appx           = tf.equal(appx, 0) # true: splashing, false: non-splashing
		y_b_corr           = tf.equal(corr, 0)

		N_corr_non_splas    = tf.reduce_sum(corr)
		N_corr_splashing    = tf.cast(tf.shape(y_)[0], tf.int64) - N_corr_non_splas
		diff_appx_corr      = tf.logical_xor(y_b_appx, y_b_corr)
		N_corr_spl_appx_non = tf.reduce_sum(tf.cast(tf.logical_and(y_b_corr, diff_appx_corr), tf.float32))
		N_corr_non_appx_spl = tf.reduce_sum(tf.cast(tf.logical_and(y_b_appx, diff_appx_corr), tf.float32))
		false_negative      = N_corr_spl_appx_non/tf.cast(N_corr_splashing, tf.float32)
		false_positive      = N_corr_non_appx_spl/tf.cast(N_corr_non_splas, tf.float32)

		log_dict['accuracy']                     = accuracy
		log_dict['false_negative_corr_T_appx_F'] = false_negative
		log_dict['false_positive_corr_F_appx_T'] = false_positive
		log_dict['splashes/corr']                = 1.0 - tf.reduce_mean(tf.cast(corr, tf.float32))
		log_dict['splashes/appx']                = 1.0 - tf.reduce_mean(tf.cast(appx, tf.float32))

	with tf.name_scope('loss'):
		loss_normalizer = 1.0/tf.cast(tf.shape(y2_)[0], tf.float32)

		with tf.name_scope('detector'):
			if pargs.nosmax: loss_detector = tf.nn.l2_loss(y - y_)
			else:            loss_detector = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(logits=y, labels=y_))
			log_dict['detector/loss'] = loss_detector*loss_normalizer if pargs.nosmax else loss_detector

		with tf.name_scope('modifier'):
			loss_modifier = tf.nn.l2_loss(y2 - y2_)
			log_dict['modifier/loss'] = loss_modifier*loss_normalizer
			if pargs.mve:
				loss_modifier_mve = 0.5*tf.reduce_sum(((y2 - y2_)**2)/(s**2 + 1e-4)) + 0.5*tf.reduce_sum(tf.log(s**2 + 1e-4)) # mean variance estimate
				log_dict['modifier_mve/loss'] = loss_modifier_mve*loss_normalizer

		loss = loss_detector + loss_modifier
		log_dict['sum_loss'] = log_dict['detector/loss'] + log_dict['modifier/loss']
		if pargs.mve:
			loss_mve = loss_detector + loss_modifier_mve
			log_dict['sum_loss_mve'] = log_dict['detector/loss'] + log_dict['modifier_mve/loss']

		if pargs.decay>0.0:
			w_detector     = tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES, "detector/")
			w_modifier     = tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES, "modifier/")
			decay_detector = tf.add_n([tf.nn.l2_loss(v) for v in w_detector])*pargs.decay
			decay_modifier = tf.add_n([tf.nn.l2_loss(v) for v in w_modifier])*pargs.decay

			loss += decay_detector + decay_modifier

			log_dict['detector/decay'] = decay_detector
			log_dict['modifier/decay'] = decay_modifier
			log_dict['sum_loss'] += decay_detector + decay_modifier

			if pargs.mve:
				w_modifier_var = tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES, "modifier_var/")
				decay_modifier_var = tf.add_n([tf.nn.l2_loss(v) for v in w_modifier_var])*pargs.decay
				loss_mve += decay_modifier_var
				log_dict['modifier_var/decay'] = decay_modifier_var
				log_dict['sum_loss_mve'] += decay_modifier_var

	for i in log_dict:
		tf.summary.scalar(name=i, tensor=log_dict[i])

with tf.name_scope('train'):
	train_step = tf.train.AdamOptimizer(1e-4).minimize(loss)
	if pargs.mve: train_step_mve = tf.train.AdamOptimizer(1e-4).minimize(loss_mve)

merged = tf.summary.merge_all()
train_writer = tf.summary.FileWriter(pargs.output + '/summary/train', sess.graph)
test_writer  = tf.summary.FileWriter(pargs.output + '/summary/test')

saver = tf.train.Saver()
tf.global_variables_initializer().run()

################################################################################
# train
batch_test = data_sets.test.get_data()
feed_data_test = { x: batch_test['inputs'], y_: batch_test['labels'], y2_: batch_test['modvel']/scale['modvel'] }
if pargs.ddrop>0.0: feed_data_test.update({keep_prob : 1.0})
if pargs.mdrop>0.0: feed_data_test.update({keep_prob2: 1.0})

losses_key = sorted([k for k in log_dict if 'loss' in k])
losses_fetch = [log_dict[k] for k in losses_key]
for i in range(pargs.steps):
	if (i%10==0):           # test
		summary, acc, loss_i = sess.run([merged, accuracy, losses_fetch], feed_dict=feed_data_test)
		test_writer.add_summary(summary, i)
		print('Step {}/{}: accuracy={:.5f}, {}'.format(i,pargs.steps, acc, ', '.join('{}={:.5f}'.format(*t) for t in zip(losses_key, loss_i))))

	else:                   # train
		batch = data_sets.train.next_batch(pargs.batch)

		fetches = [merged, train_step]
		feed_data = { x: batch['inputs'], y_: batch['labels'], y2_: batch['modvel']/scale['modvel'] }

		if pargs.mve and i>int(pargs.steps/2):
			fetches = [merged, train_step_mve]
			if pargs.mdrop>0.0: y2_appx = sess.run(y2, feed_dict={ x: batch['inputs'], keep_prob2: 1.0 })
			else:               y2_appx = sess.run(y2, feed_dict={ x: batch['inputs'] })
			feed_data = { x: batch['inputs'], y_: batch['labels'], y2_: batch['modvel']/scale['modvel'], y2: y2_appx }

		if pargs.ddrop>0.0: feed_data.update({keep_prob : 1.0-pargs.ddrop})
		if pargs.mdrop>0.0: feed_data.update({keep_prob2: 1.0-pargs.mdrop})

		params_sess_run = dict(fetches=fetches, feed_dict=feed_data)
		if (i%100==99):
			run_metadata = tf.RunMetadata()
			params_sess_run.update(options=tf.RunOptions(trace_level=tf.RunOptions.FULL_TRACE), run_metadata=run_metadata)

		summary = sess.run(**params_sess_run)[0]

		if (i%100==99):
			print('Adding run metadata for step {}'.format(i))
			train_writer.add_run_metadata(run_metadata, 'step{:03d}'.format(i))

		train_writer.add_summary(summary, i)

# the last test
summary, acc_summary = sess.run([merged, list(log_dict.values())], feed_dict=feed_data_test)
test_writer.add_summary(summary, i)

train_writer.close()
test_writer.close()

logging.info(dict(zip(log_dict.keys(), acc_summary)))

# save the trained model
model_file = saver.save(sess, '{}/model.ckpt'.format(pargs.output))
logging.info('Trained model saved to {}'.format(model_file))

sess.close()
