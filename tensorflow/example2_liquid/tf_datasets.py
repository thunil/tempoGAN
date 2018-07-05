# ----------------------------------------------------------------------------
#
# MantaFlow fluid solver framework
# Copyright 2017 Kiwon Um, Nils Thuerey
#
# This program is free software, distributed under the terms of the
# Apache License, Version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Utilities for loading training/test data of TensorFlow
#
# ----------------------------------------------------------------------------

import os, glob
import numpy as np

class DataSet(object):

	def __init__(self, data):

		self._num_examples = data['labels'].shape[0]
		self._data = data

		self._epochs_completed = 0
		self._index_in_epoch = 0

	def get_data(self):
		return self._data

	def next_batch(self, batch_size):
		start = self._index_in_epoch
		self._index_in_epoch += batch_size
		if self._index_in_epoch > self._num_examples:
			# finished epoch
			self._epochs_completed += 1
			# shuffle the data
			perm = np.arange(self._data['labels'].shape[0])
			np.random.shuffle(perm)
			for ikey in self._data: self._data[ikey] = self._data[ikey][perm]
			# start next epoch
			start = 0
			self._index_in_epoch = batch_size
			assert batch_size <= self._num_examples

		end = self._index_in_epoch
		batch = {k: self._data[k][start:end] for k in self._data}
		return batch

def read_data_sets(dirs, use_softmax=True):
	class DataSets(object):
		pass
	data_sets = DataSets()

	np.random.seed(1)

	files_0, files_1 = [], []
	for i in dirs:
		files_0.extend(sorted(glob.glob('{}/*_p0.npz'.format(i))))
		files_1.extend(sorted(glob.glob('{}/*_p1.npz'.format(i))))

	set_0 = {}
	for ifile in files_0:
		data = np.load(ifile)
		for ikey in data:
			set_0[ikey] = np.concatenate((set_0[ikey], data[ikey]), axis=0) if ikey in set_0 else data[ikey]

	set_1 = {}
	for ifile in files_1:
		data = np.load(ifile)
		for ikey in data:
			set_1[ikey] = np.concatenate((set_1[ikey], data[ikey]), axis=0) if ikey in set_1 else data[ikey]

	print(format(dirs))
	if len(set_0)==0 or len(set_1)==0:
		print("Error - no data directories found!")
		exit(1)
	if set_0['labels'].shape[0]==0 or set_1['labels'].shape[0]==0:
		print("Error - data directories conainted no data!")
		exit(1)

	larger_set = set_0 if (set_0['labels'].shape[0]>set_1['labels'].shape[0]) else set_1

	min_size = min(set_0['labels'].shape[0], set_1['labels'].shape[0])
	larger_set = set_0 if (set_0['labels'].shape[0]>set_1['labels'].shape[0]) else set_1
	perm = np.arange(larger_set['labels'].shape[0])
	np.random.shuffle(perm)
	for ikey in larger_set:
		larger_set[ikey] = larger_set[ikey][perm]
		larger_set[ikey] = larger_set[ikey][-min_size:]

	dataset = {}
	perm = np.arange(set_0['labels'].shape[0] + set_1['labels'].shape[0])
	np.random.shuffle(perm)
	for ikey in set_0:
		arr = np.concatenate((set_0[ikey], set_1[ikey]), axis=0)
		arr[0::2], arr[1::2] = set_0[ikey], set_1[ikey]
		dataset[ikey] = arr[perm]

	if use_softmax:
		# softmax-model; use two values
		nlabels = np.empty([dataset['labels'].shape[0], 2])
		nlabels[:,0] = dataset['labels'].reshape((dataset['labels'].shape[0]))
		nlabels[:,1] = np.logical_not(nlabels[:,0])
		dataset['labels'] = nlabels

	set_train, set_test, test_size = {}, {}, int(dataset['labels'].shape[0]*0.25)
	for ikey in dataset:
		set_train[ikey] = dataset[ikey][:-test_size]
		set_test[ikey]  = dataset[ikey][-test_size:]

	data_sets.train = DataSet(set_train)
	data_sets.test  = DataSet(set_test)

	return data_sets, dataset['labels'].shape[0]
