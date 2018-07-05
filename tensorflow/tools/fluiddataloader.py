#******************************************************************************
#
# tempoGAN: A Temporally Coherent, Volumetric GAN for Super-resolution Fluid Flow
# Copyright 2018 You Xie, Erik Franz, Mengyu Chu, Nils Thuerey
#
# flexible fluid data loader
#
# This program is free software, distributed under the terms of the
# Apache License, Version 2.0 
# http://www.apache.org/licenses/LICENSE-2.0
#
#******************************************************************************

import os, glob, re, math, threading 
import numpy as np
import scipy.ndimage 

# necessary for loading uni files , could easily be disabled if necessary
import uniio 
# data type for new python arrays
FDG_DTYPE = np.float32


class FluidDataLoader(object):
	""" Fluid Data Loader - load npz files from a collection of directories

		label/GT data can be passed in different ways: an array with 1 data per dir, 
		a filename for npz data, or a generic function called for each loaded input
	"""

	def __init__(self, print_info=1, base_path="../data/", simdirname="sim_%04d/", indices=[], 
				filename=None, filename_index_min=0, filename_index_max=200, wildcard=None,
				array_y=None, filename_y=None, func_y=None, data_fraction=1.,
				shape=None, shape_y=None, collapse_z=False, shuffle_on_load=False,
				multi_file_list=None, multi_file_list_y=None, multi_file_idxOff=None, multi_file_idxOff_y=None,
				postproc_func=None, postproc_func_y=None,
				np_load_string=None , np_load_string_y=None , oldNamingScheme=False):
		""" Constructor , arguments:
			print_info: debugging info , <=0 off, 1 some info, 2 full
			base_path: path prefix for all sim dirs
			simdirname: sim directory name with printf placeholders (eg %04d) for indices
			indices: list of directory numbers to load
			filename: filename with printf placeholders (eg %04d) for numbered input data x
					  typical example string: "density_%04d.npz"
			          currently uni and npz files are supported
			filename_index_min: start index for filenames, controls index range (min to max)
			filename_index_max: maximal index for filenames 
			wildcard: optional, use specified wildcard for regexp matching filenames in sim dir for x.
			          has to contain a group (...) to extract ID string for y, group(1) used by default.
			          note, if wildcard string is given, this overrides filename with index range
			multi_file_list: list of file name prefixes, if given, the loader will load and concatenate
					all correspnding files
			multi_file_list_y: " analogous for y
			multi_file_idxOff: list of file index offsets for files multi_file_list
					can be used to load files with different index into same data entry
			multi_file_idxOff_y: " analogous for y
			postproc_func: function to be called for every data sample, can be used to post-process
					data in a custom way
			postproc_func_y: " analogous for y
			array_y: optional, label data as array, 1 entry per dir
			filename_y: optional, filenames for label data; needs # placeholder if used with wildcard
			func_y: optional, labelling func, called with sim index, filename and file index for every loeaded entry
			data_fraction: don't load all files, but also a fraction of it
			shape, shape_y: target shapes for x,y data; input shapes are determined from files to load; 
				warning - can easily rescale & interpolate channels target channels dont match source channels.
				note - only used for scaling. no reshapes; if no shape is given, pass through data unmodified
			collapse_z: remove Z-axis for 2d data sets, ie, transform ZYXc to YXc when z==1
			            also removes Z-component for pure velocity 2d data sets, ie, changes c=3 to c=2; TODO , make separate switch? (see removeZComponent)
			shuffle_on_load: randomize order on load? definitely not recommended for validation set
				note: better use flow(... , shuffle=True) in most cases
			np_load_string: by default, load "arr_0" from npz files, change using this string
			np_load_string_y: same as np_load_string but for loading y data; if none is given, 
				the string for x is used (np_load_string)
			oldNamingScheme: revert to old scheme with double indices for dir & filename
				by default the loader expects: data/sim_XXXX/density_low_YYYY.sth
				the old naming scheme was: data/sim_XXXX/frame_YYYY/density_low_XXXX_YYYY.sth
			
		"""
		# path basics
		self.base_path = base_path
		self.simdirname = simdirname
		self.indices = indices

		# x data files
		self.filename = filename
		self.filename_index_min = filename_index_min
		self.filename_index_max = filename_index_max
		self.wildcard = wildcard

		self.multi_file_list   = multi_file_list  
		self.multi_file_list_y = multi_file_list_y
		self.multi_file_idxOff   = multi_file_idxOff  
		self.multi_file_idxOff_y = multi_file_idxOff_y
		self.postproc_func   = postproc_func  
		self.postproc_func_y = postproc_func_y

		# y data for labeling x
		self.filename_y = filename_y
		self.array_y = array_y
		self.func_y = func_y

		# further options
		self.data_fraction = data_fraction
		self.shape = shape
		self.shape_y = shape_y
		self.collapse_z = collapse_z
		self.shuffle_on_load = shuffle_on_load

		# initialize npz load
		if np_load_string is not None:
			self.np_load_string = np_load_string
		else:
			self.np_load_string = "arr_0" 
		# input data
		inCnt = 0
		if self.filename is not None: inCnt += 1
		if self.wildcard is not None: inCnt += 1
		if inCnt>1:
			# sanity check - only one of those allowed
			raise FluidDataLoaderError("FluidDataLoader error: for input data loading, only specify one of: input filename, or wildcard")
		# label data
		inCnt = 0
		if self.filename_y is not None: inCnt += 1
		if self.array_y is not None: inCnt += 1
		if self.func_y is not None: inCnt += 1
		if inCnt>1:
			# sanity check - only one of those allowed
			raise FluidDataLoaderError("FluidDataLoader error:  for label data loading, only specify one of: input filename, array or function")

		if np_load_string_y is not None:
			self.np_load_string_y = np_load_string_y
		else:
			self.np_load_string_y = self.np_load_string

		self.print_info = print_info
		if self.print_info:
			print("FluidDataLoader init, path %s, filename %s" % (self.base_path,self.filename) )
		self.oldNamingScheme = oldNamingScheme

		# sanity check file lists 
		if self.multi_file_idxOff is not None and self.multi_file_list is not None:
			if len(self.multi_file_list) != len(self.multi_file_idxOff):
				raise FluidDataLoaderError("FluidDataLoader error: multi file list and idxOff lists have to match " + format([len(self.multi_file_list) , len(self.multi_file_idxOff)]) )
		if self.multi_file_idxOff_y is not None and self.multi_file_list_y is not None:
			if len(self.multi_file_list_y) != len(self.multi_file_idxOff_y):
				raise FluidDataLoaderError("FluidDataLoader error: multi file list and idxOff lists for y have to match " + format([len(self.multi_file_list_y) , len(self.multi_file_idxOff_y)]) )

		# all initialized upon load:
		self.x = None
		self.y = None
		self.xfn = None
		self.have_y_npz = False # does y contain numpy array data? 

		self.loadDirs()
		self.printStats()

	def getFilename(self, sim_index, fnbase, frame_index): 
		if not self.oldNamingScheme:
			fn = os.path.join( self.base_path, os.path.join((self.simdirname % sim_index), (fnbase % frame_index)) )
		else:
			# both parts simdir & file have both indices!
			fn = os.path.join( self.base_path, os.path.join( (self.simdirname % (sim_index,frame_index)), (fnbase % (sim_index,frame_index) ) ))
		return fn

	def collectFilenamesFromDir(self, list_index): 
		""" Build filename list from single dir
			list_index: number in index list (or alternatively label list)
		"""

		sim_index = self.indices[list_index] # get simulation directory index from list
		labelstr = "" # debug info only
		foundCnt = 0

		if self.wildcard is not None:
			search_dir = os.path.join( self.base_path, (self.simdirname % sim_index) )
			os.chdir(search_dir)
			allFiles = [f for f in glob.glob("*") if os.path.isfile(f)] # list all files
			files = []
			for f in allFiles:
				match = re.search(self.wildcard, f) # note, matched again below...
				if match:
					files.append(f)
	
			if len(files)<1:
				raise FluidDataLoaderError("Error - no files found in directory '%s' with wildcard '%s' " %(search_dir, self.wildcard) )

			files = sorted(files) # sort by name

			n  = max(1, int(len(files)*self.data_fraction))
			tf = float(len(files))/n # spread over time range (eg 200 frames) 
			fcnt = 0
			for t in range(0,n):
				filelist_index = int(t*tf) # full range
				fn = files[filelist_index]
				self.xfn.append(os.path.join(search_dir, fn))
				foundCnt += 1

				# construct label, closely follows index version below
				if self.filename_y is not None:
					mx = re.search(self.wildcard, fn) 
					listy = self.filename_y.split("$")
					if(len(listy)!=2):
						raise FluidDataLoaderError("Error - when using a wildcard for x, filename_y needs to contain exactly one '$' where the file id string from x will be inserted to build the filename for y. Current, invalid, filename_y is '%s' " %(self.filename_y) )
					fny = listy[0] + mx.group(1) + listy[1]
	
					if not os.path.isfile(fny): # make sure file for y exists
						raise FluidDataLoaderError("Error - y file '%s' for x file '%s' doesnt exist in search dir '%s' " %(fny, fn, search_dir ) )

					fny = os.path.join(search_dir, fny) 
					self.yfn.append(fny)
					self.have_y_npz = True # flag to indicate we have np arrays in y

				if self.array_y is not None:
					if self.y is None:
						self.y = []
					self.y.append( self.array_y[list_index] )
					labelstr = " with label " + format( self.array_y[list_index] )

				if self.func_y is not None:
					print("NYI! test...")

		else:
			# "simple" index range 
			n  = max(1, int((self.filename_index_max-self.filename_index_min)*self.data_fraction))
			tf = float(self.filename_index_max-self.filename_index_min)/n
			for t in range(0,n):
				filelist_index = int(self.filename_index_min + t*tf) # full range

				fn = self.getFilename(sim_index, self.filename, filelist_index)
				self.xfn.append(fn)
				foundCnt += 1

				if self.filename_y is not None:
					fny = self.getFilename(sim_index, self.filename_y, filelist_index)
					self.yfn.append(fny)
					self.have_y_npz = True # flag to indicate we have np arrays in y

				if self.array_y is not None:
					if self.y is None:
						self.y = []
					self.y.append( self.array_y[list_index] )
					labelstr = " with label " + format( self.array_y[list_index] )

				if self.func_y is not None:
					print("NYI! test...")
					if self.y is None:
						self.y = []
					self.y.append( self.func_y(list_index, sim_index, t, fn) )

		if self.print_info:
			print("Found " +format(foundCnt) +" files from sim ID "+format(sim_index) + labelstr )


	def getDim(self,shape):
		""" small helper to compute dimensionality of data from shape 
		"""
		dim = -1 
		if len(shape)==4: # probably ZYXc
			if(shape[0]==1): 
				dim = 2
			else:
				dim = 3
		if len(shape)==5: # probably 4d, TZYXc
			dim = 4
		#print("Dim "+format(dim)+ " for " + format(shape) )
		return dim

	def removeZComponent(self,x):
		""" Optional, and 2D only: remove Z entry from 3d vec fields
		"""
		if not self.collapse_z: return x
		if not self.getDim( x.shape )==2: return x
		if not x.shape[3]==3: return x  # only apply for pure velocity grids with 3 channels
		x2d = np.zeros( (1,x.shape[1],x.shape[2],2), dtype=FDG_DTYPE )
		x2d[:,:,:,0] = x[:,:,:,0] # only keep x,y
		x2d[:,:,:,1] = x[:,:,:,1]
		return x2d

	def mogrifyFilenameIndex(self, fn, idxOffset):
		""" Parse, determine index, and change
		"""
		match = re.search("(.*_)([\d]+)\.([\w]+)", fn)  # split into groups: path/name_ , %04d , ext
		if match:
			if len(match.groups())!=3:
				raise FluidDataLoaderError("FluidDataLoader error: got filename %s, but could not fully split up into name,4-digit and extension " % (fn))
			#print "A "  + format(match.groups())
			idx = int(match.group(2))
			idx = max(self.filename_index_min, min(self.filename_index_max-1, idx+idxOffset) )
			#print "A "  + format(match.group(2))
			fn = "%s%04d.%s" % (match.group(1), idx, match.group(3))
			#print "fn "  + fn
		else:
			raise FluidDataLoaderError("FluidDataLoader error: got filename %s, but could not split up into name,4-digit and extension " % (fn))
			
		#exit()
		# density1_([\w]+).npz
		return fn

	def loadSingleDatum(self, fn, lstr, idxOffset=0):
		""" Determine file type and load
		"""
		if idxOffset!=0:
			fn = self.mogrifyFilenameIndex(fn,idxOffset)
		if self.print_info>1:
			print("Loading: "+fn+", "+lstr)
		# detect file type
		if fn.endswith( ".npz" ):
			ar = np.load(fn)[ lstr ]
		elif fn.endswith( ".uni" ):
			_, ar = uniio.readUni(fn) # load-string lstr not needed for uni files
			#ar = ar[::-1] # make a copy of the array in reverse order
		else:
			raise FluidDataLoaderError("FluidDataLoader error: got filename %s, but only .uni or .npz supported at the moment " % (fn))

		return ar

	def loadFiles(self):
		""" Load all NPZs from list.
			Note, data always has to have shape ZYXc (3d) or YXc (2d), 
			where c is channels (eg 3 for vels, 1 for scalar data).
		""" 
		n = len(self.xfn)
		for t in range(n):
			fof = 0 if self.multi_file_idxOff is None else self.multi_file_idxOff[0]
			fx = self.loadSingleDatum(self.xfn[t], self.np_load_string , fof) 

			if self.multi_file_list is not None:
				# concat multiple files...
				basename = self.xfn[t]
				if not basename.find(self.multi_file_list[0])>=0:
					raise FluidDataLoaderError("Error, input filename '%s' doesnt contain given string '%s'"%(basename,self.multi_file_list[0]))
				for i in range(1,len(self.multi_file_list)):
					fnr = basename.replace(self.multi_file_list[0] , self.multi_file_list[i])
					fof = 0 if self.multi_file_idxOff is None else self.multi_file_idxOff[i]
					_fx = self.loadSingleDatum(fnr, self.np_load_string , fof ) 
					fx = np.append( fx, _fx , axis=len(fx.shape)-1 )

			# apply post-processing function (if given)
			if self.postproc_func is not None:
				fx = self.postproc_func(fx, self)

			# ... and the same again for y
			if self.have_y_npz:
				fofy = 0 if self.multi_file_idxOff_y is None else self.multi_file_idxOff_y[0]
				fy = self.loadSingleDatum(self.yfn[t], self.np_load_string_y , fofy ) 

				if self.multi_file_list_y is not None:
					basename = self.yfn[t]
					if not basename.find(self.multi_file_list_y[0])>=0:
						raise FluidDataLoaderError("Error, input filename y '%s' doesnt contain given string '%s'"%(basename,self.multi_file_list_y[0]))
					for i in range(1,len(self.multi_file_list_y)):
						fnr = basename.replace(self.multi_file_list_y[0] , self.multi_file_list_y[i])
						fofy = 0 if self.multi_file_idxOff_y is None else self.multi_file_idxOff_y[i]
						_fy = self.loadSingleDatum(fnr, self.np_load_string_y , fofy ) 
						fy = np.append( fy, _fy , axis=len(fy.shape)-1 )

				if self.postproc_func_y is not None:
					fy = self.postproc_func_y(fy, self)

			fx = self.removeZComponent(fx) # optional!

			# intialize x/y arrays upon first use
			if self.x is None:
				self.data_shape = fx.shape

				if self.shape is None: # no target shape? use data res
					self.shape = fx.shape
					self.do_zoom = False
				else:
					self.do_zoom = True
					self.zoom_shape = []
					for i in range(len(self.shape)):
						self.zoom_shape.append( float(self.shape[i]) / self.data_shape[i] )
					#if self.collapse_z and self.dim==2: self.zoom_shape[ len(self.zoom_shape)-1 ] = 1. # old, dont zoom channels
					if self.print_info: print("Zoom for x by "+format(self.zoom_shape) )

				if self.print_info: print("Allocating x data for "+format(n)+" entries of size "+format(self.shape) )
				self.x = np.zeros( tuple([n]+list(self.shape)) , dtype=FDG_DTYPE )

			# optional zoom, is initialized with original array
			if self.do_zoom:
				fx = scipy.ndimage.zoom( fx, self.zoom_shape, order=1 ) 

			# finally store t-th data sample
			self.x[t,:]  = fx

			# and again for y ...
			if self.have_y_npz:
				fy = self.removeZComponent(fy)

				if self.y is None:
					self.data_shape_y = fy.shape
					if self.shape_y is None: # no target shape? use data res
						self.shape_y = fy.shape
						self.do_zoom = False
					else:
						self.do_zoom = True
						self.zoom_shape_y = []
						for i in range(len(self.shape_y)):
							self.zoom_shape_y.append( float(self.shape_y[i]) / self.data_shape_y[i] )
						if self.print_info: print("Zoom for y by "+format(self.zoom_shape_y) )

					if self.print_info: print("Allocating y data for "+format(n)+" entries of size "+format(self.shape_y) )
					self.y = np.zeros( tuple([n]+list(self.shape_y)) , dtype=FDG_DTYPE )

				if self.do_zoom:
					fy = scipy.ndimage.zoom( fy, self.zoom_shape_y, order=1 ) 

				self.y[t,:]  = fy

			if self.print_info and t==0: print("loadFiles: data size x "+ format(self.x.shape) + ((", y " + format(self.y.shape)) if self.filename_y is not None else "") ) 

		# x (and optionally y) arrays complete now, retrieve with get() later on



	def loadDirs(self):
		""" Main load function: collect all files in multiple directories,
			and load the necessary fraction; potentially rescale (zoom) data, if enabled
		"""
		self.xfn = []
		self.yfn = []
		currDir = os.getcwd()

		for i in range(len(self.indices)):
			self.collectFilenamesFromDir( i )
			os.chdir( currDir )
			
		# debug info, print full lists
		if self.print_info>1:
			#print("Full list x: "+format(self.xfn)) print("Full list y: "+format(self.yfn))
			print( "\nfilenames x:" ); print( ("\n".join(self.xfn)) )
			if self.filename_y is not None:
				print( "\nfilenames y:" ); print( ("\n".join(self.yfn)) )

		self.loadFiles()
		os.chdir( currDir )

		# remove z axis of all 3D data fields for whole data vector
		if self.collapse_z:
			if self.getDim(self.x[0].shape)==2:
				self.x = np.reshape( self.x, [self.x.shape[0], self.shape[1],self.shape[2],self.shape[3]] ) # remove z-axis for x
			if self.have_y_npz and self.getDim(self.y[0].shape)==2:
				self.y = np.reshape( self.y, [self.y.shape[0], self.shape_y[1],self.shape_y[2],self.shape_y[3]] )

		# do manual shuffling once (needs to reorder x,y and filenames for x,y)
		if self.shuffle_on_load:
			idxr = np.random.permutation(self.x.shape[0])
			self.x = self.x[idxr]
			if self.have_y_npz: self.y = self.y[idxr] # y is np array , reorder...

			xfn2,yfn2,y2 = [],[],[]
			for i in range(len(self.xfn)):
				xfn2.append( self.xfn[idxr[i]] )
				if not self.have_y_npz and self.y is not None: y2.append( self.y[idxr[i]] ) # non np array y
				if self.filename_y is not None: yfn2.append( self.yfn[idxr[i]] )
			self.xfn, self.yfn = xfn2,yfn2
			if not self.have_y_npz and self.y is not None: self.y = y2
		# loading done


	def arrayStats(self, values, weights=None):
		average = np.average(values) #, weights=weights)
		variance = np.average((values-average)**2) #, weights=weights)  # Fast and numerically precise
		return (average, math.sqrt(variance))

	def perChannelStats(self, values, info=None):
		if values.shape[-1]>1:
			if info:
				print(format(info))
			for c in range(values.shape[-1]):
				print("\t\t"+format(c)+": "+format(self.arrayStats(values[...,c]) ))

	def printStats(self):
		""" General info about loaded data sets """
		if self.print_info:
			print("Loaded "+format(self.x.shape[0])+" datasets" + (", shuffled" if self.shuffle_on_load else "") )
			print("\tData shape x " + format(self.x.shape))
			m,s = self.arrayStats(self.x)
			print("\tx mean & var: " + format([m,s]))
			if m<1e-10 and s<1e-10: # sanity check, any non-zero values?
				raise FluidDataLoaderError("FluidDataLoader error: aborting, input data x is all zero")
			self.perChannelStats(self.x, "\tPer channel mean & var for x: ")
			if self.have_y_npz: 
				print("\tData shape y " + format(self.y.shape))
				m,s = self.arrayStats(self.y)
				print("\tmean & var for y: " + format([m,s]))
				if m<1e-10 and s<1e-10: # sanity check, any non-zero values?
					raise FluidDataLoaderError("FluidDataLoader error: aborting, input data y is all zero")

	def get(self):
		""" After loading, return arrays 
		"""
		return self.x , self.y , self.xfn

	def getFullInfo(self):
		""" Summarize full data set as string
		"""
		ret = ""
		printMean = True
		for i in range(len(self.xfn)):
			ret = ret + ("%d/%d, file %s, shape %s" % (i, len(self.xfn), self.xfn[i], format(self.x[i].shape) ))
			if printMean:
				ret = ret + (", x mean %s " % (format(np.mean(self.x[i])) ))
			if self.filename_y is not None:
				ret = ret + (", file_y %s " % (self.yfn[i]) )
			if self.have_y_npz:
				ret = ret + (", shape_y %s " % (format(self.y[i].shape)) )
				if printMean:
					ret = ret + (", y mean %s " % (format(np.mean(self.y[i])) ))
			if self.array_y is not None:
				ret = ret + (", y %s " % (format(self.y[i])) )
			ret = ret + "\n"
		return ret

class FluidDataLoaderError(Exception):
	''' FDL errors '''


