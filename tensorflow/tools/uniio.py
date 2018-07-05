#******************************************************************************
#
# MantaFlow fluid solver framework
# Copyright 2017 Nils Thuerey, Boris Bonev
#
# This program is free software, distributed under the terms of the
# Apache License, Version 2.0 
# http://www.apache.org/licenses/LICENSE-2.0
#
# Read mantaflow uni files into numpy arrays
# note - only supports 3D grids for now
# (python2 , switch to python3 below)
#
#******************************************************************************

import gzip
import struct
import sys
import os
import shutil
from datetime import date
from collections import namedtuple
import numpy as np

PY3K = sys.version_info >= (3, 0)

# read content of grid
def RU_read_content(bytestream, header):
	assert (header['bytesPerElement'] == 12 and header['elementType'] == 2) or (header['bytesPerElement'] == 4 and (header['elementType'] == 0 or header['elementType'] == 1))

	if (header['elementType'] == 0):
		data = np.frombuffer(bytestream.read(), dtype="int32") # int grid
	else:
		data = np.frombuffer(bytestream.read(), dtype="float32") # float grid , scalar or vec3
	
	channels = 1
	if (header['elementType'] == 2):
		channels = 3

	dimensions = [header['dimT'], header['dimZ'], header['dimY'], header['dimX'], channels]
	if header['dimT']<=1:
		dimensions = [header['dimZ'], header['dimY'], header['dimX'], channels]

	return data.reshape( *dimensions, order='C')

# read uni file header (v3)
def RU_read_header(bytestream):
	ID = bytestream.read(4)
	# in python3, ID == b'MNT3' or b'MNT2' or ..., have to decode
	if(PY3K): ID = ID.decode("utf-8") 
	if ID=="MNT2":
		# unpack header struct object
		header = namedtuple('HeaderV3', 'dimX, dimY, dimZ, gridType, elementType, bytesPerElement, info, timestamp')
		# convert to namedtuple and then directly to a dict
		header = header._asdict(header._make(struct.unpack('iiiiii256sQ', bytestream.read(288))))

		# when writing, we'll need a v4 header field, re-pack...
		header['dimT'] = 0 
		header['info'] = header['info'][0:252]
		head4 = namedtuple('HeaderV4', 'dimX, dimY, dimZ, gridType, elementType, bytesPerElement, info, dimT, timestamp')(**header)
		header = head4._asdict()

	elif ID=="MNT3":
		# unpack header struct object
		header = namedtuple('HeaderV4', 'dimX, dimY, dimZ, gridType, elementType, bytesPerElement, info, dimT, timestamp')
		# convert to namedtuple and then directly to a dict
		# header is shorter for v3!
		header = header._asdict(header._make(struct.unpack('iiiiii252siQ', bytestream.read(288))))

	elif ID=="M4T2" or ID=="M4T3":
		print("read_header error - 4D grids not yet supported")
		exit(1)

	else:
		print("read_header error - unknown header '%s' " % ID)
		exit(1)
	
	return header

# use this to read the .uni file. It will return the header as dictionary and the content as np-array
def readUni(filename):
	#print("Reading '%s'" % filename) # debug
	with gzip.open(filename, 'rb') as bytestream:
		header = RU_read_header(bytestream)
		content = RU_read_content(bytestream, header)
		#print("Strides "+format(content.strides))

		return header, content

# use this to write a .uni file. The header has to be supplied in the same dictionary format as the output of readuni
def writeUni(filename, header, content):
	#print("Writing '%s'" % filename) # debug
	#print("Strides "+format(content.strides))
	with gzip.open(filename, 'wb') as bytestream:

		# write the header of the uni file (old v3 header)
		#bytestream.write(b'MNT2') # v3
		#head_tuple = namedtuple('GenericDict', header.keys())(**header)
		#head_buffer = struct.pack('iiiiii256sQ', *head_tuple)

		# current header
		bytestream.write(b'MNT3') # new, v4
		head_tuple = namedtuple('HeaderV4', header.keys())(**header)
		head_buffer = struct.pack('iiiiii252siQ', *head_tuple)
		bytestream.write(head_buffer)

		# always convert to single precision floats
		if content.dtype!="float32":
			content = np.asarray(content, dtype="float32") 

		# write grid content
		if (header['elementType'] == 2):
			# vec3 grid
			content = content.reshape(header['dimX']*header['dimY']*header['dimZ']*3, order='C')
		else:
			# int or scalar grid
			content = content.reshape(header['dimX']*header['dimY']*header['dimZ'], order='C')

		if sys.version_info >= (3,0):
			# changed for Python3
			bytestream.write(memoryview(content))
		else:
			bytestream.write(np.getbuffer(content))

# backup code to test folder
def backupFile(name, test_path):
	code_path = os.path.dirname(name) + '/' + os.path.basename(name)
	if len(os.path.dirname(name))==0:
		code_path = ".%s" % code_path
	shutil.copy(code_path, test_path + os.path.basename(name))

#******************************************************************************
# particle data

def RP_read_header(bytestream):
    ID = bytestream.read(4)     # NOTE: useless
    # unpack header struct object
    head = namedtuple('UniPartHeader', 'dim, dimX, dimY, dimZ, elementType, bytesPerElement, info, timestamp')
    # convert to namedtuple and then directly to a dict
    head = head._asdict(head._make(struct.unpack('iiiiii256sQ', bytestream.read(288))))

    return head

def RP_read_content(bytestream, head, data_type=None): # data_type = {None: BasicParticleSystem; "float32": Real; "int32": Int}
    assert(head['bytesPerElement']==16 or head['bytesPerElement']==12 or head['bytesPerElement']==4)

    if(head['elementType']==0): # BasicParticleSystem
        print('(BasicParticleSystem) ' )
        data = np.frombuffer(bytestream.read(), dtype=np.dtype([('f1',(np.float32,3)),('f2',(np.int32,1))]))['f1']
    else:                       # head['elementType']==1: ParticleDataImpl<T>, where T = {float32: Real(4) or Vec3(12); int32: Int(4)}
        print('(ParticleDataImpl<T={}{}>) '.format(data_type, 'x3' if (head['bytesPerElement']==12) else '') )
        data = np.reshape(np.frombuffer(bytestream.read(), dtype=data_type), (-1, 3 if (head['bytesPerElement']==12) else 1))

    return data

def readParticles(filename, data_type=None):
    print('Reading {} ... '.format(filename) )
    with gzip.open(filename, 'rb') as bytestream:
        head = RP_read_header(bytestream)
        data = RP_read_content(bytestream, head, data_type)

        print('Done.')
        return head, data

#******************************************************************************
# numpy array files

npBuf = {} # store arrays
npCnt = {} # filename counter
# FIXME , todo - add byte size limit per file at some point, to prevent them from getting too large

# buffer arrays, and write multiple to single file
def writeNumpyBuf(filename, content):
	global npBuf,npCnt
	if not filename in npBuf:
		npBuf[filename] = []
		npCnt[filename] = 0
	npBuf[filename].append(content)
	#print("writing buffered, arrays "+format( len(npBuf[filename]) ) + ", size "+ format(content.size) )
	if len(npBuf[filename])>10:
		#print("writing buffered "+filename)
		np.savez_compressed( filename+("_%04d.npz"%(npCnt[filename])), *npBuf[filename] )
		npCnt[filename] += 1
		npBuf[filename] = []

# write all remaining ones
def finalizeNumpyBufs():
	global npBuf,npCnt
	for filename in npBuf.keys():
		if len(npBuf[filename])>0:
			#print("writing last buffered "+filename+ ", left " + format(len(npBuf[filename])))
			np.savez_compressed( filename+("_%04d.npz"%(npCnt[filename])), *npBuf[filename] )
	# reset...
	npBuf = {} 
	npCnt = {} 
		

# write a single numpy array into an npz file
def writeNumpySingle(filename, content):
	#print("writing "+filename)
	np.savez_compressed( filename, content )

def readNumpy(filename):
	#print("reading "+filename)
	npz = np.load( filename )
	return npz


