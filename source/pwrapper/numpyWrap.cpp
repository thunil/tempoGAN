/******************************************************************************
*
* MantaFlow fluid solver framework
* Copyright 2017 Steffen Wiewel, Moritz Baecher, Rachel Chu
*
* This program is free software, distributed under the terms of the
* Apache License, Version 2.0 
* http://www.apache.org/licenses/LICENSE-2.0
*
* Convert mantaflow grids to/from numpy arrays
*
******************************************************************************/

#include "manta.h"
#include "pythonInclude.h"

namespace Manta {

PyMODINIT_FUNC initNumpy() { import_array(); }

template<> PyArrayContainer fromPy<PyArrayContainer>(PyObject* obj) {
	if (PyArray_API == NULL){
		initNumpy();
	}
	
	if (!PyArray_Check(obj)) {
		errMsg("argument is not an numpy array");
	}
	
	PyArrayContainer abuf;
	PyArrayObject* obj_p = reinterpret_cast<PyArrayObject*>(PyArray_CheckFromAny( obj, NULL, 0, 0, NPY_ARRAY_C_CONTIGUOUS | NPY_ARRAY_ENSUREARRAY | NPY_ARRAY_NOTSWAPPED, NULL));

	abuf.TotalSize = PyArray_SIZE(obj_p);
	int source_typ = PyArray_TYPE(obj_p);
	abuf.pData = PyArray_DATA(obj_p);
	
	switch (source_typ) {
		case NPY_FLOAT:
			abuf.DataType = N_FLOAT;
			break;
		case NPY_DOUBLE:
			abuf.DataType = N_DOUBLE;
			break;
		case NPY_INT:
			abuf.DataType = N_INT;
			break;
		default:
			errMsg("unknown type of Numpy array");
			break;
	}
	
	return abuf;
}

template<> PyArrayContainer* fromPyPtr<PyArrayContainer>(PyObject* obj, std::vector<void*>* tmp)
{
	if (!tmp) throw Error("dynamic de-ref not supported for this type");
	void* ptr = malloc(sizeof(PyArrayContainer));
	tmp->push_back(ptr);
	
	*((PyArrayContainer*) ptr) = fromPy<PyArrayContainer>(obj);
	return (PyArrayContainer*) ptr;
}

}
