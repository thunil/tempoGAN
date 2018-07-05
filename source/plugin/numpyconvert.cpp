/******************************************************************************
*
* MantaFlow fluid solver framework
* Copyright 2017 Steffen Wiewel, Moritz Becher 
*
* This program is free software, distributed under the terms of the
* Apache License, Version 2.0 
* http://www.apache.org/licenses/LICENSE-2.0
*
* Plugins to convert mantaflow grids to/from numpy arrays, also support pdata fields  
# (only compiled if NUMPY is enabled)
*
******************************************************************************/

#include "manta.h"
#include "kernel.h"
#include "grid.h"
#include "particle.h"
#include "levelset.h"

using namespace std;
namespace Manta
{

//====================================================================================================
// Grid numpy conversion
//----------------------------------------------------------------------------------------------------

template<typename T>
void copyArrayToGridScalar(const PyArrayContainer source, T& target)
{
	target.setConst(0.0f);
	unsigned int uGridSize = target.getSizeX() * target.getSizeY() * target.getSizeZ();
	assertMsg(source.TotalSize == uGridSize, "The size of the numpy array doesn't match the size of the Grid!");
	
	NumpyTypes eDataType  = source.DataType; 

	switch (eDataType)
	{
		case NumpyTypes::N_FLOAT:
			FOR_IDX(target) { target(idx) = (reinterpret_cast<float*>(source.pData))[idx]; }
			break;
		case NumpyTypes::N_DOUBLE:
			FOR_IDX(target) { target(idx) = (reinterpret_cast<double*>(source.pData))[idx]; } 
			break;
		default:
			errMsg("unknown/unsupported type of Numpy array");
			return;
	}
}

template<typename T>
void copyGridToArrayScalar(const T& source, PyArrayContainer target)
{
	unsigned int uGridsize = source.getSizeX() * source.getSizeY() * source.getSizeZ();
	assertMsg(target.TotalSize == uGridsize, "The size of the numpy array doesn't match the size of the grid!");
	
	NumpyTypes eDataType = target.DataType;

	switch (eDataType)
	{
		case NumpyTypes::N_FLOAT:
			FOR_IDX(source) { reinterpret_cast<float*>(target.pData)[idx] = source(idx); }
			break;
		case NumpyTypes::N_DOUBLE:
			FOR_IDX(source) { reinterpret_cast<double*>(target.pData)[idx] = source(idx); }
			break;
		default:
			errMsg("unknown/unsupported type of Numpy array");
			break;
	}
}

template<typename T>
void copyArrayToGridVector(const PyArrayContainer source, T& target)
{
	unsigned int uSizeX = target.getSizeX();
	unsigned int uSizeY = target.getSizeY();
	unsigned int uSizeZ = target.getSizeZ();
	unsigned int uSizeW = 3u;
	
	assertMsg(source.TotalSize == uSizeX * uSizeY * uSizeZ * uSizeW, "The size of the numpy array doesn't match the size of the grid!");
	
	NumpyTypes eDataType = source.DataType;

	switch (eDataType)
	{
		case NumpyTypes::N_FLOAT:
			FOR_IDX(target) { for(int w = 0; w < 3; ++w) { target(idx)[w] = (reinterpret_cast<float*>(source.pData))[idx*3+w]; } }
			break;
		case NumpyTypes::N_DOUBLE:
			FOR_IDX(target) { for(int w = 0; w < 3; ++w) { target(idx)[w] = (reinterpret_cast<double*>(source.pData))[idx*3+w]; } }
			break;
		default:
			errMsg("unknown/unsupported type of Vec3 Numpy array");
			break;
	}
}

template<typename T>
void copyGridToArrayVector(const T& source, PyArrayContainer target)
{
	unsigned int uSizeX = source.getSizeX();
	unsigned int uSizeY = source.getSizeY();
	unsigned int uSizeZ = source.getSizeZ();
	unsigned int uSizeW = 3u;

	assertMsg(target.TotalSize == uSizeX * uSizeY * uSizeZ * uSizeW, "The size of the numpy array doesn't match the size of the grid!");
	
	NumpyTypes eDataType = target.DataType;
	
	switch (eDataType)
	{
		case NumpyTypes::N_FLOAT:
			FOR_IDX(source) { for(int w = 0; w < 3; ++w) { (reinterpret_cast<float*>(target.pData))[idx*3+w] = source(idx)[w]; } }
			break;
		case NumpyTypes::N_DOUBLE:
			FOR_IDX(source) { for(int w = 0; w < 3; ++w) { (reinterpret_cast<double*>(target.pData))[idx*3+w] = source(idx)[w]; } }
			break;
		default:
			errMsg("unknown/unsupported type of Vec3 Numpy array");
			break;
	}
}

//====================================================================================================
// Python interface
//----------------------------------------------------------------------------------------------------

PYTHON() void copyArrayToGridReal(const PyArrayContainer source, Grid<Real>& target) {
	copyArrayToGridScalar<Grid<Real>>(source, target);
}

PYTHON() void copyGridToArrayReal(const Grid<Real>& source, PyArrayContainer target) {
	copyGridToArrayScalar<Grid<Real>>(source, target);
}

PYTHON() void copyArrayToGridLevelset(const PyArrayContainer source, LevelsetGrid& target) {
	copyArrayToGridScalar<LevelsetGrid>(source, target);
}

PYTHON() void copyGridToArrayLevelset(const LevelsetGrid& source, PyArrayContainer target) {
	copyGridToArrayScalar<LevelsetGrid>(source, target);
}

PYTHON() void copyArrayToGridVec3(const PyArrayContainer source, Grid<Vec3>& target) {
	copyArrayToGridVector<Grid<Vec3>>(source, target);
}

PYTHON() void copyGridToArrayVec3(const Grid<Vec3>& source, PyArrayContainer target) {
	copyGridToArrayVector<Grid<Vec3>>(source, target);
}

PYTHON() void copyArrayToGridMAC(const PyArrayContainer source, MACGrid& target) {
	copyArrayToGridVector<MACGrid>(source, target);
}

PYTHON() void copyGridToArrayMAC(const MACGrid& source, PyArrayContainer target) {
	copyGridToArrayVector<MACGrid>(source, target);
}

//====================================================================================================
// pdata conversion functions
//----------------------------------------------------------------------------------------------------

template<typename T>
void numpyToParticleDataImpl(const PyArrayContainer source, ParticleDataImpl<T> &target) {
	assertMsg(source.TotalSize == target.size(), "The size of the numpy array doesn't match the size of the pdata field!");
	std::copy(reinterpret_cast<const T*>(source.pData), reinterpret_cast<const T*>(source.pData)+source.TotalSize,  &(target[0]));
}
template<typename T>
void particleDataImplToNumpy(const ParticleDataImpl<T> &source, PyArrayContainer target) {
	assertMsg(target.TotalSize == source.size(), "The size of the numpy array doesn't match the size of the pdata field!");
	std::copy(&(source[0]), &(source[0])+target.TotalSize, reinterpret_cast<T*>(target.pData));
}

// python interface

PYTHON() void copyArrayToPdataInt(const PyArrayContainer source, ParticleDataImpl<int> &target) { 
	numpyToParticleDataImpl<int>(source, target); 
}
PYTHON() void copyPdataToArrayInt(const ParticleDataImpl<int> &source, PyArrayContainer target) { 
	particleDataImplToNumpy<int>(source, target); 
}

PYTHON() void copyArrayToPdataReal(const PyArrayContainer source, ParticleDataImpl<Real> &target) { 
	numpyToParticleDataImpl<Real>(source, target); 
}
PYTHON() void copyPdataToArrayReal(const ParticleDataImpl<Real> &source, PyArrayContainer target) { 
	particleDataImplToNumpy<Real>(source, target); 
}

PYTHON() void copyArrayToPdataVec3(const PyArrayContainer source, ParticleDataImpl<Vec3> &target) {
	assertMsg(source.TotalSize == target.size()*3, "The size of the numpy array doesn't match the size of the pdata field!");
	std::copy(reinterpret_cast<const Real*>(source.pData), reinterpret_cast<const Real*>(source.pData)+source.TotalSize,  &(target[0][0]));
}
PYTHON() void copyPdataToArrayVec3(const ParticleDataImpl<Vec3> &source, PyArrayContainer target) {
	assertMsg(target.TotalSize == source.size()*3, "The size of the numpy array doesn't match the size of the pdata field!");
	std::copy(&(source[0][0]), &(source[0][0])+target.TotalSize, reinterpret_cast<Real*>(target.pData));
}

} // manta

