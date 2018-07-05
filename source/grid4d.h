/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Grid representation
 *
 ******************************************************************************/

#ifndef _GRID4D_H
#define _GRID4D_H

#include "manta.h"
#include "vectorbase.h"
#include "vector4d.h"
#include "kernel.h"


namespace Manta {

	
//! Base class for all grids
PYTHON() class Grid4dBase : public PbClass {
public:
	enum Grid4dType { TypeNone = 0, TypeReal = 1, TypeInt = 2, TypeVec3 = 4, TypeVec4 = 8 };
		
	PYTHON() Grid4dBase(FluidSolver* parent);
	
	//! Get the grids X dimension
	inline int getSizeX() const { return mSize.x; }
	//! Get the grids Y dimension
	inline int getSizeY() const { return mSize.y; }
	//! Get the grids Z dimension
	inline int getSizeZ() const { return mSize.z; }
	//! Get the grids T dimension
	inline int getSizeT() const { return mSize.t; }
	//! Get the grids dimensions
	inline Vec4i getSize() const { return mSize; }
	
	//! Get Stride in X dimension
	inline IndexInt getStrideX() const { return 1; }
	//! Get Stride in Y dimension
	inline IndexInt getStrideY() const { return mSize.x; }
	//! Get Stride in Z dimension
	inline IndexInt getStrideZ() const { return mStrideZ; }
	//! Get Stride in T dimension
	inline IndexInt getStrideT() const { return mStrideT; }
	
	inline Real getDx() { return mDx; }
	
	//! Check if indices are within bounds, otherwise error (should only be called when debugging)
	inline void checkIndex(int i, int j, int k, int t) const;
	//! Check if indices are within bounds, otherwise error (should only be called when debugging)
	inline void checkIndex(IndexInt idx) const;
	//! Check if index is within given boundaries
	inline bool isInBounds(const Vec4i& p, int bnd) const;
	//! Check if index is within given boundaries
	inline bool isInBounds(const Vec4i& p) const;
	//! Check if index is within given boundaries
	inline bool isInBounds(const Vec4& p, int bnd = 0) const { return isInBounds(toVec4i(p), bnd); }
	//! Check if linear index is in the range of the array
	inline bool isInBounds(IndexInt idx) const;
	
	//! Get the type of grid
	inline Grid4dType getType() const { return mType; }
	//! Check dimensionality
	inline bool is3D() const { return true; }
	inline bool is4D() const { return true; }

	//! 3d compatibility
	inline bool isInBounds(int i,int j, int k, int t, int bnd) const { return isInBounds( Vec4i(i,j,k,t), bnd ); }
	
	//! Get index into the data
	inline IndexInt index(int i, int j, int k, int t) const { DEBUG_ONLY(checkIndex(i,j,k,t)); return (IndexInt)i + (IndexInt)mSize.x * j + (IndexInt)mStrideZ * k + (IndexInt)mStrideT * t; }
	//! Get index into the data
	inline IndexInt index(const Vec4i& pos) const    { DEBUG_ONLY(checkIndex(pos.x,pos.y,pos.z,pos.t)); return (IndexInt)pos.x + (IndexInt)mSize.x * pos.y + (IndexInt)mStrideZ * pos.z + (IndexInt)mStrideT * pos.t; }
protected:
	
	Grid4dType mType;
	Vec4i      mSize;
	Real       mDx;
	// precomputed Z,T shift: to ensure 2D compatibility, always use this instead of sx*sy !
	IndexInt   mStrideZ; 
	IndexInt   mStrideT; 
};

//! Grid class
PYTHON() template<class T>
class Grid4d : public Grid4dBase {
public:
	//! init new grid, values are set to zero
	PYTHON() Grid4d(FluidSolver* parent, bool show = true);
	//! create new & copy content from another grid
	Grid4d(const Grid4d<T>& a);
	//! return memory to solver
	virtual ~Grid4d();
	
	typedef T BASETYPE;
	typedef Grid4dBase BASETYPE_GRID;
	
	PYTHON() void save(std::string name);
	PYTHON() void load(std::string name);
	
	//! set all cells to zero
	PYTHON() void clear();
	
	//! all kinds of access functions, use grid(), grid[] or grid.get()
	//! access data
	inline T get(int i,int j, int k, int t) const         { return mData[index(i,j,k,t)]; }
	//! access data
	inline T& get(int i,int j, int k, int t)              { return mData[index(i,j,k,t)]; }
	//! access data
	inline T get(IndexInt idx) const                           { DEBUG_ONLY(checkIndex(idx)); return mData[idx]; }
	//! access data
	inline T get(const Vec4i& pos) const                  { return mData[index(pos)]; }
	//! access data
	inline T& operator()(int i, int j, int k, int t)      { return mData[index(i, j, k,t)]; }
	//! access data
	inline T operator()(int i, int j, int k, int t) const { return mData[index(i, j, k,t)]; }
	//! access data
	inline T& operator()(IndexInt idx)                  { DEBUG_ONLY(checkIndex(idx)); return mData[idx]; }
	//! access data
	inline T operator()(IndexInt idx) const             { DEBUG_ONLY(checkIndex(idx)); return mData[idx]; }
	//! access data
	inline T& operator()(const Vec4i& pos)         { return mData[index(pos)]; }
	//! access data
	inline T operator()(const Vec4i& pos) const    { return mData[index(pos)]; }
	//! access data
	inline T& operator[](IndexInt idx)                  { DEBUG_ONLY(checkIndex(idx)); return mData[idx]; }
	//! access data
	inline const T operator[](IndexInt idx) const       { DEBUG_ONLY(checkIndex(idx)); return mData[idx]; }
	
	// interpolated access
	inline T    getInterpolated(const Vec4& pos) const { return interpol4d<T>(mData, mSize, mStrideZ, mStrideT, pos); }
	
	// assignment / copy

	//! warning - do not use "=" for grids in python, this copies the reference! not the grid content...
	//Grid4d<T>& operator=(const Grid4d<T>& a);
	//! copy content from other grid (use this one instead of operator= !)
	PYTHON() Grid4d<T>& copyFrom(const Grid4d<T>& a, bool copyType=true ); // old: { *this = a; }

	// helper functions to work with grids in scene files 

	//! add/subtract other grid
	PYTHON() void add(const Grid4d<T>& a);
	PYTHON() void sub(const Grid4d<T>& a);
	//! set all cells to constant value
	PYTHON() void setConst(T s);
	//! add constant to all grid cells
	PYTHON() void addConst(T s);
	//! add scaled other grid to current one (note, only "Real" factor, "T" type not supported here!)
	PYTHON() void addScaled(const Grid4d<T>& a, const T& factor); 
	//! multiply contents of grid
	PYTHON() void mult( const Grid4d<T>& a);
	//! multiply each cell by a constant scalar value
	PYTHON() void multConst(T s);
	//! clamp content to range (for vec3, clamps each component separately)
	PYTHON() void clamp(Real min, Real max);
	
	// common compound operators
	//! get absolute max value in grid 
	PYTHON() Real getMaxAbs();
	//! get max value in grid 
	PYTHON() Real getMax();
	//! get min value in grid 
	PYTHON() Real getMin();    
	//! set all boundary cells to constant value (Dirichlet)
	PYTHON() void setBound(T value, int boundaryWidth=1);
	//! set all boundary cells to last inner value (Neumann)
	PYTHON() void setBoundNeumann(int boundaryWidth=1);

	//! debugging helper, print grid from Python
	PYTHON() void printGrid(int zSlice=-1, int tSlice=-1,  bool printIndex=false, int bnd=0); 

	// c++ only operators
	template<class S> Grid4d<T>& operator+=(const Grid4d<S>& a);
	template<class S> Grid4d<T>& operator+=(const S& a);
	template<class S> Grid4d<T>& operator-=(const Grid4d<S>& a);
	template<class S> Grid4d<T>& operator-=(const S& a);
	template<class S> Grid4d<T>& operator*=(const Grid4d<S>& a);
	template<class S> Grid4d<T>& operator*=(const S& a);
	template<class S> Grid4d<T>& operator/=(const Grid4d<S>& a);
	template<class S> Grid4d<T>& operator/=(const S& a);
	Grid4d<T>& safeDivide(const Grid4d<T>& a);    
	
	//! Swap data with another grid (no actual data is moved)
	void swap(Grid4d<T>& other);

protected:
	T* mData;
};

// Python doesn't know about templates: explicit aliases needed
PYTHON() alias Grid4d<int>  Grid4Int;
PYTHON() alias Grid4d<Real> Grid4Real;
PYTHON() alias Grid4d<Vec3> Grid4Vec3;
PYTHON() alias Grid4d<Vec4> Grid4Vec4;


//! helper to compute grid conversion factor between local coordinates of two grids
inline Vec4 calcGridSizeFactor4d(Vec4i s1, Vec4i s2) {
	return Vec4( Real(s1[0])/s2[0], Real(s1[1])/s2[1], Real(s1[2])/s2[2] , Real(s1[3])/s2[3] );
}
inline Vec4 calcGridSizeFactor4d(Vec4 s1, Vec4 s2) {
	return Vec4( s1[0]/s2[0], s1[1]/s2[1], s1[2]/s2[2] , s1[3]/s2[3] );
}

// prototypes for grid plugins
void getComponent4d(const Grid4d<Vec4>& src, Grid4d<Real>& dst, int c);
void setComponent4d(const Grid4d<Real>& src, Grid4d<Vec4>& dst, int c);


//******************************************************************************
// Implementation of inline functions

inline void Grid4dBase::checkIndex(int i, int j, int k, int t) const {
	if ( i<0 || j<0  || i>=mSize.x || j>=mSize.y || k<0|| k>= mSize.z ||
         t<0|| t>= mSize.t ) {
		std::ostringstream s;
		s << "Grid4d " << mName << " dim " << mSize << " : index " << i << "," << j << "," << k << ","<<t<<" out of bound ";
		errMsg(s.str());
	}
}

inline void Grid4dBase::checkIndex(IndexInt idx) const {
	if (idx<0 || idx >= mSize.x * mSize.y * mSize.z * mSize.t) {
		std::ostringstream s;
		s << "Grid4d " << mName << " dim " << mSize << " : index " << idx << " out of bound ";
		errMsg(s.str());
	}
}

bool Grid4dBase::isInBounds(const Vec4i& p) const { 
	return (p.x >= 0 && p.y >= 0 && p.z >= 0 && p.t >= 0 && 
			p.x < mSize.x && p.y < mSize.y && p.z < mSize.z && p.t < mSize.t); 
}

bool Grid4dBase::isInBounds(const Vec4i& p, int bnd) const { 
	bool ret = (p.x >= bnd && p.y >= bnd && p.x < mSize.x-bnd && p.y < mSize.y-bnd);
	ret &= (p.z >= bnd && p.z < mSize.z-bnd); 
	ret &= (p.t >= bnd && p.t < mSize.t-bnd); 
	return ret;
}
//! Check if linear index is in the range of the array
bool Grid4dBase::isInBounds(IndexInt idx) const {
	if (idx<0 || idx >= mSize.x * mSize.y * mSize.z * mSize.t) {
		return false;
	}
	return true;
}

// note - ugly, mostly copied from normal GRID!

KERNEL(idx) template<class T, class S> void Grid4dAdd  (Grid4d<T>& me, const Grid4d<S>& other) { me[idx] += other[idx]; }
KERNEL(idx) template<class T, class S> void Grid4dSub  (Grid4d<T>& me, const Grid4d<S>& other) { me[idx] -= other[idx]; }
KERNEL(idx) template<class T, class S> void Grid4dMult (Grid4d<T>& me, const Grid4d<S>& other) { me[idx] *= other[idx]; }
KERNEL(idx) template<class T, class S> void Grid4dDiv  (Grid4d<T>& me, const Grid4d<S>& other) { me[idx] /= other[idx]; }
KERNEL(idx) template<class T, class S> void Grid4dAddScalar (Grid4d<T>& me, const S& other)  { me[idx] += other; }
KERNEL(idx) template<class T, class S> void Grid4dMultScalar(Grid4d<T>& me, const S& other)  { me[idx] *= other; }
KERNEL(idx) template<class T, class S> void Grid4dScaledAdd (Grid4d<T>& me, const Grid4d<T>& other, const S& factor) { me[idx] += factor * other[idx]; }

KERNEL(idx) template<class T> void Grid4dSafeDiv (Grid4d<T>& me, const Grid4d<T>& other) { me[idx] = safeDivide(me[idx], other[idx]); }
KERNEL(idx) template<class T> void Grid4dSetConst(Grid4d<T>& me, T value) { me[idx] = value; }

template<class T> template<class S> Grid4d<T>& Grid4d<T>::operator+= (const Grid4d<S>& a) {
	Grid4dAdd<T,S> (*this, a);
	return *this;
}
template<class T> template<class S> Grid4d<T>& Grid4d<T>::operator+= (const S& a) {
	Grid4dAddScalar<T,S> (*this, a);
	return *this;
}
template<class T> template<class S> Grid4d<T>& Grid4d<T>::operator-= (const Grid4d<S>& a) {
	Grid4dSub<T,S> (*this, a);
	return *this;
}
template<class T> template<class S> Grid4d<T>& Grid4d<T>::operator-= (const S& a) {
	Grid4dAddScalar<T,S> (*this, -a);
	return *this;
}
template<class T> template<class S> Grid4d<T>& Grid4d<T>::operator*= (const Grid4d<S>& a) {
	Grid4dMult<T,S> (*this, a);
	return *this;
}
template<class T> template<class S> Grid4d<T>& Grid4d<T>::operator*= (const S& a) {
	Grid4dMultScalar<T,S> (*this, a);
	return *this;
}
template<class T> template<class S> Grid4d<T>& Grid4d<T>::operator/= (const Grid4d<S>& a) {
	Grid4dDiv<T,S> (*this, a);
	return *this;
}
template<class T> template<class S> Grid4d<T>& Grid4d<T>::operator/= (const S& a) {
	S rez((S)1.0 / a);
	Grid4dMultScalar<T,S> (*this, rez);
	return *this;
}


//******************************************************************************
// Other helper functions

inline Vec4 getGradient4d(const Grid4d<Real>& data, int i, int j, int k, int t) {
	Vec4 v;
	if (i > data.getSizeX()-2) i= data.getSizeX()-2;
	if (j > data.getSizeY()-2) j= data.getSizeY()-2;
	if (k > data.getSizeZ()-2) k= data.getSizeZ()-2;
	if (t > data.getSizeT()-2) t= data.getSizeT()-2;
	if (i < 1) i = 1;
	if (j < 1) j = 1;
	if (k < 1) k = 1;
	if (t < 1) t = 1;
	v = Vec4( data(i+1,j  ,k  ,t  ) - data(i-1,j  ,k  ,t  ) ,
			  data(i  ,j+1,k  ,t  ) - data(i  ,j-1,k  ,t  ) , 
			  data(i  ,j  ,k+1,t  ) - data(i  ,j  ,k-1,t  ) , 
			  data(i  ,j  ,k  ,t+1) - data(i  ,j  ,k  ,t-1) );
	return v;
}


KERNEL(fourd) template<class S>
void KnInterpolateGrid4dTempl(Grid4d<S>& target, Grid4d<S>& source, const Vec4& sourceFactor , Vec4 offset) {
	Vec4 pos = Vec4(i,j,k,t) * sourceFactor + offset;
	if(!source.is3D()) pos[2] = 0.; // allow 2d -> 3d
	if(!source.is4D()) pos[3] = 0.; // allow 3d -> 4d
	target(i,j,k,t) = source.getInterpolated(pos);
} 

} //namespace
#endif
