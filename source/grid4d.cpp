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

#include <limits>
#include <sstream>
#include <cstring>

#include "grid4d.h"
#include "levelset.h"
#include "kernel.h"
#include "mantaio.h" 

using namespace std;
namespace Manta {


//******************************************************************************
// GridBase members

Grid4dBase::Grid4dBase (FluidSolver* parent) 
	: PbClass(parent), mType(TypeNone)
{
	checkParent();
}


//******************************************************************************
// Grid4d<T> members

// helpers to set type
template<class T> inline Grid4dBase::Grid4dType typeList()        { return Grid4dBase::TypeNone; }
template<>        inline Grid4dBase::Grid4dType typeList<Real>()  { return Grid4dBase::TypeReal; }
template<>        inline Grid4dBase::Grid4dType typeList<int>()   { return Grid4dBase::TypeInt;  }
template<>        inline Grid4dBase::Grid4dType typeList<Vec3>()  { return Grid4dBase::TypeVec3; }
template<>        inline Grid4dBase::Grid4dType typeList<Vec4>()  { return Grid4dBase::TypeVec4; }


template<class T>
Grid4d<T>::Grid4d(FluidSolver* parent, bool show)
	: Grid4dBase(parent)
{
	assertMsg( parent->is3D() && parent->supports4D(), "To use 4d grids create a 3d solver with fourthDim>0");

	mType = typeList<T>();
	Vec3i s = parent->getGridSize();
	mSize = Vec4i(s.x, s.y, s.z, parent->getFourthDim() ); 
	mData = parent->getGrid4dPointer<T>();
	assertMsg( mData, "Couldnt allocate data pointer!");
	
	mStrideZ = (mSize.x * mSize.y);
	mStrideT = (mStrideZ * mSize.z);

	Real sizemax = (Real)mSize[0];
	for(int c=1; c<3; ++c) if(mSize[c]>sizemax) sizemax = mSize[c];
	// note - the 4d component is ignored for dx! keep same scaling as for 3d...
	mDx = 1.0 / sizemax;

	clear();
	setHidden(!show);
}

template<class T>
Grid4d<T>::Grid4d(const Grid4d<T>& a) : Grid4dBase(a.getParent()) {
	mSize = a.mSize;
	mType = a.mType;
	mStrideZ = a.mStrideZ;
	mStrideT = a.mStrideT;
	mDx = a.mDx;
	FluidSolver *gp = a.getParent();
	mData = gp->getGrid4dPointer<T>();
	assertMsg( mData, "Couldnt allocate data pointer!");

	memcpy(mData, a.mData, sizeof(T) * a.mSize.x * a.mSize.y * a.mSize.z * a.mSize.t);
}

template<class T>
Grid4d<T>::~Grid4d() {
	mParent->freeGrid4dPointer<T>(mData);
}

template<class T>
void Grid4d<T>::clear() {
	memset(mData, 0, sizeof(T) * mSize.x * mSize.y * mSize.z * mSize.t);
}

template<class T>
void Grid4d<T>::swap(Grid4d<T>& other) {
	if (other.getSizeX() != getSizeX() || other.getSizeY() != getSizeY() || other.getSizeZ() != getSizeZ() || other.getSizeT() != getSizeT())
		errMsg("Grid4d::swap(): Grid4d dimensions mismatch.");
	
	T* dswap = other.mData;
	other.mData = mData;
	mData = dswap;
}

template<class T>
void Grid4d<T>::load(string name) {
	if (name.find_last_of('.') == string::npos)
		errMsg("file '" + name + "' does not have an extension");
	string ext = name.substr(name.find_last_of('.'));
	if (ext == ".uni")
		readGrid4dUni(name, this);
	else if (ext == ".raw")
		readGrid4dRaw(name, this);
	else
		errMsg("file '" + name +"' filetype not supported");
}

template<class T>
void Grid4d<T>::save(string name) {
	if (name.find_last_of('.') == string::npos)
		errMsg("file '" + name + "' does not have an extension");
	string ext = name.substr(name.find_last_of('.'));
	if (ext == ".uni")
		writeGrid4dUni(name, this);
	else if (ext == ".raw")
		writeGrid4dRaw(name, this);
	else
		errMsg("file '" + name +"' filetype not supported");
}

//******************************************************************************
// Grid4d<T> operators

//! Kernel: Compute min value of Real Grid4d
KERNEL(idx, reduce=min) returns(Real minVal=std::numeric_limits<Real>::max())
Real kn4dMinReal(Grid4d<Real>& val) {
	if (val[idx] < minVal)
		minVal = val[idx];
}

//! Kernel: Compute max value of Real Grid4d
KERNEL(idx, reduce=max) returns(Real maxVal=-std::numeric_limits<Real>::max())
Real kn4dMaxReal(Grid4d<Real>& val) {
	if (val[idx] > maxVal)
		maxVal = val[idx];
}

//! Kernel: Compute min value of int Grid4d
KERNEL(idx, reduce=min) returns(int minVal=std::numeric_limits<int>::max())
int kn4dMinInt(Grid4d<int>& val) {
	if (val[idx] < minVal)
		minVal = val[idx];
}

//! Kernel: Compute max value of int Grid4d
KERNEL(idx, reduce=max) returns(int maxVal=std::numeric_limits<int>::min())
int kn4dMaxInt(Grid4d<int>& val) {
	if (val[idx] > maxVal)
		maxVal = val[idx];
}

//! Kernel: Compute min norm of vec Grid4d
KERNEL(idx, reduce=min) returns(Real minVal=std::numeric_limits<Real>::max())
template<class VEC> Real kn4dMinVec(Grid4d<VEC>& val) {
	const Real s = normSquare(val[idx]);
	if (s < minVal)
		minVal = s;
}

//! Kernel: Compute max norm of vec Grid4d
KERNEL(idx, reduce=max) returns(Real maxVal=-std::numeric_limits<Real>::max())
template<class VEC> Real kn4dMaxVec(Grid4d<VEC>& val) {
	const Real s = normSquare(val[idx]);
	if (s > maxVal)
		maxVal = s;
}


template<class T> Grid4d<T>& Grid4d<T>::safeDivide (const Grid4d<T>& a) {
	Grid4dSafeDiv<T> (*this, a);
	return *this;
}
template<class T> Grid4d<T>& Grid4d<T>::copyFrom (const Grid4d<T>& a, bool copyType ) {
	assertMsg (a.mSize.x == mSize.x && a.mSize.y == mSize.y && a.mSize.z == mSize.z && a.mSize.t == mSize.t, "different Grid4d resolutions "<<a.mSize<<" vs "<<this->mSize );
	memcpy(mData, a.mData, sizeof(T) * mSize.x * mSize.y * mSize.z * mSize.t);
	if(copyType) mType = a.mType; // copy type marker
	return *this;
}
/*template<class T> Grid4d<T>& Grid4d<T>::operator= (const Grid4d<T>& a) {
	note: do not use , use copyFrom instead
}*/

KERNEL(idx) template<class T> void kn4dSetConstReal (Grid4d<T>& me, T val) { me[idx]  = val; }
KERNEL(idx) template<class T> void kn4dAddConstReal (Grid4d<T>& me, T val) { me[idx] += val; }
KERNEL(idx) template<class T> void kn4dMultConst (Grid4d<T>& me, T val) { me[idx] *= val; }
KERNEL(idx) template<class T> void kn4dClamp (Grid4d<T>& me, T min, T max) { me[idx] = clamp( me[idx], min, max); }

template<class T> void Grid4d<T>::add(const Grid4d<T>& a) {
	Grid4dAdd<T,T>(*this, a);
}
template<class T> void Grid4d<T>::sub(const Grid4d<T>& a) {
	Grid4dSub<T,T>(*this, a);
}
template<class T> void Grid4d<T>::addScaled(const Grid4d<T>& a, const T& factor) { 
	Grid4dScaledAdd<T,T> (*this, a, factor); 
}
template<class T> void Grid4d<T>::setConst(T a) {
	kn4dSetConstReal<T>( *this, T(a) );
}
template<class T> void Grid4d<T>::addConst(T a) {
	kn4dAddConstReal<T>( *this, T(a) );
}
template<class T> void Grid4d<T>::multConst(T a) {
	kn4dMultConst<T>( *this, a );
}

template<class T> void Grid4d<T>::mult(const Grid4d<T>& a) {
	Grid4dMult<T,T> (*this, a);
}

template<class T> void Grid4d<T>::clamp(Real min, Real max) {
	kn4dClamp<T> (*this, T(min), T(max) );
}

template<> Real Grid4d<Real>::getMax() {
	return kn4dMaxReal (*this);
}
template<> Real Grid4d<Real>::getMin() {
	return kn4dMinReal (*this);
}
template<> Real Grid4d<Real>::getMaxAbs() {
	Real amin = kn4dMinReal (*this);
	Real amax = kn4dMaxReal (*this);
	return max( fabs(amin), fabs(amax));
}
template<> Real Grid4d<Vec4>::getMax() {
	return sqrt(kn4dMaxVec<Vec4> (*this));
}
template<> Real Grid4d<Vec4>::getMin() { 
	return sqrt(kn4dMinVec<Vec4> (*this));
}
template<> Real Grid4d<Vec4>::getMaxAbs() {
	return sqrt(kn4dMaxVec<Vec4> (*this));
}
template<> Real Grid4d<int>::getMax() {
	return (Real) kn4dMaxInt (*this);
}
template<> Real Grid4d<int>::getMin() {
	return (Real) kn4dMinInt (*this);
}
template<> Real Grid4d<int>::getMaxAbs() {
	int amin = kn4dMinInt (*this);
	int amax = kn4dMaxInt (*this);
	return max( fabs((Real)amin), fabs((Real)amax));
}
template<> Real Grid4d<Vec3>::getMax() {
	return sqrt(kn4dMaxVec<Vec3> (*this));
}
template<> Real Grid4d<Vec3>::getMin() { 
	return sqrt(kn4dMinVec<Vec3> (*this));
}
template<> Real Grid4d<Vec3>::getMaxAbs() {
	return sqrt(kn4dMaxVec<Vec3> (*this));
}


template<class T> void Grid4d<T>::printGrid(int zSlice, int tSlice, bool printIndex, int bnd) {
	std::ostringstream out;
	out << std::endl;
	FOR_IJKT_BND(*this,bnd) {
		IndexInt idx = (*this).index(i,j,k,t);
		if (  ( (zSlice>=0 && k==zSlice) || (zSlice<0) ) &&
		  	  ( (tSlice>=0 && t==tSlice) || (tSlice<0) ) ) {
			out << " ";
			if(printIndex) out << "  "<<i<<","<<j<<","<<k<<","<<t <<":";
			out << (*this)[idx]; 
			if(i==(*this).getSizeX()-1 -bnd) {
				out << std::endl; 
				if(j==(*this).getSizeY()-1 -bnd) {
					out << std::endl; 
					if(k==(*this).getSizeZ()-1 -bnd) { out << std::endl; }
			} }
		}
	}
	out << endl; debMsg("Printing '" << this->getName() <<"' "<< out.str().c_str()<<" " , 1);
}


// helper to set/get components of vec4 Grids
KERNEL(idx) void knGetComp4d(const Grid4d<Vec4>& src, Grid4d<Real>& dst, int c) { dst[idx]    = src[idx][c]; };
KERNEL(idx) void knSetComp4d(const Grid4d<Real>& src, Grid4d<Vec4>& dst, int c) { dst[idx][c] = src[idx];    };
PYTHON() void getComp4d(const Grid4d<Vec4>& src, Grid4d<Real>& dst, int c) { knGetComp4d(src,dst,c); };
PYTHON() void setComp4d(const Grid4d<Real>& src, Grid4d<Vec4>& dst, int c) { knSetComp4d(src,dst,c); };


KERNEL(fourd) template<class T> void knSetBnd4d (Grid4d<T>& grid, T value, int w) { 
	bool bnd = 
		(i<=w || i>=grid.getSizeX()-1-w || 
		 j<=w || j>=grid.getSizeY()-1-w || 
		 k<=w || k>=grid.getSizeZ()-1-w ||
		 t<=w || t>=grid.getSizeT()-1-w );
	if (bnd) 
		grid(i,j,k,t) = value;
}

template<class T> void Grid4d<T>::setBound(T value, int boundaryWidth) {
	knSetBnd4d<T>( *this, value, boundaryWidth );
}

KERNEL(fourd) template<class T> void knSetBnd4dNeumann (Grid4d<T>& grid, int w) { 
	bool set = false;
	int  si=i, sj=j, sk=k, st=t;
	if( i<=w) {
		si = w+1; set=true;
	}
	if( i>=grid.getSizeX()-1-w){
		si = grid.getSizeX()-1-w-1; set=true;
	}
	if( j<=w){
		sj = w+1; set=true;
	}
	if( j>=grid.getSizeY()-1-w){
		sj = grid.getSizeY()-1-w-1; set=true;
	}
	if( k<=w ) {
		sk = w+1; set=true;
	}
	if( k>=grid.getSizeZ()-1-w ) {
		sk = grid.getSizeZ()-1-w-1; set=true;
	}
	if( t<=w ) {
		st = w+1; set=true;
	}
	if( t>=grid.getSizeT()-1-w ) {
		st = grid.getSizeT()-1-w-1; set=true;
	}
	if(set)
		grid(i,j,k,t) = grid(si, sj, sk, st);
}

template<class T> void Grid4d<T>::setBoundNeumann(int boundaryWidth) {
	knSetBnd4dNeumann<T>( *this, boundaryWidth );
}

//******************************************************************************
// testing helpers

//! compute maximal diference of two cells in the grid, needed for testing system
PYTHON() Real grid4dMaxDiff(Grid4d<Real>& g1, Grid4d<Real>& g2 )
{
	double maxVal = 0.;
	FOR_IJKT_BND(g1,0) {
		maxVal = std::max(maxVal, (double)fabs( g1(i,j,k,t)-g2(i,j,k,t) ));
	}
	return maxVal; 
}
PYTHON() Real grid4dMaxDiffInt(Grid4d<int>& g1, Grid4d<int>& g2 )
{
	double maxVal = 0.;
	FOR_IJKT_BND(g1,0) {
		maxVal = std::max(maxVal, (double)fabs( (double)g1(i,j,k,t)-g2(i,j,k,t) ));
	}
	return maxVal; 
}
PYTHON() Real grid4dMaxDiffVec3(Grid4d<Vec3>& g1, Grid4d<Vec3>& g2 )
{
	double maxVal = 0.;
	FOR_IJKT_BND(g1,0) {
		double d = 0.;
		for(int c=0; c<3; ++c) { 
			d += fabs( (double)g1(i,j,k,t)[c] - (double)g2(i,j,k,t)[c] );
		}
		maxVal = std::max(maxVal, d );
	}
	return maxVal; 
}
PYTHON() Real grid4dMaxDiffVec4(Grid4d<Vec4>& g1, Grid4d<Vec4>& g2 )
{
	double maxVal = 0.;
	FOR_IJKT_BND(g1,0) {
		double d = 0.;
		for(int c=0; c<4; ++c) { 
			d += fabs( (double)g1(i,j,k,t)[c] - (double)g2(i,j,k,t)[c] );
		}
		maxVal = std::max(maxVal, d );
	}
	return maxVal; 
}

// set a region to some value
KERNEL(fourd) template<class S>
void knSetRegion4d (Grid4d<S>& dst, Vec4 start, Vec4 end, S value )
{
	Vec4 p(i,j,k,t);
	for(int c=0; c<4; ++c) if(p[c]<start[c] || p[c]>end[c]) return;
	dst(i,j,k,t) = value;
}
//! simple init functions in 4d
PYTHON() void setRegion4d    (Grid4d<Real>& dst, Vec4 start, Vec4 end, Real value) { knSetRegion4d<Real>(dst,start,end,value); }
//! simple init functions in 4d, vec4
PYTHON() void setRegion4dVec4(Grid4d<Vec4>& dst, Vec4 start, Vec4 end, Vec4 value) { knSetRegion4d<Vec4>(dst,start,end,value); }

//! slow helper to visualize tests, get a 3d slice of a 4d grid
PYTHON() void getSliceFrom4d(Grid4d<Real>& src, int srct, Grid<Real>& dst) { 
	const int bnd = 0;
	if(! src.isInBounds(Vec4i(bnd,bnd,bnd,srct)) ) return;

	for(int k=bnd; k<src.getSizeZ()-bnd; k++) 
	for(int j=bnd; j<src.getSizeY()-bnd; j++) 
	for(int i=bnd; i<src.getSizeX()-bnd; i++)
	{
		if(!dst.isInBounds(Vec3i(i,j,k))) continue;
		dst(i,j,k) = src(i,j,k,srct);
	}
}
//! slow helper to visualize tests, get a 3d slice of a 4d vec4 grid
PYTHON() void getSliceFrom4dVec(Grid4d<Vec4>& src, int srct, Grid<Vec3>& dst, Grid<Real>* dstt=NULL) { 
	const int bnd = 0;
	if(! src.isInBounds(Vec4i(bnd,bnd,bnd,srct)) ) return;

	for(int k=bnd; k<src.getSizeZ()-bnd; k++) 
	for(int j=bnd; j<src.getSizeY()-bnd; j++) 
	for(int i=bnd; i<src.getSizeX()-bnd; i++)
	{
		if(!dst.isInBounds(Vec3i(i,j,k))) continue;
		for(int c=0; c<3; ++c) 
			dst(i,j,k)[c] = src(i,j,k,srct)[c];
		if(dstt) (*dstt)(i,j,k) = src(i,j,k,srct)[3];
	}
}


//******************************************************************************
// interpolation

//! same as in grid.h , but takes an additional optional "desired" size
static inline void gridFactor4d(Vec4 s1, Vec4 s2, Vec4 optSize, Vec4 scale, Vec4& srcFac, Vec4& retOff ) {
	for(int c=0; c<4; c++) { if(optSize[c] > 0.){ s2[c] = optSize[c]; } }
	srcFac = calcGridSizeFactor4d( s1, s2) / scale;
	retOff       = -retOff * srcFac + srcFac*0.5;
}

//! interpolate 4d grid from one size to another size
// real valued offsets & scale
KERNEL(fourd) template<class S>
void knInterpol4d(Grid4d<S>& target, Grid4d<S>& source, const Vec4& srcFac, const Vec4& offset)
{
	Vec4 pos = Vec4(i,j,k,t) * srcFac + offset;
	target(i,j,k,t) = source.getInterpolated(pos);
} 
//! linearly interpolate data of a 4d grid
PYTHON() void interpolateGrid4d( Grid4d<Real>& target, Grid4d<Real>& source , Vec4 offset=Vec4(0.), Vec4 scale=Vec4(1.), Vec4 size=Vec4(-1.) )
{
	Vec4 srcFac(1.), off2 = offset;
	gridFactor4d( toVec4(source.getSize()), toVec4(target.getSize()), size,scale,   srcFac,off2   );
	knInterpol4d<Real> (target, source, srcFac, off2 );
}
//! linearly interpolate vec4 data of a 4d grid
PYTHON() void interpolateGrid4dVec( Grid4d<Vec4>& target, Grid4d<Vec4>& source, Vec4 offset=Vec4(0.), Vec4 scale=Vec4(1.), Vec4 size=Vec4(-1.) )
{
	Vec4 srcFac(1.), off2 = offset;
	gridFactor4d( toVec4(source.getSize()), toVec4(target.getSize()), size,scale,   srcFac,off2   );
	knInterpol4d<Vec4> (target, source, srcFac, off2 );
}


// explicit instantiation
template class Grid4d<int>;
template class Grid4d<Real>;
template class Grid4d<Vec3>;
template class Grid4d<Vec4>;

} //namespace
