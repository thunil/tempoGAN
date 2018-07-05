/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Functions for calculating wavelet turbulence,
 * plus helpers to compute vorticity, and strain rate magnitude 
 *
 ******************************************************************************/

#include "vectorbase.h"
#include "shapes.h"
#include "commonkernels.h"
#include "noisefield.h"

using namespace std;

namespace Manta {

//*****************************************************************************

// first some fairly generic interpolation functions for grids with multiple sizes

//! same as in grid.h , but takes an additional optional "desired" size
inline void calcGridSizeFactorMod(Vec3i s1, Vec3i s2, Vec3i optSize, Vec3 scale, Vec3& sourceFactor, Vec3& retOff ) {
	for(int c=0; c<3; c++) {
		if(optSize[c] > 0){ s2[c] = optSize[c]; }
	}
	sourceFactor = calcGridSizeFactor(s1,s2) / scale;
	retOff       = -retOff * sourceFactor + sourceFactor*0.5;
}

PYTHON() void interpolateGrid( Grid<Real>& target, Grid<Real>& source , Vec3 scale=Vec3(1.), Vec3 offset=Vec3(0.), Vec3i size=Vec3i(-1,-1,-1) , int orderSpace=1 )
{
	Vec3 sourceFactor(1.), off2 = offset;
	calcGridSizeFactorMod(source.getSize(), target.getSize(), size, scale, sourceFactor, off2);

	// a brief note on a mantaflow specialty: the target grid has to be the first argument here!
	// the parent fluidsolver object is taken from the first grid, and it determines the size of the
	// loop for the kernel call. as we're writing into target, it's important to loop exactly over
	// all cells of the target grid... (note, when calling the plugin in python, it doesnt matter anymore).

	// sourceFactor offset necessary to shift eval points by half a small cell width
	knInterpolateGridTempl<Real>(target, source, sourceFactor, off2, orderSpace);
}

PYTHON() void interpolateGridVec3( Grid<Vec3>& target, Grid<Vec3>& source , Vec3 scale=Vec3(1.), Vec3 offset=Vec3(0.), Vec3i size=Vec3i(-1,-1,-1) , int orderSpace=1 )
{
	Vec3 sourceFactor(1.), off2 = offset;
	calcGridSizeFactorMod(source.getSize(), target.getSize(), size, scale, sourceFactor, off2);
	knInterpolateGridTempl<Vec3>(target, source, sourceFactor, off2, orderSpace);
}


//!interpolate a mac velocity grid from one size to another size
KERNEL() 
void KnInterpolateMACGrid(MACGrid& target, MACGrid& source, const Vec3& sourceFactor, const Vec3& off, int orderSpace)
{
	Vec3 pos = Vec3(i,j,k) * sourceFactor + off;

	Real vx = source.getInterpolatedHi(pos - Vec3(0.5,0,0), orderSpace)[0];
	Real vy = source.getInterpolatedHi(pos - Vec3(0,0.5,0), orderSpace)[1];
	Real vz = 0.f;
	if(source.is3D()) vz = source.getInterpolatedHi(pos - Vec3(0,0,0.5), orderSpace)[2];

	target(i,j,k) = Vec3(vx,vy,vz);
}

PYTHON() void interpolateMACGrid(MACGrid& target, MACGrid& source, Vec3 scale=Vec3(1.), Vec3 offset=Vec3(0.), Vec3i size=Vec3i(-1,-1,-1) , int orderSpace=1)
{
	Vec3 sourceFactor(1.), off2 = offset;
	calcGridSizeFactorMod(source.getSize(), target.getSize(), size, scale, sourceFactor, off2);
	KnInterpolateMACGrid(target, source, sourceFactor, off2, orderSpace);
}



//*****************************************************************************

//! Apply vector noise to grid, this is a simplified version - no position scaling or UVs
KERNEL() 
void knApplySimpleNoiseVec3(FlagGrid& flags, Grid<Vec3>& target, WaveletNoiseField& noise, 
					  Real scale, Grid<Real>* weight ) 
{
	if ( !flags.isFluid(i,j,k) ) return; 
	Real factor = 1;
	if(weight) factor = (*weight)(i,j,k);
	target(i,j,k) += noise.evaluateCurl( Vec3(i,j,k)+Vec3(0.5) ) * scale * factor;
}
PYTHON() void applySimpleNoiseVec3(FlagGrid& flags, Grid<Vec3>& target, WaveletNoiseField& noise, 
							Real scale=1.0 , Grid<Real>* weight=NULL )
{
	// note - passing a MAC grid here is slightly inaccurate, we should evaluate each component separately
	knApplySimpleNoiseVec3(flags, target, noise, scale , weight );
}


//! Simple noise for a real grid , follows applySimpleNoiseVec3
KERNEL() 
void knApplySimpleNoiseReal(FlagGrid& flags, Grid<Real>& target, WaveletNoiseField& noise, 
					  Real scale, Grid<Real>* weight ) 
{
	if ( !flags.isFluid(i,j,k) ) return; 
	Real factor = 1;
	if(weight) factor = (*weight)(i,j,k);
	target(i,j,k) += noise.evaluate( Vec3(i,j,k)+Vec3(0.5) ) * scale * factor;
}
PYTHON() void applySimpleNoiseReal(FlagGrid& flags, Grid<Real>& target, WaveletNoiseField& noise, 
							Real scale=1.0 , Grid<Real>* weight=NULL )
{
	knApplySimpleNoiseReal(flags, target, noise, scale , weight );
}



//! Apply vector-based wavelet noise to target grid
//! This is the version with more functionality - supports uv grids, and on-the-fly interpolation
//! of input grids.
KERNEL() 
void knApplyNoiseVec3(FlagGrid& flags, Grid<Vec3>& target, WaveletNoiseField& noise, 
					  Real scale, Real scaleSpatial, Grid<Real>* weight, Grid<Vec3>* uv, bool uvInterpol, const Vec3& sourceFactor ) 
{
	if ( !flags.isFluid(i,j,k) ) return;

	// get weighting, interpolate if necessary
	Real w = 1;
	if(weight) {
		if(!uvInterpol) {
			w = (*weight)(i,j,k);
		} else {
			w = weight->getInterpolated( Vec3(i,j,k) * sourceFactor );
		}
	}

	// compute position where to evaluate the noise
	Vec3 pos = Vec3(i,j,k)+Vec3(0.5);
	if(uv) {
		if(!uvInterpol) {
			pos = (*uv)(i,j,k);
		} else {
			pos = uv->getInterpolated( Vec3(i,j,k) * sourceFactor );
			// uv coordinates are in local space - so we need to adjust the values of the positions
			pos /= sourceFactor;
		}
	}
	pos *= scaleSpatial;

	Vec3 noiseVec3 = noise.evaluateCurl( pos ) * scale * w; 
	//noiseVec3=pos; // debug , show interpolated positions
	target(i,j,k) += noiseVec3;
} 
PYTHON() void applyNoiseVec3(FlagGrid& flags, Grid<Vec3>& target, WaveletNoiseField& noise, 
							Real scale=1.0 , Real scaleSpatial=1.0 , Grid<Real>* weight=NULL , Grid<Vec3>* uv=NULL )
{
	// check whether the uv grid has a different resolution
	bool uvInterpol = false; 
	// and pre-compute conversion (only used if uvInterpol==true)
	// used for both uv and weight grid...
	Vec3 sourceFactor = Vec3(1.);
	if(uv) {
		uvInterpol = (target.getSize() != uv->getSize());
		sourceFactor = calcGridSizeFactor( uv->getSize(), target.getSize() );
	} else if(weight) {
		uvInterpol = (target.getSize() != weight->getSize());
		sourceFactor = calcGridSizeFactor( weight->getSize(), target.getSize() );
	}
	if(uv && weight) assertMsg( uv->getSize() == weight->getSize(), "UV and weight grid have to match!");

	// note - passing a MAC grid here is slightly inaccurate, we should evaluate each component separately
	knApplyNoiseVec3(flags, target, noise, scale, scaleSpatial, weight , uv,uvInterpol,sourceFactor );
}



//! Compute energy of a staggered velocity field (at cell center)
KERNEL() 
void KnApplyComputeEnergy( FlagGrid& flags, MACGrid& vel, Grid<Real>& energy ) 
{
	Real e = 0.f;
	if ( flags.isFluid(i,j,k) ) {
		Vec3 v = vel.getCentered(i,j,k);
		e = 0.5 * (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
	}
	energy(i,j,k) = e;
}

PYTHON() void computeEnergy( FlagGrid& flags, MACGrid& vel, Grid<Real>& energy )
{
	KnApplyComputeEnergy( flags, vel, energy );
}


PYTHON() void computeWaveletCoeffs(Grid<Real>& input)
{
	Grid<Real> temp1(input.getParent()), temp2(input.getParent());
	WaveletNoiseField::computeCoefficients(input, temp1, temp2);
}

// note - alomst the same as for vorticity confinement
PYTHON() void computeVorticity(MACGrid& vel, Grid<Vec3>& vorticity, Grid<Real>* norm) {
	Grid<Vec3> velCenter(vel.getParent());
	GetCentered(velCenter, vel);
	CurlOp(velCenter, vorticity);
	if(norm) GridNorm( *norm, vorticity);
}

// note - very similar to KnComputeProductionStrain, but for use as wavelet turb weighting
KERNEL(bnd=1) 
void KnComputeStrainRateMag(const MACGrid& vel, const Grid<Vec3>& velCenter, Grid<Real>& prod ) 
{
	// compute Sij = 1/2 * (dU_i/dx_j + dU_j/dx_i)
	Vec3 diag = Vec3(vel(i+1,j,k).x, vel(i,j+1,k).y, 0. ) - vel(i,j,k);
	if(vel.is3D()) diag[2] += vel(i,j,k+1).z;
	else           diag[2]  = 0.;

	Vec3 ux =         0.5*(velCenter(i+1,j,k)-velCenter(i-1,j,k));
	Vec3 uy =         0.5*(velCenter(i,j+1,k)-velCenter(i,j-1,k));
	Vec3 uz;
	if(vel.is3D()) uz=0.5*(velCenter(i,j,k+1)-velCenter(i,j,k-1));

	Real S12 = 0.5*(ux.y+uy.x);
	Real S13 = 0.5*(ux.z+uz.x);
	Real S23 = 0.5*(uy.z+uz.y);
	Real S2 = square(diag.x) + square(diag.y) + square(diag.z) +
		2.0*square(S12) + 2.0*square(S13) + 2.0*square(S23);
	prod(i,j,k) = S2;
}
PYTHON() void computeStrainRateMag(MACGrid& vel, Grid<Real>& mag) {
	Grid<Vec3> velCenter(vel.getParent());
	GetCentered(velCenter, vel);
	KnComputeStrainRateMag(vel, velCenter, mag);
}


// extrapolate a real grid into a flagged region (based on initial flags)
// by default extrapolates from fluid to obstacle cells
template<class T> 
void extrapolSimpleFlagsHelper (FlagGrid& flags, Grid<T>& val, int distance = 4, 
									int flagFrom=FlagGrid::TypeFluid, int flagTo=FlagGrid::TypeObstacle ) 
{
	Grid<int> tmp( flags.getParent() );
	int dim = (flags.is3D() ? 3:2);
	const Vec3i nb[6] = { 
		Vec3i(1 ,0,0), Vec3i(-1,0,0),
		Vec3i(0,1 ,0), Vec3i(0,-1,0),
		Vec3i(0,0,1 ), Vec3i(0,0,-1) };

	// remove all fluid cells (set to 1)
	tmp.clear();
	bool foundTarget = false;
	FOR_IJK_BND(flags,0) {
		if (flags(i,j,k) & flagFrom) 
			tmp( Vec3i(i,j,k) ) = 1;
		if (!foundTarget && (flags(i,j,k) & flagTo)) foundTarget=true;
	}
	// optimization, skip extrapolation if we dont have any cells to extrapolate to
	if(!foundTarget) {
		debMsg("No target cells found, skipping extrapolation", 1);
		return;
	}

	// extrapolate for given distance
	for(int d=1; d<1+distance; ++d) {

		// TODO, parallelize
		FOR_IJK_BND(flags,1) {
			if (tmp(i,j,k) != 0)          continue;
			if (!(flags(i,j,k) & flagTo)) continue;

			// copy from initialized neighbors
			Vec3i p(i,j,k);
			int nbs = 0;
			T avgVal = 0.;
			for (int n=0; n<2*dim; ++n) {
				if (tmp(p+nb[n]) == d) {
					avgVal += val(p+nb[n]);
					nbs++;
				}
			}

			if(nbs>0) {
				tmp(p) = d+1;
				val(p) = avgVal / nbs;
			}
		}

	} // distance 
}
PYTHON() void extrapolateSimpleFlags (FlagGrid& flags, GridBase* val, int distance = 4, 
									int flagFrom=FlagGrid::TypeFluid, int flagTo=FlagGrid::TypeObstacle ) 
{
	if (val->getType() & GridBase::TypeReal) {
		extrapolSimpleFlagsHelper<Real>(flags,*((Grid<Real>*) val),distance,flagFrom,flagTo);
	}
	else if (val->getType() & GridBase::TypeInt) {    
		extrapolSimpleFlagsHelper<int >(flags,*((Grid<int >*) val),distance,flagFrom,flagTo);
	}
	else if (val->getType() & GridBase::TypeVec3) {    
		extrapolSimpleFlagsHelper<Vec3>(flags,*((Grid<Vec3>*) val),distance,flagFrom,flagTo);
	}
	else
		errMsg("extrapolateSimpleFlags: Grid Type is not supported (only int, Real, Vec3)");    
}

//! convert vel to a centered grid, then compute its curl
PYTHON() void getCurl(MACGrid& vel, Grid<Real>& vort, int comp) {
	Grid<Vec3> velCenter(vel.getParent()), curl(vel.getParent());
	
	GetCentered(velCenter, vel);
	CurlOp(velCenter, curl);
	GetComponent(curl, vort, comp);
}
PYTHON() void setObsFromLevelset(FlagGrid& target, LevelsetGrid& levelset, Real threshold = 1.0, int flagv = 2) {
	FOR_IDX(target) {
		//if (target.isObstacle(idx)) continue;
		const Real phi = levelset[idx];
		if (phi <= threshold) {
			//target(idx) &= ~(FlagGrid::TypeEmpty | FlagGrid::TypeFluid); // clear empty/fluid flags
			target(idx) = flagv; // set resepctive flag
		}
	}
}
PYTHON() void clampRealGridFromLevelset(Grid<Real>& target, LevelsetGrid& levelset, Real threshold = 1.5, Real v = 0.0) {

	FOR_IDX(target) {

		//if (target.isObstacle(idx)) continue;

		const Real phi = levelset[idx];

		if (phi <= threshold) {

			//target(idx) &= ~(FlagGrid::TypeEmpty | FlagGrid::TypeFluid); // clear empty/fluid flags

			target(idx) = v; // set resepctive flag

		}

	}

}

} // namespace
