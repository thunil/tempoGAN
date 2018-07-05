/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Conjugate gradient solver
 *
 ******************************************************************************/

#ifndef _CONJUGATEGRADIENT_H
#define _CONJUGATEGRADIENT_H

#include "vectorbase.h"
#include "grid.h"
#include "kernel.h"
#include "multigrid.h"

namespace Manta { 

static const bool CG_DEBUG = false;

//! Basic CG interface 
class GridCgInterface {
	public:
		enum PreconditionType { PC_None=0, PC_ICP, PC_mICP, PC_MGP };
		
		GridCgInterface() : mUseL2Norm(true) {};
		virtual ~GridCgInterface() {};

		// solving functions
		virtual bool iterate() = 0;
		virtual void solve(int maxIter) = 0;

		// precond
		virtual void setICPreconditioner(PreconditionType method, Grid<Real> *A0, Grid<Real> *Ai, Grid<Real> *Aj, Grid<Real> *Ak) = 0;
		virtual void setMGPreconditioner(PreconditionType method, GridMg* MG) = 0;

		// access
		virtual Real getSigma() const = 0;
		virtual Real getIterations() const = 0;
		virtual Real getResNorm() const = 0;
		virtual void setAccuracy(Real set) = 0;
		virtual Real getAccuracy() const = 0;

		//! force reinit upon next iterate() call, can be used for doing multiple solves
		virtual void forceReinit() = 0;

		void setUseL2Norm(bool set) { mUseL2Norm = set; }

	protected:

		// use l2 norm of residualfor threshold? (otherwise uses max norm)
		bool mUseL2Norm; 
};


//! Run single iteration of the cg solver
/*! the template argument determines the type of matrix multiplication,
	typically a ApplyMatrix kernel, another one is needed e.g. for the
	mesh-based wave equation solver */
template<class APPLYMAT>
class GridCg : public GridCgInterface {
	public:
		//! constructor
		GridCg(Grid<Real>& dst, Grid<Real>& rhs, Grid<Real>& residual, Grid<Real>& search, const FlagGrid& flags, Grid<Real>& tmp, 
				Grid<Real>* A0, Grid<Real>* pAi, Grid<Real>* pAj, Grid<Real>* pAk);
		~GridCg() {}
		
		void doInit();
		bool iterate();
		void solve(int maxIter);
		//! init pointers, and copy values from "normal" matrix
		void setICPreconditioner(PreconditionType method, Grid<Real> *A0, Grid<Real> *Ai, Grid<Real> *Aj, Grid<Real> *Ak);
		void setMGPreconditioner(PreconditionType method, GridMg* MG);
		void forceReinit() { mInited = false; }
		
		// Accessors        
		Real getSigma() const { return mSigma; }
		Real getIterations() const { return mIterations; }

		Real getResNorm() const { return mResNorm; }

		void setAccuracy(Real set) { mAccuracy=set; }
		Real getAccuracy() const { return mAccuracy; }

	protected:
		bool mInited;
		int mIterations;
		// grids
		Grid<Real>& mDst;
		Grid<Real>& mRhs;
		Grid<Real>& mResidual;
		Grid<Real>& mSearch;
		const FlagGrid& mFlags;
		Grid<Real>& mTmp;

		Grid<Real> *mpA0, *mpAi, *mpAj, *mpAk;

		PreconditionType mPcMethod;
		//! preconditioning grids
		Grid<Real> *mpPCA0, *mpPCAi, *mpPCAj, *mpPCAk;
		GridMg* mMG;

		//! sigma / residual
		Real mSigma;
		//! accuracy of solver (max. residuum)
		Real mAccuracy;
		//! norm of the residual
		Real mResNorm;
}; // GridCg


//! Kernel: Apply symmetric stored Matrix
KERNEL(idx) 
void ApplyMatrix (const FlagGrid& flags, Grid<Real>& dst, const Grid<Real>& src, 
				  Grid<Real>& A0, Grid<Real>& Ai, Grid<Real>& Aj, Grid<Real>& Ak)
{
	if (!flags.isFluid(idx)) {
		dst[idx] = src[idx]; return;
	}    

	dst[idx] =  src[idx] * A0[idx]
				+ src[idx-X] * Ai[idx-X]
				+ src[idx+X] * Ai[idx]
				+ src[idx-Y] * Aj[idx-Y]
				+ src[idx+Y] * Aj[idx]
				+ src[idx-Z] * Ak[idx-Z] 
				+ src[idx+Z] * Ak[idx];
}

//! Kernel: Apply symmetric stored Matrix. 2D version
KERNEL(idx) 
void ApplyMatrix2D (const FlagGrid& flags, Grid<Real>& dst, const Grid<Real>& src, 
					Grid<Real>& A0, Grid<Real>& Ai, Grid<Real>& Aj, Grid<Real>& Ak)
{
	unusedParameter(Ak); // only there for parameter compatibility with ApplyMatrix
	
	if (!flags.isFluid(idx)) {
		dst[idx] = src[idx]; return;
	}    

	dst[idx] =  src[idx] * A0[idx]
				+ src[idx-X] * Ai[idx-X]
				+ src[idx+X] * Ai[idx]
				+ src[idx-Y] * Aj[idx-Y]
				+ src[idx+Y] * Aj[idx];
}

//! Kernel: Construct the matrix for the poisson equation
KERNEL (bnd=1) 
void MakeLaplaceMatrix(const FlagGrid& flags, Grid<Real>& A0, Grid<Real>& Ai, Grid<Real>& Aj, Grid<Real>& Ak, const MACGrid* fractions = 0) {
	if (!flags.isFluid(i,j,k))
		return;
	
	if(!fractions) {
		// diagonal, A0
		if (!flags.isObstacle(i-1,j,k))                 A0(i,j,k) += 1.;
		if (!flags.isObstacle(i+1,j,k))                 A0(i,j,k) += 1.;
		if (!flags.isObstacle(i,j-1,k))                 A0(i,j,k) += 1.;
		if (!flags.isObstacle(i,j+1,k))                 A0(i,j,k) += 1.;
		if (flags.is3D() && !flags.isObstacle(i,j,k-1)) A0(i,j,k) += 1.;
		if (flags.is3D() && !flags.isObstacle(i,j,k+1)) A0(i,j,k) += 1.;
		
		// off-diagonal entries
		if (flags.isFluid(i+1,j,k))                 Ai(i,j,k) = -1.;
		if (flags.isFluid(i,j+1,k))                 Aj(i,j,k) = -1.;
		if (flags.is3D() && flags.isFluid(i,j,k+1)) Ak(i,j,k) = -1.;
	} else {
		// diagonal
		A0(i,j,k)                   += fractions->get(i  ,j,k).x;
		A0(i,j,k)                   += fractions->get(i+1,j,k).x;
		A0(i,j,k)                   += fractions->get(i,j  ,k).y;
		A0(i,j,k)                   += fractions->get(i,j+1,k).y;
		if (flags.is3D()) A0(i,j,k) += fractions->get(i,j,k  ).z;
		if (flags.is3D()) A0(i,j,k) += fractions->get(i,j,k+1).z;

		// off-diagonal entries
		if (flags.isFluid(i+1,j,k))                 Ai(i,j,k) = -fractions->get(i+1,j,k).x;
		if (flags.isFluid(i,j+1,k))                 Aj(i,j,k) = -fractions->get(i,j+1,k).y;
		if (flags.is3D() && flags.isFluid(i,j,k+1)) Ak(i,j,k) = -fractions->get(i,j,k+1).z;
	}

}




} // namespace

#endif 
