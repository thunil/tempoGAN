/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Wave equation
 *
 ******************************************************************************/

#include "levelset.h"
#include "commonkernels.h"
#include "particle.h"
#include "conjugategrad.h"
#include <cmath>

using namespace std;

namespace Manta {


/******************************************************************************
 *
 * explicit integration
 *
 ******************************************************************************/


KERNEL(bnd=1) 
void knCalcSecDeriv2d(const Grid<Real>& v, Grid<Real>& ret) { 
    ret(i,j,k) = 
		( -4. * v(i,j,k) + v(i-1,j,k) + v(i+1,j,k) + v(i,j-1,k) + v(i,j+1,k) );
};

//! calculate a second derivative for the wave equation
PYTHON() void calcSecDeriv2d(const Grid<Real>& v, Grid<Real>& curv) {
	knCalcSecDeriv2d(v,curv);
}


// mass conservation 

KERNEL(bnd=1, reduce=+) returns(double sum=0)
double knTotalSum(Grid<Real>& h) { sum += h(i,j,k); }

//! calculate the sum of all values in a grid (for wave equation solves)
PYTHON() Real totalSum(Grid<Real>& height) {
	knTotalSum ts(height);
	return ts.sum;
}

//! normalize all values in a grid (for wave equation solves)
PYTHON() void normalizeSumTo(Grid<Real>& height, Real target) {
	knTotalSum ts(height);
	Real factor = target / ts.sum;
	height.multConst(factor);
}


/******************************************************************************
 *
 * implicit time integration
 *
 ******************************************************************************/



//! Kernel: Construct the right-hand side of the poisson equation
KERNEL(bnd=1)
void MakeRhsWE(const FlagGrid& flags, Grid<Real>& rhs, const Grid<Real>& ut, const Grid<Real>& utm1,
			Real s, bool crankNic=false) 
{
	rhs(i,j,k) = ( 2.*ut(i,j,k) - utm1(i,j,k) );
	if(crankNic) {
		rhs(i,j,k) += s * ( -4.*ut(i,j,k) + 1.*ut(i-1,j,k) + 1.*ut(i+1,j,k) + 1.*ut(i,j-1,k) + 1.*ut(i,j+1,k) );
	} 
}





//! do a CG solve for the wave equation (note, out grid only there for debugging... could be removed)
PYTHON() void cgSolveWE(const FlagGrid& flags, Grid<Real>& ut, Grid<Real>& utm1, Grid<Real>& out,
						bool crankNic     = false,
						Real cSqr         = 0.25,
						Real cgMaxIterFac = 1.5,
						Real cgAccuracy   = 1e-5 )
{
	// reserve temp grids
	FluidSolver* parent = flags.getParent();
	Grid<Real> rhs(parent);
	Grid<Real> residual(parent);
	Grid<Real> search(parent);
	Grid<Real> A0(parent);
	Grid<Real> Ai(parent);
	Grid<Real> Aj(parent);
	Grid<Real> Ak(parent);
	Grid<Real> tmp(parent);
	// solution...
	out.clear();
		
	// setup matrix and boundaries
	MakeLaplaceMatrix (flags, A0, Ai, Aj, Ak);
	Real dt   = parent->getDt();
	Real s    = dt*dt*cSqr * 0.5;
	FOR_IJK(flags) {
		Ai(i,j,k) *= s;
		Aj(i,j,k) *= s;
		Ak(i,j,k) *= s;
		A0(i,j,k) *= s;
		A0(i,j,k) += 1.;
	}
	
	// compute divergence and init right hand side
	rhs.clear();
	// h=dt
	// rhs:   = 2 ut - ut-1 
   	// A:    (h2 c2/ dx)=s   ,  (1+4s)uij + s ui-1j + ...
	// Cr.Nic.
	// rhs:  cr nic = 2 ut - ut-1 + h^2c^2/2 b 
   	// A:    (h2 c2/2 dx)=s   ,  (1+4s)uij + s ui-1j + ...
	MakeRhsWE kernMakeRhs(flags, rhs, ut,utm1, s, crankNic);
	
	const int maxIter = (int)(cgMaxIterFac * flags.getSize().max()) * (flags.is3D() ? 1 : 4);
	GridCgInterface *gcg;
	if (flags.is3D())
		gcg = new GridCg<ApplyMatrix  >(out, rhs, residual, search, flags, tmp, &A0, &Ai, &Aj, &Ak );
	else
		gcg = new GridCg<ApplyMatrix2D>(out, rhs, residual, search, flags, tmp, &A0, &Ai, &Aj, &Ak );
	
	gcg->setAccuracy( cgAccuracy ); 

	// no preconditioning for now...
	for (int iter=0; iter<maxIter; iter++) {
		if (!gcg->iterate()) iter=maxIter;
	} 
	debMsg("cgSolveWaveEq iterations:"<<gcg->getIterations()<<", res:"<<gcg->getSigma(), 1);

	utm1.swap( ut );
	ut.copyFrom( out );

	delete gcg;
}




} //namespace

