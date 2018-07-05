/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Turbulence modeling plugins
 *
 ******************************************************************************/
 
#include "grid.h"
#include "commonkernels.h"
#include "vortexsheet.h"
#include "conjugategrad.h"

using namespace std;

namespace Manta {

// k-epsilon model constants
const Real keCmu = 0.09;
const Real keC1 = 1.44;
const Real keC2 = 1.92;
const Real keS1 = 1.0;
const Real keS2 = 1.3;

// k-epsilon limiters
const Real keU0 = 1.0;
const Real keImin = 2e-3;
const Real keImax = 1.0;
const Real keNuMin = 1e-3;
const Real keNuMax = 5.0;

//! clamp k and epsilon to limits    
KERNEL(idx) 
void KnTurbulenceClamp(Grid<Real>& kgrid, Grid<Real>& egrid, Real minK, Real maxK, Real minNu, Real maxNu) {
	Real eps = egrid[idx];
	Real ke = clamp(kgrid[idx],minK,maxK);
	Real nu = keCmu*square(ke)/eps;
	if (nu > maxNu) 
		eps = keCmu*square(ke)/maxNu;
	if (nu < minNu) 
		eps = keCmu*square(ke)/minNu;

	kgrid[idx] = ke;
	egrid[idx] = eps;
}

//! Compute k-epsilon production term P = 2*nu_T*sum_ij(Sij^2) and the turbulent viscosity nu_T=C_mu*k^2/eps
KERNEL(bnd=1) 
void KnComputeProduction(const MACGrid& vel, const Grid<Vec3>& velCenter, const Grid<Real>& ke, const Grid<Real>& eps, 
							   Grid<Real>& prod, Grid<Real>& nuT, Grid<Real>* strain, Real pscale = 1.0f) 
{
	Real curEps = eps(i,j,k);
	if (curEps > 0) {
		// turbulent viscosity: nu_T = C_mu * k^2/eps
		Real curNu = keCmu * square(ke(i,j,k)) / curEps;
		
		// compute Sij = 1/2 * (dU_i/dx_j + dU_j/dx_i)
		Vec3 diag = Vec3(vel(i+1,j,k).x, vel(i,j+1,k).y, vel(i,j,k+1).z) - vel(i,j,k);
		Vec3 ux = 0.5*(velCenter(i+1,j,k)-velCenter(i-1,j,k));
		Vec3 uy = 0.5*(velCenter(i,j+1,k)-velCenter(i,j-1,k));
		Vec3 uz = 0.5*(velCenter(i,j,k+1)-velCenter(i,j,k-1));
		Real S12 = 0.5*(ux.y+uy.x);
		Real S13 = 0.5*(ux.z+uz.x);
		Real S23 = 0.5*(uy.z+uz.y);
		Real S2 = square(diag.x) + square(diag.y) + square(diag.z) +
				  2.0*square(S12) + 2.0*square(S13) + 2.0*square(S23);
		
		// P = 2*nu_T*sum_ij(Sij^2)
		prod(i,j,k) = 2.0 * curNu * S2 * pscale;
		nuT(i,j,k) = curNu;
		if (strain) (*strain)(i,j,k) = sqrt(S2);
	} 
	else {
		prod(i,j,k) = 0;
		nuT(i,j,k) = 0;
		if (strain) (*strain)(i,j,k) = 0;
	}
}
	
//! Compute k-epsilon production term P = 2*nu_T*sum_ij(Sij^2) and the turbulent viscosity nu_T=C_mu*k^2/eps
PYTHON() void KEpsilonComputeProduction(const MACGrid& vel, Grid<Real>& k, Grid<Real>& eps, Grid<Real>& prod, Grid<Real>& nuT, Grid<Real>* strain=0, Real pscale = 1.0f)
{
	// get centered velocity grid
	Grid<Vec3> vcenter(k.getParent());
	GetCentered(vcenter, vel);
	FillInBoundary(vcenter,1);
	
	// compute limits
	const Real minK = 1.5*square(keU0)*square(keImin);
	const Real maxK = 1.5*square(keU0)*square(keImax);    
	KnTurbulenceClamp(k, eps, minK, maxK, keNuMin, keNuMax);
	
	KnComputeProduction(vel, vcenter, k, eps, prod, nuT, strain, pscale);    
}

//! Integrate source terms of k-epsilon equation
KERNEL(idx) 
void KnAddTurbulenceSource(Grid<Real>& kgrid, Grid<Real>& egrid, const Grid<Real>& pgrid, Real dt) {
	Real eps = egrid[idx], prod = pgrid[idx], ke = kgrid[idx];
	if (ke <= 0) ke = 1e-3; // pre-clamp to avoid nan
	
	Real newK = ke + dt*(prod - eps);
	Real newEps = eps + dt*(prod * keC1 - eps * keC2) * (eps / ke);
	if (newEps <= 0) newEps = 1e-4; // pre-clamp to avoid nan

	kgrid[idx] = newK;
	egrid[idx] = newEps;
}


//! Integrate source terms of k-epsilon equation
PYTHON() void KEpsilonSources(Grid<Real>& k, Grid<Real>& eps, Grid<Real>& prod) {
	Real dt = k.getParent()->getDt();
		
	KnAddTurbulenceSource(k, eps, prod, dt);
	
	// compute limits
	const Real minK = 1.5*square(keU0)*square(keImin);
	const Real maxK = 1.5*square(keU0)*square(keImax);
	KnTurbulenceClamp(k, eps, minK, maxK, keNuMin, keNuMax);    
}

//! Initialize the domain or boundary conditions
PYTHON() void KEpsilonBcs(const FlagGrid& flags, Grid<Real>& k, Grid<Real>& eps, Real intensity, Real nu, bool fillArea) {
	// compute limits
	const Real vk = 1.5*square(keU0)*square(intensity);
	const Real ve = keCmu*square(vk) / nu;
	
	FOR_IDX(k) {
		if (fillArea || flags.isObstacle(idx)) {
			k[idx] = vk;
			eps[idx] = ve;
		}
	}
}

//! Gradient diffusion smoothing. Not unconditionally stable -- should probably do substepping etc.
void ApplyGradDiff(const Grid<Real>& grid, Grid<Real>& res, const Grid<Real>& nu, Real dt, Real sigma) {
	// should do this (but requires better boundary handling)
	/*MACGrid grad(grid.getParent());
	GradientOpMAC(grad, grid);
	grad *= nu;
	DivergenceOpMAC(res, grad);
	res *= dt/sigma;  */
	
	LaplaceOp(res, grid);
	res *= nu;
	res *= dt/sigma;
}

//! Compute k-epsilon turbulent viscosity
PYTHON() void KEpsilonGradientDiffusion(Grid<Real>& k, Grid<Real>& eps, Grid<Real>& nuT, Real sigmaU=4.0, MACGrid* vel=0) {
	Real dt = k.getParent()->getDt();
	Grid<Real> res(k.getParent());
	
	// gradient diffusion of k
	ApplyGradDiff(k, res, nuT, dt, keS1);
	k += res;

	// gradient diffusion of epsilon
	ApplyGradDiff(eps, res, nuT, dt, keS2);
	eps += res;
	
	// gradient diffusion of velocity
	if (vel) {
		Grid<Real> vc(k.getParent());
		for (int c=0; c<3; c++) {
			GetComponent(*vel, vc, c);
			ApplyGradDiff(vc, res, nuT, dt, sigmaU);
			vc += res;
			SetComponent(*vel, vc, c);    
		}
	}
}



} // namespace
