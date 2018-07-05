/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Plugins for pressure correction: solve_pressure, and ghost fluid helpers
 *
 ******************************************************************************/
#include "vectorbase.h"
#include "grid.h"
#include "kernel.h"
#include "conjugategrad.h"
#include "rcmatrix.h"

using namespace std;
namespace Manta {

// only supports a single blur size for now, globals stored here
bool gBlurPrecomputed = false;
int  gBlurKernelRadius = -1;
Matrix gBlurKernel;

// *****************************************************************************
// Helper functions for fluid guiding

//! creates a 1D (horizontal) Gaussian blur kernel of size n and standard deviation sigma 
Matrix get1DGaussianBlurKernel(const int n, const int sigma) {
	Matrix x(n), y(n);
	for (int j = 0; j < n; j++) {
		x.add_to_element(0, j, - (n - 1)*0.5);
		y.add_to_element(0, j, j - (n - 1)*0.5);
	}
	Matrix G(n);
	Real sumG = 0;
	for (int j = 0; j < n; j++) {
		G.add_to_element(0, j, 1 / (2 * M_PI*sigma*sigma)*exp(-(x(0, j)*x(0, j) + y(0, j)*y(0, j)) / (2 * sigma*sigma)));
		sumG += G(0, j);
	}
	G = G*(1.0 / sumG);
	return G;
}

//! convolves in with 1D kernel (centred at the kernel's midpoint) in the x-direction
//! (out must be a grid of zeros) 
KERNEL() void apply1DKernelDirX(const MACGrid &in, MACGrid &out, const Matrix &kernel) {
	int nx = in.getSizeX();
	int kn = kernel.n;
	int kCentre = kn / 2;
	for (int m = 0, ind = kn - 1, ii = i - kCentre; m < kn; m++, ind--, ii++) {
		if (ii < 0) continue;
		else if (ii >= nx) break;
		else out(i, j, k) += in(ii, j, k)*kernel(0,ind);
	}
}

//! convolves in with 1D kernel (centred at the kernel's midpoint) in the y-direction
//! (out must be a grid of zeros)
KERNEL() void apply1DKernelDirY(const MACGrid &in, MACGrid &out, const Matrix &kernel) {
	int ny = in.getSizeY();
	int kn = kernel.n;
	int kCentre = kn / 2;
	for (int m = 0, ind = kn - 1, jj = j - kCentre; m < kn; m++, ind--, jj++) {
		if (jj < 0) continue;
		else if (jj >= ny) break;
		else out(i, j, k) += in(i, jj, k)*kernel(0, ind);
	}
}

//! convolves in with 1D kernel (centred at the kernel's midpoint) in the z-direction
//! (out must be a grid of zeros)
KERNEL() void apply1DKernelDirZ(const MACGrid &in, MACGrid &out, const Matrix &kernel) {
	int nz = in.getSizeZ();
	int kn = kernel.n;
	int kCentre = kn / 2;
	for (int m = 0, ind = kn - 1, kk = k - kCentre; m < kn; m++, ind--, kk++) {
		if (kk < 0) continue;
		else if (kk >= nz) break;
		else out(i, j, k) += in(i, j, kk)*kernel(0, ind);
	}
}

//! Apply separable Gaussian blur in 2D
void applySeparableKernel2D(MACGrid &grid, const FlagGrid &flags, const Matrix &kernel) {
	//int nx = grid.getSizeX(), ny = grid.getSizeY();
	//int kn = kernel.n;
	//int kCentre = kn / 2;
	FluidSolver* parent = grid.getParent();
	MACGrid orig = MACGrid(parent);
	orig.copyFrom(grid);
	MACGrid gridX = MACGrid(parent);
	apply1DKernelDirX(grid, gridX, kernel);
	MACGrid gridXY = MACGrid(parent);
	apply1DKernelDirY(gridX, gridXY, kernel);
	grid.copyFrom(gridXY);
	FOR_IJK(grid) {
		if ((i>0 && flags.isObstacle(i - 1, j, k)) || (j>0 && flags.isObstacle(i, j - 1, k)) || flags.isObstacle(i, j, k)) {
			grid(i, j, k).x = orig(i, j, k).x;
			grid(i, j, k).y = orig(i, j, k).y;
			grid(i, j, k).z = orig(i, j, k).z;
		}
	}
}

//! Apply separable Gaussian blur in 3D
void applySeparableKernel3D(MACGrid &grid, const FlagGrid &flags, const Matrix &kernel) {
	//int nx = grid.getSizeX(), ny = grid.getSizeY(), nz = grid.getSizeZ();
	//int kn = kernel.n;
	//int kCentre = kn / 2;
	FluidSolver* parent = grid.getParent();
	MACGrid orig = MACGrid(parent);
	orig.copyFrom(grid);
	MACGrid gridX = MACGrid(parent);
	apply1DKernelDirX(grid, gridX, kernel);
	MACGrid gridXY = MACGrid(parent);
	apply1DKernelDirY(gridX, gridXY, kernel);
	MACGrid gridXYZ = MACGrid(parent);
	apply1DKernelDirZ(gridXY, gridXYZ, kernel);
	grid.copyFrom(gridXYZ);
	FOR_IJK(grid) {
		if ((i>0 && flags.isObstacle(i - 1, j, k)) || (j>0 && flags.isObstacle(i, j - 1, k)) || (k>0 && flags.isObstacle(i, j, k - 1)) || flags.isObstacle(i, j, k)) {
			grid(i, j, k).x = orig(i, j, k).x;
			grid(i, j, k).y = orig(i, j, k).y;
			grid(i, j, k).z = orig(i, j, k).z;
		}
	}
}

//! Apply separable Gaussian blur in 2D or 3D depending on input dimensions
void applySeparableKernel(MACGrid &grid, const FlagGrid &flags, const Matrix &kernel) {
	if (!grid.is3D()) applySeparableKernel2D(grid, flags, kernel);
	else applySeparableKernel3D(grid, flags, kernel);
}


//! Compute r-norm for the stopping criterion
Real getRNorm(const MACGrid &x, const MACGrid &z) {
	MACGrid r = MACGrid(x.getParent());
	r.copyFrom(x);
	r.sub(z);
	return r.getMaxAbs();
}

//! Compute s-norm for the stopping criterion
Real getSNorm(const Real rho, const MACGrid &z, const MACGrid &z_prev) {
	MACGrid s = MACGrid(z_prev.getParent());
	s.copyFrom(z_prev);
	s.sub(z);
	s.multConst(rho);
	return s.getMaxAbs();
}

//! Compute primal eps for the stopping criterion
Real getEpsPri(const Real eps_abs, const Real eps_rel,
        const MACGrid &x, const MACGrid &z) {
	Real max_norm = max(x.getMaxAbs(), z.getMaxAbs());
	Real eps_pri = sqrt(x.is3D() ? 3.0 : 2.0)*eps_abs + eps_rel*max_norm;
	return eps_pri;
}

//! Compute dual eps for the stopping criterion
Real getEpsDual(const Real eps_abs, const Real eps_rel, const MACGrid &y) {
	Real eps_dual = sqrt(y.is3D() ? 3.0 : 2.0)*eps_abs + eps_rel*y.getMaxAbs();
	return eps_dual;
}

//! Create a spiral velocity field in 2D as a test scene
PYTHON() void getSpiralVelocity2D(const FlagGrid &flags, MACGrid &vel, Real strength = 1.0) {
	int nx = flags.getSizeX(), ny = flags.getSizeY();
	Real midX = 0.5*(Real)(nx - 1);
	Real midY = 0.5*(Real)(ny - 1);
	int k = 0;
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			int idx = flags.index(i, j, k);
			Real diffX = midX - i;
			Real diffY = midY - j;
			Real hypotenuse = sqrt(diffX*diffX + diffY*diffY);
			if (hypotenuse > 0) {
				vel[idx].x = diffY / hypotenuse;
				vel[idx].y = -diffX / hypotenuse;
			}
		}
	}
	vel.multConst(strength);
}

//! Set the guiding weight W as a gradient in the y-direction
PYTHON() void setGradientYWeight(Grid<Real> &W, const int minY, const int maxY, const Real valAtMin, const Real valAtMax) {
	FOR_IJK(W) {
		if (minY <= j && j <= maxY) {
			Real val = valAtMin;
			if (valAtMax != valAtMin) {
				Real ratio = (Real)(j - minY) / (Real)(maxY - minY);
				val = ratio*valAtMax + (1.0 - ratio)*valAtMin;
			}
			W(i, j, k) = val;
		}
	}
}

// *****************************************************************************
// More helper functions for fluid guiding
	
//! Apply Gaussian blur (either 2D or 3D) in a separable way
void applySeparableGaussianBlur(MACGrid &grid, const FlagGrid &flags, const Matrix &kernel1D) {
	assertMsg(gBlurPrecomputed, "Error - blue kernel not precomputed");
	applySeparableKernel(grid, flags, kernel1D);
}

//! Precomputation performed before the first PD iteration
void ADMM_precompute_Separable(int blurRadius) {
	if (gBlurPrecomputed) {
		assertMsg( gBlurKernelRadius == blurRadius, "More than a single blur radius not supported at the moment." );
		return;
	}
	int kernelSize = 2 * blurRadius + 1;
	gBlurKernel = get1DGaussianBlurKernel(kernelSize, kernelSize);
	gBlurPrecomputed = true;
	gBlurKernelRadius = blurRadius;
}

//! Apply approximate multiplication of inverse(M)
void applyApproxInvM(MACGrid& v, const FlagGrid &flags, const MACGrid& invA) {
	MACGrid v_new = MACGrid(v.getParent());
	v_new.copyFrom(v); 
	v_new.mult(invA);
	applySeparableGaussianBlur(v_new, flags, gBlurKernel); 
	applySeparableGaussianBlur(v_new, flags, gBlurKernel);
	v_new.multConst(2.0); 
	v_new.mult(invA);
	v.mult(invA);
	v.sub(v_new);
}

//! Precompute Q, a reused quantity in the PD iterations
//! Q = 2*G*G*(velT-velC)-sigma*velC
void precomputeQ(MACGrid &Q, const FlagGrid &flags, const MACGrid &velT_region, const MACGrid &velC, const Matrix &gBlurKernel, const Real sigma) {
	Q.copyFrom(velT_region);
	Q.sub(velC);
	applySeparableGaussianBlur(Q, flags, gBlurKernel); 
	applySeparableGaussianBlur(Q, flags, gBlurKernel);
	Q.multConst(2.0);
	Q.addScaled(velC, -sigma);
}

//! Precompute inverse(A), a reused quantity in the PD iterations
//! A = 2*S^2 + p*I, invA = elementwise 1/A
void precomputeInvA(MACGrid &invA, const Grid<Real> &weight, const Real sigma) {
	FOR_IJK(invA) {
		Real val = 2 * weight(i, j, k)*weight(i, j, k) + sigma;
		if (val<0.01) val = 0.01;
		Real invVal = 1.0 / val;
		invA(i, j, k).x = invVal;
		invA(i, j, k).y = invVal;
		invA(i, j, k).z = invVal;
	}
}

//! proximal operator of f , guiding 
void prox_f(MACGrid& v, const FlagGrid &flags, const MACGrid& Q, const MACGrid& velC, const Real sigma, const MACGrid& invA) {
	v.multConst(sigma);
	v.add(Q);
	applyApproxInvM(v, flags, invA);
	v.add(velC);
}

// *****************************************************************************

// re-uses main pressure solve from pressure.cpp
void solvePressure( 
        MACGrid& vel, Grid<Real>& pressure, const FlagGrid& flags, Real cgAccuracy = 1e-3,
    const Grid<Real>* phi = 0,
    const Grid<Real>* perCellCorr = 0,
    const MACGrid* fractions = 0,
    Real gfClamp = 1e-04,
    Real cgMaxIterFac = 1.5,
    bool precondition = true, 
	int preconditioner = 1,
	bool enforceCompatibility = false,
    bool useL2Norm = false, 
	bool zeroPressureFixing = false,
	const Grid<Real> *curv = NULL,
	const Real surfTens = 0.0,
	Grid<Real>* retRhs = NULL );

//! Main function for fluid guiding , includes "regular" pressure solve
PYTHON() void PD_fluid_guiding(MACGrid& vel, MACGrid& velT,
	Grid<Real>& pressure, FlagGrid& flags, Grid<Real>& weight, int blurRadius = 5,
	Real theta = 1.0, Real tau = 1.0, Real sigma = 1.0,
	Real epsRel = 1e-3, Real epsAbs = 1e-3, int maxIters = 200,
	// duplicated for pressure solve
	Grid<Real>* phi = 0, Grid<Real>* perCellCorr = 0, MACGrid* fractions = 0, Real gfClamp = 1e-04, Real cgMaxIterFac = 1.5, Real cgAccuracy = 1e-3,
	int preconditioner = 1, bool zeroPressureFixing = false)
{
	FluidSolver* parent = vel.getParent();

	// initialize dual/slack variables
	MACGrid velC = MACGrid(parent); velC.copyFrom(vel);
	MACGrid x = MACGrid(parent);
	MACGrid y = MACGrid(parent);
	MACGrid z = MACGrid(parent);
	MACGrid x0 = MACGrid(parent);
	MACGrid z0 = MACGrid(parent);

	// precomputation
	ADMM_precompute_Separable(blurRadius);
	MACGrid Q = MACGrid(parent);
	precomputeQ(Q, flags, velT, velC, gBlurKernel, sigma);
	MACGrid invA = MACGrid(parent);
	precomputeInvA(invA, weight, sigma);

	// loop
	int iter = 0;
	for (iter = 0; iter < maxIters; iter++) {
		// x-update
		x0.copyFrom(x);
		x.multConst(1.0 / sigma); 
		x.add(y); 
		prox_f(x, flags, Q, velC, sigma, invA);
		x.multConst(-sigma); x.addScaled(y, sigma); x.add(x0);

		// z-update
		z0.copyFrom(z);
		z.addScaled(x, -tau);
		Real cgAccuracyAdaptive = cgAccuracy;

		solvePressure (z, pressure, flags, cgAccuracyAdaptive, phi, perCellCorr, fractions, gfClamp,
		    cgMaxIterFac, true, preconditioner, false, false, zeroPressureFixing );

		// y-update
		y.copyFrom(z);
		y.sub(z0);
		y.multConst(theta);
		y.add(z);

		// stopping criterion
		bool stop = (iter > 0 && getRNorm(z, z0) < getEpsDual(epsAbs, epsRel, z));

		if (stop || (iter == maxIters - 1)) break;
	}

	// vel_new = z
	vel.copyFrom(z);

	debMsg("PD_fluid_guiding iterations:" << iter, 1);
}


} // end namespace
