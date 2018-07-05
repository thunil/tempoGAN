/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Common grid kernels
 *
 ******************************************************************************/

#ifndef _COMMONKERNELS_H
#define _COMMONKERNELS_H

#include "general.h"
#include "kernel.h"
#include "grid.h"

namespace Manta {
   
//! Kernel: Invert real values, if positive and fluid
KERNEL(idx) 
void InvertCheckFluid (const FlagGrid& flags, Grid<Real>& grid)
{
	if (flags.isFluid(idx) && grid[idx] > 0)
		grid[idx] = 1.0 / grid[idx];
}

//! Kernel: Squared sum over grid
KERNEL(idx, reduce=+) returns(double sum=0)
double GridSumSqr (const Grid<Real>& grid) {
	sum += square((double)grid[idx]);
}

//! Kernel: rotation operator \nabla x v for centered vector fields
KERNEL(bnd=1) 
void CurlOp (const Grid<Vec3>& grid, Grid<Vec3>& dst) {
	Vec3 v = Vec3(0. , 0. , 
			   0.5*((grid(i+1,j,k).y - grid(i-1,j,k).y) - (grid(i,j+1,k).x - grid(i,j-1,k).x)) );
	if(dst.is3D()) {
		v[0] = 0.5*((grid(i,j+1,k).z - grid(i,j-1,k).z) - (grid(i,j,k+1).y - grid(i,j,k-1).y));
		v[1] = 0.5*((grid(i,j,k+1).x - grid(i,j,k-1).x) - (grid(i+1,j,k).z - grid(i-1,j,k).z));
	}
	dst(i,j,k) = v;
};

//! Kernel: divergence operator (from MAC grid)
KERNEL(bnd=1) 
void DivergenceOpMAC(Grid<Real>& div, const MACGrid& grid) {
	Vec3 del = Vec3(grid(i+1,j,k).x, grid(i,j+1,k).y, 0.) - grid(i,j,k); 
	if(grid.is3D()) del[2] += grid(i,j,k+1).z;
	else            del[2]  = 0.;
	div(i,j,k) = del.x + del.y + del.z;
}

//! Kernel: gradient operator for MAC grid
KERNEL(bnd=1)void GradientOpMAC(MACGrid& gradient, const Grid<Real>& grid) {
	Vec3 grad = (Vec3(grid(i,j,k)) - Vec3(grid(i-1,j,k), grid(i,j-1,k), 0. ));
	if(grid.is3D()) grad[2] -= grid(i,j,k-1);
	else            grad[2]  = 0.;
	gradient(i,j,k) = grad;
}

//! Kernel: centered gradient operator 
KERNEL(bnd=1) void GradientOp(Grid<Vec3>& gradient, const Grid<Real>& grid) {
	Vec3 grad = 0.5 * Vec3(        grid(i+1,j,k)-grid(i-1,j,k), 
								   grid(i,j+1,k)-grid(i,j-1,k), 0.);
	if(grid.is3D()) grad[2]= 0.5*( grid(i,j,k+1)-grid(i,j,k-1) );
	gradient(i,j,k) = grad;
}

//! Kernel: Laplace operator
KERNEL (bnd=1) void LaplaceOp(Grid<Real>& laplace, const Grid<Real>& grid) {
	laplace(i, j, k)  = grid(i+1, j, k) - 2.0*grid(i, j, k) + grid(i-1, j, k); 
	laplace(i, j, k) += grid(i, j+1, k) - 2.0*grid(i, j, k) + grid(i, j-1, k); 
	if(grid.is3D()) {
	laplace(i, j, k) += grid(i, j, k+1) - 2.0*grid(i, j, k) + grid(i, j, k-1); }
}

//! Kernel: get component at MAC positions
KERNEL(bnd=1) void GetShiftedComponent(const Grid<Vec3>& grid, Grid<Real>& comp, int dim) {
	Vec3i ishift(i,j,k);
	ishift[dim]--;
	comp(i,j,k) = 0.5*(grid(i,j,k)[dim] + grid(ishift)[dim]);
};

//! Kernel: get component (not shifted)
KERNEL(idx) void GetComponent(const Grid<Vec3>& grid, Grid<Real>& comp, int dim) {
	comp[idx] = grid[idx][dim];
};

//! Kernel: get norm of centered grid
KERNEL(idx) void GridNorm(Grid<Real>& n, const Grid<Vec3>& grid) {
	n[idx] = norm(grid[idx]);
};

//! Kernel: set component (not shifted)
KERNEL(idx) void SetComponent(Grid<Vec3>& grid, const Grid<Real>& comp, int dim) {
	grid[idx][dim] = comp[idx];
};

//! Kernel: compute centered velocity field from MAC
KERNEL(bnd=1) void GetCentered(Grid<Vec3>& center, const MACGrid& vel) {
	Vec3 v = 0.5 * ( vel(i,j,k) + Vec3(vel(i+1,j,k).x, vel(i,j+1,k).y, 0. ) );
	if(vel.is3D()) v[2] += 0.5 * vel(i,j,k+1).z;
	else           v[2]  = 0.;
	center(i,j,k) = v;
};

//! Kernel: compute MAC from centered velocity field
KERNEL(bnd=1) void GetMAC(MACGrid& vel, const Grid<Vec3>& center) {
	Vec3 v = 0.5*(center(i,j,k) + Vec3(center(i-1,j,k).x, center(i,j-1,k).y, 0. ));
	if(vel.is3D()) v[2] += 0.5 * center(i,j,k-1).z; 
	else           v[2]  = 0.;
	vel(i,j,k) = v;
};

//! Fill in the domain boundary cells (i,j,k=0/size-1) from the neighboring cells
KERNEL() void FillInBoundary(Grid<Vec3>& grid, int g) {
	if (i==0) grid(i,j,k) = grid(i+1,j,k);
	if (j==0) grid(i,j,k) = grid(i,j+1,k);
	if (k==0) grid(i,j,k) = grid(i,j,k+1);
	if (i==grid.getSizeX()-1) grid(i,j,k) = grid(i-1,j,k);
	if (j==grid.getSizeY()-1) grid(i,j,k) = grid(i,j-1,k);
	if (k==grid.getSizeZ()-1) grid(i,j,k) = grid(i,j,k-1);
}


// ****************************************************************************

// helper functions for converting mex data to manta grids and back (for matlab integration)

// MAC grids
KERNEL() void kn_conv_mex_in_to_MAC(const double *p_lin_array, MACGrid *p_result)
{
	int ijk = i+j*p_result->getSizeX()+k*p_result->getSizeX()*p_result->getSizeY();
	const int n = p_result->getSizeX() * p_result->getSizeY()*p_result->getSizeZ();

	p_result->get(i,j,k).x = p_lin_array[ijk];
	p_result->get(i,j,k).y = p_lin_array[ijk+n];
	p_result->get(i,j,k).z = p_lin_array[ijk+2*n];
}

KERNEL() void kn_conv_MAC_to_mex_out(const MACGrid *p_mac, double *p_result)
{
	int ijk = i+j*p_mac->getSizeX()+k*p_mac->getSizeX()*p_mac->getSizeY();
	const int n = p_mac->getSizeX() * p_mac->getSizeY()*p_mac->getSizeZ();

	p_result[ijk]     = p_mac->get(i,j,k).x;
	p_result[ijk+n]   = p_mac->get(i,j,k).y;
	p_result[ijk+2*n] = p_mac->get(i,j,k).z;
}

// Vec3 Grids
KERNEL() void kn_conv_mex_in_to_Vec3(const double *p_lin_array, Grid<Vec3> *p_result)
{
	int ijk = i+j*p_result->getSizeX()+k*p_result->getSizeX()*p_result->getSizeY();
	const int n = p_result->getSizeX() * p_result->getSizeY()*p_result->getSizeZ();

	p_result->get(i,j,k).x = p_lin_array[ijk];
	p_result->get(i,j,k).y = p_lin_array[ijk+n];
	p_result->get(i,j,k).z = p_lin_array[ijk+2*n];
}

KERNEL() void kn_conv_Vec3_to_mex_out(const Grid<Vec3> *p_Vec3, double *p_result)
{
	int ijk = i+j*p_Vec3->getSizeX()+k*p_Vec3->getSizeX()*p_Vec3->getSizeY();
	const int n = p_Vec3->getSizeX() * p_Vec3->getSizeY()*p_Vec3->getSizeZ();

	p_result[ijk]     = p_Vec3->get(i,j,k).x;
	p_result[ijk+n]   = p_Vec3->get(i,j,k).y;
	p_result[ijk+2*n] = p_Vec3->get(i,j,k).z;
}

// Real Grids
KERNEL() void kn_conv_mex_in_to_Real(const double *p_lin_array, Grid<Real> *p_result)
{
	int ijk = i+j*p_result->getSizeX()+k*p_result->getSizeX()*p_result->getSizeY();

	p_result->get(i,j,k) = p_lin_array[ijk];
}

KERNEL() void kn_conv_Real_to_mex_out(const Grid<Real> *p_grid, double *p_result)
{
	int ijk = i+j*p_grid->getSizeX()+k*p_grid->getSizeX()*p_grid->getSizeY();

	p_result[ijk] = p_grid->get(i,j,k);
}


} // namespace
#endif
