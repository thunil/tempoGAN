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
#include "kernel.h"
#include "conjugategrad.h"
#include "multigrid.h"

using namespace std;
namespace Manta {

//! Preconditioner for CG solver
// - None: Use standard CG
// - MIC: Modified incomplete Cholesky preconditioner
// - MGDynamic: Multigrid preconditioner, rebuilt for each solve
// - MGStatic: Multigrid preconditioner, built only once (faster than
//       MGDynamic, but works only if Poisson equation does not change)
enum Preconditioner { PcNone = 0, PcMIC = 1, PcMGDynamic = 2, PcMGStatic = 3 };

inline static Real surfTensHelper(const IndexInt idx, const int offset, const Grid<Real> &phi, const Grid<Real> &curv, const Real surfTens, const Real gfClamp);

//! Kernel: Construct the right-hand side of the poisson equation
KERNEL(bnd=1, reduce=+) returns(int cnt=0) returns(double sum=0)
void MakeRhs (const FlagGrid& flags, Grid<Real>& rhs, const MACGrid& vel,
                          const Grid<Real>* perCellCorr, const MACGrid* fractions,
			  // note - all of the following are necessary for surface tension
              const Grid<Real> *phi, const Grid<Real> *curv, const Real surfTens, const Real gfClamp)
{
	if (!flags.isFluid(i,j,k)) {
		rhs(i,j,k) = 0;
		return;
	}

	// compute negative divergence 
	// no flag checks: assumes vel at obstacle interfaces is set to zero
	Real set(0);
	if(!fractions) {
		set =               vel(i,j,k).x - vel(i+1,j,k).x + 
				 			vel(i,j,k).y - vel(i,j+1,k).y; 
		if(vel.is3D()) set+=vel(i,j,k).z - vel(i,j,k+1).z;
	}else{
		set =               (*fractions)(i,j,k).x * vel(i,j,k).x - (*fractions)(i+1,j,k).x * vel(i+1,j,k).x + 
							(*fractions)(i,j,k).y * vel(i,j,k).y - (*fractions)(i,j+1,k).y * vel(i,j+1,k).y; 
		if(vel.is3D()) set+=(*fractions)(i,j,k).z * vel(i,j,k).z - (*fractions)(i,j,k+1).z * vel(i,j,k+1).z;
	}

	// compute surface tension effect (optional)
	if(phi && curv) {
		const IndexInt idx = flags.index(i,j,k);
		const int X = flags.getStrideX(), Y = flags.getStrideY(), Z = flags.getStrideZ();
		if(flags.isEmpty(i-1,j,k))     set += surfTensHelper(idx, -X, *phi, *curv, surfTens, gfClamp);
		if(flags.isEmpty(i+1,j,k))	   set += surfTensHelper(idx, +X, *phi, *curv, surfTens, gfClamp);
		if(flags.isEmpty(i,j-1,k))	   set += surfTensHelper(idx, -Y, *phi, *curv, surfTens, gfClamp);
		if(flags.isEmpty(i,j+1,k))	   set += surfTensHelper(idx, +Y, *phi, *curv, surfTens, gfClamp);
		if(vel.is3D()) {
			if(flags.isEmpty(i,j,k-1)) set += surfTensHelper(idx, -Z, *phi, *curv, surfTens, gfClamp);
			if(flags.isEmpty(i,j,k+1)) set += surfTensHelper(idx, +Z, *phi, *curv, surfTens, gfClamp);
		}
	}
	
	// per cell divergence correction (optional)
	if(perCellCorr) 
		set += perCellCorr->get(i,j,k);
	
	// obtain sum, cell count
	sum += set;
	cnt++;
	
	rhs(i,j,k) = set;
}

//! Kernel: make velocity divergence free by subtracting pressure gradient
KERNEL(bnd = 1)
void knCorrectVelocity(const FlagGrid& flags, MACGrid& vel, const Grid<Real>& pressure)
{
	IndexInt idx = flags.index(i,j,k);
	if (flags.isFluid(idx))
	{
		if (flags.isFluid(i-1,j,k)) vel[idx].x -= (pressure[idx] - pressure(i-1,j,k));
		if (flags.isFluid(i,j-1,k)) vel[idx].y -= (pressure[idx] - pressure(i,j-1,k));
		if (flags.is3D() && flags.isFluid(i,j,k-1)) vel[idx].z -= (pressure[idx] - pressure(i,j,k-1));
 
		if (flags.isEmpty(i-1,j,k)) vel[idx].x -= pressure[idx];
		if (flags.isEmpty(i,j-1,k)) vel[idx].y -= pressure[idx];
		if (flags.is3D() && flags.isEmpty(i,j,k-1)) vel[idx].z -= pressure[idx];
	}
	else if (flags.isEmpty(idx)&&!flags.isOutflow(idx)) // don't change velocities in outflow cells
	{
		if (flags.isFluid(i-1,j,k)) vel[idx].x += pressure(i-1,j,k);
		else                        vel[idx].x  = 0.f;
		if (flags.isFluid(i,j-1,k)) vel[idx].y += pressure(i,j-1,k);
		else                        vel[idx].y  = 0.f;
		if (flags.is3D() ) {
		if (flags.isFluid(i,j,k-1)) vel[idx].z += pressure(i,j,k-1);
		else                        vel[idx].z  = 0.f;
		}
	}
}

// *****************************************************************************
// Ghost fluid helpers

// calculate fraction filled with liquid (note, assumes inside value is < outside!)
inline static Real thetaHelper(Real inside, Real outside)
{
	Real denom = inside-outside;
	if (denom > -1e-04) return 0.5; // should always be neg, and large enough...
	return std::max(Real(0), std::min(Real(1), inside/denom));
}

// calculate ghost fluid factor, cell at idx should be a fluid cell
inline static Real ghostFluidHelper(IndexInt idx, int offset, const Grid<Real> &phi, Real gfClamp)
{
	Real alpha = thetaHelper(phi[idx], phi[idx+offset]);
	if (alpha < gfClamp) return alpha = gfClamp;
	return (1-(1/alpha)); 
}

inline static Real surfTensHelper(const IndexInt idx, const int offset, const Grid<Real> &phi, const Grid<Real> &curv, const Real surfTens, const Real gfClamp)
{
	return surfTens*(curv[idx+offset] - ghostFluidHelper(idx, offset, phi, gfClamp) * curv[idx]);
}

//! Kernel: Adapt A0 for ghost fluid
KERNEL(bnd=1) 
void ApplyGhostFluidDiagonal(Grid<Real> &A0, const FlagGrid &flags, const Grid<Real> &phi, Real gfClamp)
{
	const int X = flags.getStrideX(), Y = flags.getStrideY(), Z = flags.getStrideZ();
	IndexInt idx = flags.index(i,j,k);
	if (!flags.isFluid(idx)) return;

	if (flags.isEmpty(i-1,j,k)) A0[idx] -= ghostFluidHelper(idx, -X, phi, gfClamp);
	if (flags.isEmpty(i+1,j,k)) A0[idx] -= ghostFluidHelper(idx, +X, phi, gfClamp);
	if (flags.isEmpty(i,j-1,k)) A0[idx] -= ghostFluidHelper(idx, -Y, phi, gfClamp);
	if (flags.isEmpty(i,j+1,k)) A0[idx] -= ghostFluidHelper(idx, +Y, phi, gfClamp);
	if (flags.is3D()) {
		if (flags.isEmpty(i,j,k-1)) A0[idx] -= ghostFluidHelper(idx, -Z, phi, gfClamp);
		if (flags.isEmpty(i,j,k+1)) A0[idx] -= ghostFluidHelper(idx, +Z, phi, gfClamp);
	}
}

//! Kernel: Apply velocity update: ghost fluid contribution
KERNEL(bnd=1)
void knCorrectVelocityGhostFluid(MACGrid &vel, const FlagGrid &flags, const Grid<Real> &pressure, const Grid<Real> &phi, Real gfClamp,
                               const Grid<Real> *curv, const Real surfTens)
{
	const IndexInt X = flags.getStrideX(), Y = flags.getStrideY(), Z = flags.getStrideZ();
	const IndexInt idx = flags.index(i,j,k);
	if (flags.isFluid(idx)) {
		if (flags.isEmpty(i-1,j,k)) vel[idx][0] += pressure[idx] * ghostFluidHelper(idx, -X, phi, gfClamp);
		if (flags.isEmpty(i,j-1,k)) vel[idx][1] += pressure[idx] * ghostFluidHelper(idx, -Y, phi, gfClamp);
		if (flags.is3D() && flags.isEmpty(i,j,k-1)) vel[idx][2] += pressure[idx] * ghostFluidHelper(idx, -Z, phi, gfClamp);
	}
	else if (flags.isEmpty(idx) && !flags.isOutflow(idx)) // do not change velocities in outflow cells
	{
		if (flags.isFluid(i-1,j,k)) vel[idx][0] -= pressure(i-1,j,k) * ghostFluidHelper(idx-X, +X, phi, gfClamp);
		else                        vel[idx].x  = 0.f;
		if (flags.isFluid(i,j-1,k)) vel[idx][1] -= pressure(i,j-1,k) * ghostFluidHelper(idx-Y, +Y, phi, gfClamp);
		else                        vel[idx].y  = 0.f;
		if (flags.is3D() ) {
		if (flags.isFluid(i,j,k-1)) vel[idx][2] -= pressure(i,j,k-1) * ghostFluidHelper(idx-Z, +Z, phi, gfClamp);
		else                        vel[idx].z  = 0.f;
		}
	}

	if(curv) {
		if(flags.isFluid(idx))	{
			if(flags.isEmpty(i-1,j,k))                 vel[idx].x += surfTensHelper(idx, -X, phi, *curv, surfTens, gfClamp);
			if(flags.isEmpty(i,j-1,k))                 vel[idx].y += surfTensHelper(idx, -Y, phi, *curv, surfTens, gfClamp);
			if(flags.is3D() && flags.isEmpty(i,j,k-1)) vel[idx].z += surfTensHelper(idx, -Z, phi, *curv, surfTens, gfClamp);
		} 
		else if (flags.isEmpty(idx) && !flags.isOutflow(idx)) // do not change velocities in outflow cells
		{
			vel[idx].x -= (flags.isFluid(i-1,j,k)) ? surfTensHelper(idx-X, +X, phi, *curv, surfTens, gfClamp) : 0.f;
			vel[idx].y -= (flags.isFluid(i,j-1,k)) ? surfTensHelper(idx-Y, +Y, phi, *curv, surfTens, gfClamp) : 0.f;
			if(flags.is3D()) {
			vel[idx].z -= (flags.isFluid(i,j,k-1)) ? surfTensHelper(idx-Z, +Z, phi, *curv, surfTens, gfClamp) : 0.f;
			}
		}
	}
}


// improve behavior of clamping for large time steps:
inline static Real ghostFluidWasClamped(IndexInt idx, int offset, const Grid<Real> &phi, Real gfClamp)
{
	Real alpha = thetaHelper(phi[idx], phi[idx+offset]);
	if (alpha < gfClamp) return true;
	return false;
}

KERNEL(bnd=1)
void knReplaceClampedGhostFluidVels(MACGrid &vel, const FlagGrid &flags,
		const Grid<Real> &pressure, const Grid<Real> &phi, Real gfClamp )
{
	const IndexInt idx = flags.index(i,j,k);
	const IndexInt X   = flags.getStrideX(), Y = flags.getStrideY(), Z = flags.getStrideZ();
	if (!flags.isEmpty(idx)) return;

	if( (flags.isFluid(i-1,j,k)) && ( ghostFluidWasClamped(idx-X, +X, phi, gfClamp) ) )
		vel[idx][0] = vel[idx-X][0];
	if( (flags.isFluid(i,j-1,k)) && ( ghostFluidWasClamped(idx-Y, +Y, phi, gfClamp) ) )
		vel[idx][1] = vel[idx-Y][1];
	if( flags.is3D() &&
	  ( (flags.isFluid(i,j,k-1)) && ( ghostFluidWasClamped(idx-Z, +Z, phi, gfClamp) ) ))
		vel[idx][2] = vel[idx-Z][2];

	if( (flags.isFluid(i+1,j,k)) && ( ghostFluidWasClamped(idx+X, -X, phi, gfClamp)) )
		vel[idx][0] = vel[idx+X][0];
	if( (flags.isFluid(i,j+1,k)) && ( ghostFluidWasClamped(idx+Y, -Y, phi, gfClamp)) )
		vel[idx][1] = vel[idx+Y][1];
	if( flags.is3D() && 
	   (flags.isFluid(i,j,k+1))  && ( ghostFluidWasClamped(idx+Z, -Z, phi, gfClamp)) )
		vel[idx][2] = vel[idx+Z][2];
}

//! Kernel: Compute min value of Real grid
KERNEL(idx, reduce=+) returns(int numEmpty=0)
int CountEmptyCells(const FlagGrid& flags) {
	if (flags.isEmpty(idx) ) numEmpty++;
}

// *****************************************************************************
// Misc helpers

//! Change 'A' and 'rhs' such that pressure at 'fixPidx' is fixed to 'value'
void fixPressure (int fixPidx, Real value, Grid<Real>& rhs, Grid<Real>& A0, Grid<Real>& Ai, Grid<Real>& Aj, Grid<Real>& Ak)
{
	// Bring to rhs at neighbors
	rhs[fixPidx + Ai.getStrideX()] -= Ai[fixPidx] * value;
	rhs[fixPidx + Aj.getStrideY()] -= Aj[fixPidx] * value;
	rhs[fixPidx - Ai.getStrideX()] -= Ai[fixPidx - Ai.getStrideX()] * value;
	rhs[fixPidx - Aj.getStrideY()] -= Aj[fixPidx - Aj.getStrideY()] * value;
	if (rhs.is3D()) {
		rhs[fixPidx + Ak.getStrideZ()] -= Ak[fixPidx] * value;
		rhs[fixPidx - Ak.getStrideZ()] -= Ak[fixPidx - Ak.getStrideZ()] * value;
	}
	
	// Trivialize equation at 'fixPidx' to: pressure[fixPidx] = value
	rhs[fixPidx] = value;
	A0[fixPidx] = Real(1);
	Ai[fixPidx] = Aj[fixPidx] = Ak[fixPidx] = Real(0);
	Ai[fixPidx - Ai.getStrideX()] = Real(0);
	Aj[fixPidx - Aj.getStrideY()] = Real(0);
	if (rhs.is3D()) { Ak[fixPidx - Ak.getStrideZ()] = Real(0); }
}


// for "static" MG mode, keep one MG data structure per fluid solver
// leave cleanup to OS/user if nonzero at program termination (PcMGStatic mode)
// alternatively, manually release in scene file with releaseMG
static std::map<FluidSolver*, GridMg*> gMapMG;

PYTHON() void releaseMG(FluidSolver* solver=nullptr) {
	// release all?
	if(!solver) {
		for( std::map<FluidSolver*, GridMg*>::iterator it = gMapMG.begin(); it != gMapMG.end(); it++) {
			if(it->first != nullptr) releaseMG(it->first);
		}
		return;
	}

	GridMg* mg = gMapMG[solver];
	if(mg) {
		delete mg; 
		gMapMG[solver] = nullptr;
	}
}


// *****************************************************************************
// Main pressure solve

// Note , all three pressure solve helper functions take
// identical parameters, apart from the RHS grid (and different const values)


//! Compute rhs for pressure solve
PYTHON() void computePressureRhs(Grid<Real>& rhs, const MACGrid& vel, 
	const Grid<Real>& pressure, const FlagGrid& flags, Real cgAccuracy = 1e-3,
    const Grid<Real>* phi = 0,
    const Grid<Real>* perCellCorr = 0,
    const MACGrid* fractions = 0,
    Real gfClamp = 1e-04,
    Real cgMaxIterFac = 1.5,
    bool precondition = true, // Deprecated, use preconditioner instead
	int preconditioner = PcMIC,
	bool enforceCompatibility = false,
    bool useL2Norm = false, 
	bool zeroPressureFixing = false,
	const Grid<Real> *curv = NULL,
	const Real surfTens = 0. )
{
	// compute divergence and init right hand side
	MakeRhs kernMakeRhs (flags, rhs, vel, perCellCorr, fractions,  phi, curv, surfTens, gfClamp );
	
	if (enforceCompatibility)
		rhs += (Real)(-kernMakeRhs.sum / (Real)kernMakeRhs.cnt);
}
	
//! Build and solve pressure system of equations
//! perCellCorr: a divergence correction for each cell, optional
//! fractions: for 2nd order obstacle boundaries, optional
//! gfClamp: clamping threshold for ghost fluid method
//! cgMaxIterFac: heuristic to determine maximal number of CG iteations, increase for more accurate solutions
//! preconditioner: MIC, or MG (see Preconditioner enum)
//! useL2Norm: use max norm by default, can be turned to L2 here
//! zeroPressureFixing: remove null space by fixing a single pressure value, needed for MG 
//! curv: curvature for surface tension effects
//! surfTens: surface tension coefficient
//! retRhs: return RHS divergence, e.g., for debugging; optional
PYTHON() void solvePressureSystem(Grid<Real>& rhs, MACGrid& vel, 
	Grid<Real>& pressure, const FlagGrid& flags, Real cgAccuracy = 1e-3,
	const Grid<Real>* phi = 0,
    const Grid<Real>* perCellCorr = 0,
    const MACGrid* fractions = 0,
    Real gfClamp = 1e-04,
    Real cgMaxIterFac = 1.5,
    bool precondition = true, // Deprecated, use preconditioner instead
	int preconditioner = PcMIC,
	bool enforceCompatibility = false,
    bool useL2Norm = false, 
	bool zeroPressureFixing = false,
	const Grid<Real> *curv = NULL,
	const Real surfTens = 0.)
{
	if (precondition==false) preconditioner = PcNone; // for backwards compatibility

	// reserve temp grids
	FluidSolver* parent = flags.getParent();
	Grid<Real> residual(parent);
	Grid<Real> search(parent);
	Grid<Real> A0(parent);
	Grid<Real> Ai(parent);
	Grid<Real> Aj(parent);
	Grid<Real> Ak(parent);
	Grid<Real> tmp(parent);
		
	// setup matrix and boundaries 
	MakeLaplaceMatrix(flags, A0, Ai, Aj, Ak, fractions);

	if (phi) {
		ApplyGhostFluidDiagonal(A0, flags, *phi, gfClamp);
	}
	
	// check whether we need to fix some pressure value...
	// (manually enable, or automatically for high accuracy, can cause asymmetries otherwise)
	if(zeroPressureFixing || cgAccuracy<1e-07) 
	{
		if(FLOATINGPOINT_PRECISION==1) debMsg("Warning - high CG accuracy with single-precision floating point accuracy might not converge...", 2);

		int numEmpty = CountEmptyCells(flags);
		IndexInt fixPidx = -1;
		if(numEmpty==0) {
			// Determine appropriate fluid cell for pressure fixing
			// 1) First check some preferred positions for approx. symmetric zeroPressureFixing
			Vec3i topCenter(flags.getSizeX() / 2, flags.getSizeY() - 1, flags.is3D() ? flags.getSizeZ() / 2 : 0);
			Vec3i preferredPos [] = { topCenter, 
				                      topCenter - Vec3i(0,1,0), 
				                      topCenter - Vec3i(0,2,0) };
			
			for (Vec3i pos : preferredPos) {
				if(flags.isFluid(pos)) {
					fixPidx = flags.index(pos);
					break;
				}
			}

			// 2) Then search whole domain
			if (fixPidx == -1) {
				FOR_IJK_BND(flags,1) {
					if(flags.isFluid(i,j,k)) {
						fixPidx = flags.index(i,j,k);
						// break FOR_IJK_BND loop
						i = flags.getSizeX()-1; 
						j = flags.getSizeY()-1;
						k = __kmax;
					}
				}
			}
			//debMsg("No empty cells! Fixing pressure of cell "<<fixPidx<<" to zero",1);
		}
		if(fixPidx>=0) {
			fixPressure(fixPidx, Real(0), rhs, A0, Ai, Aj, Ak);
			static bool msgOnce = false;
			if(!msgOnce) { debMsg("Pinning pressure of cell "<<fixPidx<<" to zero", 2); msgOnce=true; }
		}
	}

	// CG setup
	// note: the last factor increases the max iterations for 2d, which right now can't use a preconditioner 
	GridCgInterface *gcg;
	if (vel.is3D())
		gcg = new GridCg<ApplyMatrix>  (pressure, rhs, residual, search, flags, tmp, &A0, &Ai, &Aj, &Ak );
	else
		gcg = new GridCg<ApplyMatrix2D>(pressure, rhs, residual, search, flags, tmp, &A0, &Ai, &Aj, &Ak );
	
	gcg->setAccuracy( cgAccuracy ); 
	gcg->setUseL2Norm( useL2Norm );

	int maxIter = 0;
	
	Grid<Real> *pca0 = nullptr, *pca1 = nullptr, *pca2 = nullptr, *pca3 = nullptr;
	GridMg* pmg = nullptr;

	// optional preconditioning	
	if (preconditioner == PcNone || preconditioner == PcMIC) {			
		maxIter = (int)(cgMaxIterFac * flags.getSize().max()) * (flags.is3D() ? 1 : 4);

		pca0 = new Grid<Real>(parent);
		pca1 = new Grid<Real>(parent);
		pca2 = new Grid<Real>(parent);
		pca3 = new Grid<Real>(parent);

		gcg->setICPreconditioner( preconditioner == PcMIC ? GridCgInterface::PC_mICP : GridCgInterface::PC_None, 
			pca0, pca1, pca2, pca3);
	} else if (preconditioner == PcMGDynamic || preconditioner == PcMGStatic) {
		maxIter = 100;

		pmg = gMapMG[parent];
		if (!pmg) {
			pmg = new GridMg(pressure.getSize());
			gMapMG[parent] = pmg;
		}

		gcg->setMGPreconditioner( GridCgInterface::PC_MGP, pmg);
	}

	// CG solve
	for (int iter=0; iter<maxIter; iter++) {
		if (!gcg->iterate()) iter=maxIter;
		debMsg("FluidSolver::solvePressure iteration "<<iter<<", residual: "<<gcg->getResNorm(), 9);
	} 
	debMsg("FluidSolver::solvePressure iterations:"<<gcg->getIterations()<<", residual:"<<gcg->getResNorm(), 2);

	// Cleanup
	if (gcg)  delete gcg;
	if (pca0) delete pca0;
	if (pca1) delete pca1;
	if (pca2) delete pca2;
	if (pca3) delete pca3;

	// PcMGDynamic: always delete multigrid solver after use
	// PcMGStatic: keep multigrid solver for next solve
	if (pmg && preconditioner==PcMGDynamic) releaseMG(parent);
}

//! Apply pressure gradient to make velocity field divergence free
PYTHON() void correctVelocity(MACGrid& vel, Grid<Real>& pressure, const FlagGrid& flags, Real cgAccuracy = 1e-3,
    const Grid<Real>* phi = 0,
    const Grid<Real>* perCellCorr = 0,
    const MACGrid* fractions = 0,
    Real gfClamp = 1e-04,
    Real cgMaxIterFac = 1.5,
    bool precondition = true, // Deprecated, use preconditioner instead
	int preconditioner = PcMIC,
	bool enforceCompatibility = false,
    bool useL2Norm = false, 
	bool zeroPressureFixing = false,
	const Grid<Real> *curv = NULL,
	const Real surfTens = 0. )
{
	knCorrectVelocity(flags, vel, pressure ); 
	if (phi) {
		knCorrectVelocityGhostFluid (vel, flags, pressure, *phi, gfClamp,  curv, surfTens );
		// improve behavior of clamping for large time steps:
		knReplaceClampedGhostFluidVels (vel, flags, pressure, *phi, gfClamp);
	}
}

//! Perform pressure projection of the velocity grid, calls
//! all three pressure helper functions in a row.
PYTHON() void solvePressure(MACGrid& vel, Grid<Real>& pressure, const FlagGrid& flags, Real cgAccuracy = 1e-3,
    const Grid<Real>* phi = 0,
    const Grid<Real>* perCellCorr = 0,
    const MACGrid* fractions = 0,
    Real gfClamp = 1e-04,
    Real cgMaxIterFac = 1.5,
    bool precondition = true, // Deprecated, use preconditioner instead
	int preconditioner = PcMIC,
	bool enforceCompatibility = false,
    bool useL2Norm = false, 
	bool zeroPressureFixing = false,
	const Grid<Real> *curv = NULL,
	const Real surfTens = 0.,
	Grid<Real>* retRhs = NULL )
{
	Grid<Real> rhs(vel.getParent());

	computePressureRhs(rhs, vel, pressure, flags, cgAccuracy,
		phi, perCellCorr, fractions, gfClamp,
		cgMaxIterFac, precondition, preconditioner, enforceCompatibility,
		useL2Norm, zeroPressureFixing, curv, surfTens);

	solvePressureSystem(rhs, vel, pressure, flags, cgAccuracy,
		phi, perCellCorr, fractions, gfClamp,
		cgMaxIterFac, precondition, preconditioner, enforceCompatibility,
		useL2Norm, zeroPressureFixing, curv, surfTens);

	correctVelocity(vel, pressure, flags, cgAccuracy,
		phi, perCellCorr, fractions, gfClamp,
		cgMaxIterFac, precondition, preconditioner, enforceCompatibility,
		useL2Norm, zeroPressureFixing, curv, surfTens);

	// optionally , return RHS
	if(retRhs) {
		retRhs->copyFrom( rhs );
	}
}

} // end namespace

