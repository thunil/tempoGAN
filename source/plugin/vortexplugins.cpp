/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Plugins for using vortex sheet meshes 
 *
 ******************************************************************************/
 
#include <iostream>
#include "vortexsheet.h"
#include "vortexpart.h"
#include "shapes.h"
#include "commonkernels.h"
#include "conjugategrad.h"
#include "randomstream.h"
#include "levelset.h"

using namespace std;

namespace Manta {
	
//! Mark area of mesh inside shape as fixed nodes. 
//! Remove all other fixed nodes if 'exclusive' is set
PYTHON() void markAsFixed(Mesh& mesh, const Shape* shape, bool exclusive=true)
{
	for (int i=0; i<mesh.numNodes(); i++) {
		if (shape->isInside(mesh.nodes(i).pos))
			mesh.nodes(i).flags |= Mesh::NfFixed;
		else if (exclusive)
			mesh.nodes(i).flags &= ~Mesh::NfFixed;
	}
}

//! Adapt texture coordinates of mesh inside shape
//! to obtain an effective inflow effect
PYTHON() void texcoordInflow(VortexSheetMesh& mesh, const Shape* shape, const MACGrid& vel)
{
	static Vec3 t0 = Vec3::Zero;
	
	// get mean velocity
	int cnt=0;
	Vec3 meanV(0.0);
	FOR_IJK(vel) {
		if (shape->isInsideGrid(i,j,k)) {
			cnt++;
			meanV += vel.getCentered(i,j,k);
		}
	}
	meanV /= (Real) cnt;
	t0 -= mesh.getParent()->getDt() * meanV;
	mesh.setReferenceTexOffset(t0);

	// apply mean velocity
	for (int i=0; i<mesh.numNodes(); i++) {
		if (shape->isInside(mesh.nodes(i).pos)) {
			Vec3 tc = mesh.nodes(i).pos + t0;
			mesh.tex1(i) = tc;
			mesh.tex2(i) = tc;
		}
	}
};

//! Init smoke density values of the mesh surface inside source shape
PYTHON() void meshSmokeInflow(VortexSheetMesh& mesh, const Shape* shape, Real amount)
{
	for (int t=0; t<mesh.numTris(); t++) {
		if (shape->isInside(mesh.getFaceCenter(t)))
			mesh.sheet(t).smokeAmount = amount;
	}    
}

KERNEL(idx) 
void KnAcceleration(MACGrid& a, const MACGrid& v1, const MACGrid& v0, const Real idt) { 
	a[idx] = (v1[idx]-v0[idx])*idt; 
}

//! Add vorticity to vortex sheets based on buoyancy
PYTHON() void vorticitySource(VortexSheetMesh& mesh, Vec3 gravity, 
                                                        const MACGrid* vel=NULL, const MACGrid* velOld=NULL,
							Real scale = 0.1, Real maxAmount = 0, Real mult = 1.0)
{
	Real dt = mesh.getParent()->getDt();
	Real dx = mesh.getParent()->getDx();
	MACGrid acceleration(mesh.getParent());
	if (vel)
		KnAcceleration(acceleration, *vel, *velOld, 1.0/dt);
	const Real A= -1.0;
	Real maxV = 0, meanV = 0;
	
	for (int t=0; t<mesh.numTris(); t++) {
		Vec3 fn = mesh.getFaceNormal(t);
		Vec3 source;
		if (vel) {
			Vec3 a = acceleration.getInterpolated(mesh.getFaceCenter(t));        
			source = A*cross(fn, a-gravity) * scale;
		} else {
			source = A*cross(fn, -gravity) * scale;
		}
		
		if (mesh.isTriangleFixed(t)) source = 0;
	
		mesh.sheet(t).vorticity *= mult;
		mesh.sheet(t).vorticity += dt * source / dx;
		// upper limit
		Real v = norm(mesh.sheet(t).vorticity);
		if (maxAmount>0 && v > maxAmount)
			mesh.sheet(t).vorticity *= maxAmount/v;
		
		//stats
		if (v > maxV) maxV = v;
		meanV += v;
	}
	
	cout << "vorticity: max " << maxV << " / mean " << meanV/mesh.numTris() << endl;
}

PYTHON() void smoothVorticity(VortexSheetMesh& mesh, int iter=1, Real sigma=0.2, Real alpha=0.8)
{
	const Real mult = -0.5 / sigma / sigma;
	
	// pre-calculate positions and weights
	vector<Vec3> vort(mesh.numTris()), pos(mesh.numTris());    
	vector<Real> weights(3*mesh.numTris());
	vector<int> index(3*mesh.numTris());
	for(int i=0; i<mesh.numTris(); i++) {
		pos[i] = mesh.getFaceCenter(i);
		mesh.sheet(i).vorticitySmoothed = mesh.sheet(i).vorticity;
	}
	for(int i=0; i<mesh.numTris(); i++) {
		for (int c=0; c<3; c++) {
			int oc = mesh.corners(i,c).opposite;
			if (oc>=0) {
				int t = mesh.corners(oc).tri;
				weights[3*i+c] = exp(normSquare(pos[t]-pos[i])*mult);
				index[3*i+c] = t;
			}
			else {
				weights[3*i+c] = 0;
				index[3*i+c] = 0;
			}
		}        
	}
		
	for (int it=0; it<iter; ++it) {
		// first, preload
		for(int i=0; i<mesh.numTris(); i++) vort[i] = mesh.sheet(i).vorticitySmoothed;
			
		for(int i=0,idx=0; i<mesh.numTris(); i++) {            
			// loop over adjacent tris
			Real sum=1.0f;
			Vec3 v=vort[i];
			for (int c=0;c<3;c++,idx++) {
				Real w = weights[index[idx]];
				v += w*vort[index[idx]];
				sum += w;
			}
			mesh.sheet(i).vorticitySmoothed = v/sum;
		}
	}
	for(int i=0; i<mesh.numTris(); i++) mesh.sheet(i).vorticitySmoothed *= alpha;
}

//! Seed Vortex Particles inside shape with K41 characteristics
PYTHON() void VPseedK41(VortexParticleSystem& system, const Shape* shape, Real strength=0, Real sigma0=0.2, Real sigma1=1.0, Real probability=1.0, Real N=3.0) {
	Grid<Real> temp(system.getParent());
	const Real dt = system.getParent()->getDt();
	static RandomStream rand(3489572);
	Real s0 = pow( (Real)sigma0, (Real)(-N+1.0) );
	Real s1 = pow( (Real)sigma1, (Real)(-N+1.0) );
	
	FOR_IJK(temp) {
		if (shape->isInsideGrid(i,j,k)) {
			if (rand.getReal() < probability*dt) {
				Real p = rand.getReal();
				Real sigma = pow( (1.0-p)*s0 + p*s1, 1./(-N+1.0) );
				Vec3 randDir (rand.getReal(), rand.getReal(), rand.getReal());
				Vec3 posUpd (i+rand.getReal(), j+rand.getReal(), k+rand.getReal());
				normalize(randDir);
				Vec3 vorticity = randDir * strength * pow( (Real)sigma, (Real)(-10./6.+N/2.0) );
				system.add(VortexParticleData(posUpd, vorticity, sigma)); 
			}
		}
	}    
}
		
//! Vortex-in-cell integration
PYTHON() void VICintegration(VortexSheetMesh& mesh, Real sigma, Grid<Vec3>& vel, const FlagGrid& flags,
					  Grid<Vec3>* vorticity=NULL, Real cgMaxIterFac=1.5, Real cgAccuracy=1e-3, Real scale = 0.01, int precondition=0) {
	
	MuTime t0;
	const Real fac = 16.0; // experimental factor to balance out regularization
	
	// if no vort grid is given, use a temporary one
	Grid<Vec3> vortTemp(mesh.getParent());    
	Grid<Vec3>& vort = (vorticity) ? (*vorticity) : (vortTemp);
	vort.clear();
	
	// map vorticity to grid using Peskin kernel
	int sgi = ceil(sigma);
	Real pkfac=M_PI/sigma;
	const int numTris = mesh.numTris();
	for (int t=0; t<numTris; t++) {
		Vec3 pos = mesh.getFaceCenter(t);
		Vec3 v = mesh.sheet(t).vorticity * mesh.getFaceArea(t) * fac;
					
		// inner kernel
		// first, summate 
		Real sum=0;
		for (int i=-sgi; i<sgi; i++) {
			if (pos.x+i < 0 || (int)pos.x+i >= vort.getSizeX()) continue;
			for (int j=-sgi; j<sgi; j++) {
				if (pos.y+j < 0 || (int)pos.y+j >= vort.getSizeY()) continue;            
				for (int k=-sgi; k<sgi; k++) {
					if (pos.z+k < 0 || (int)pos.z+k >= vort.getSizeZ()) continue;                                
					Vec3i cell(pos.x+i, pos.y+j, pos.z+k);
					if (!flags.isFluid(cell)) continue;
					Vec3 d = pos - Vec3(i+0.5+floor(pos.x), j+0.5+floor(pos.y), k+0.5+floor(pos.z));
					Real dl = norm(d);
					if (dl > sigma) continue;
					// precalc Peskin kernel
					sum += 1.0 + cos(dl * pkfac);
				}
			}
		}
		// then, apply normalized kernel
		Real wnorm = 1.0/sum;
		for (int i=-sgi; i<sgi; i++) {
			if (pos.x+i < 0 || (int)pos.x+i >= vort.getSizeX()) continue;
			for (int j=-sgi; j<sgi; j++) {
				if (pos.y+j < 0 || (int)pos.y+j >= vort.getSizeY()) continue;            
				for (int k=-sgi; k<sgi; k++) {
					if (pos.z+k < 0 || (int)pos.z+k >= vort.getSizeZ()) continue;                                
					Vec3i cell(pos.x+i, pos.y+j, pos.z+k);  
					if (!flags.isFluid(cell)) continue;                    
					Vec3 d = pos - Vec3(i+0.5+floor(pos.x), j+0.5+floor(pos.y), k+0.5+floor(pos.z));
					Real dl = norm(d);
					if (dl > sigma) continue;
					Real w = (1.0 + cos(dl * pkfac))*wnorm;
					vort(cell) += v * w;
				}
			}
		}
	}
	
	// Prepare grids for poisson solve
	Grid<Vec3> vortexCurl(mesh.getParent());
	Grid<Real> rhs(mesh.getParent());
	Grid<Real> solution(mesh.getParent());
	Grid<Real> residual(mesh.getParent());
	Grid<Real> search(mesh.getParent());
	Grid<Real> temp1(mesh.getParent());
	Grid<Real> A0(mesh.getParent());
	Grid<Real> Ai(mesh.getParent());
	Grid<Real> Aj(mesh.getParent());
	Grid<Real> Ak(mesh.getParent());
	Grid<Real> pca0(mesh.getParent());
	Grid<Real> pca1(mesh.getParent());
	Grid<Real> pca2(mesh.getParent());
	Grid<Real> pca3(mesh.getParent());
	
	MakeLaplaceMatrix (flags, A0, Ai, Aj, Ak);    
	CurlOp(vort, vortexCurl);    
	
	// Solve vector poisson equation
	for (int c=0; c<3; c++) {
		// construct rhs    
		if (vel.getType() & GridBase::TypeMAC)
			GetShiftedComponent(vortexCurl, rhs, c);
		else
			GetComponent(vortexCurl, rhs, c);
				
		// prepare CG solver
		const int maxIter = (int)(cgMaxIterFac * vel.getSize().max());    
		GridCgInterface *gcg = new GridCg<ApplyMatrix>(solution, rhs, residual, search, flags, temp1, &A0, &Ai, &Aj, &Ak );
		gcg->setAccuracy(cgAccuracy); 
		gcg->setUseL2Norm(true);
		gcg->setICPreconditioner( (GridCgInterface::PreconditionType)precondition, &pca0, &pca1, &pca2, &pca3); 
		
		// iterations
		for (int iter=0; iter<maxIter; iter++) {
			if (!gcg->iterate()) iter=maxIter;
		} 
		debMsg("VICintegration CG iterations:"<<gcg->getIterations()<<", res:"<<gcg->getSigma(), 1);
		delete gcg;
	
		// copy back
		solution *= scale;
		SetComponent(vel, solution, c);
	}
}

//! Obtain density field from levelset with linear gradient of size sigma over the interface
PYTHON() void densityFromLevelset(const LevelsetGrid& phi, Grid<Real>& density, Real value=1.0, Real sigma=1.0) {
	FOR_IJK(phi) {
		// remove boundary
		if (i<2 || j<2 || k<2 || i>=phi.getSizeX()-2 || j>=phi.getSizeY()-2 || k>=phi.getSizeZ()-2)
			density(i,j,k) = 0;
		else if (phi(i,j,k) < -sigma)
			density(i,j,k) = value;
		else if (phi(i,j,k) > sigma)
			density(i,j,k) = 0;
		else
			density(i,j,k) = clamp((Real)(0.5*value/sigma*(1.0-phi(i,j,k))), (Real)0.0, value);
	}    
}

} // namespace
