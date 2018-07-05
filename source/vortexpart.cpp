/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Vortex particles
 * (warning, the vortex methods are currently experimental, and not fully supported!)
 *
 ******************************************************************************/

#include "vortexpart.h"
#include "integrator.h"
#include "mesh.h"

using namespace std;
namespace Manta {

// vortex particle effect: (cyl coord around wp)
// u = -|wp|*rho*exp( (-rho^2-z^2)/(2sigma^2) ) e_phi
inline Vec3 VortexKernel(const Vec3& p, const vector<VortexParticleData>& vp, Real scale) {
	Vec3 u(0.0);
	for (size_t i=0; i<vp.size(); i++) {
		if (vp[i].flag & ParticleBase::PDELETE) continue;
		
		// cutoff radius
		const Vec3 r = p - vp[i].pos;
		const Real rlen2 = normSquare(r);   
		const Real sigma2 = square(vp[i].sigma);
		if (rlen2 > 6.0 * sigma2 || rlen2 < 1e-8) continue;
		
		// split vortex strength
		Vec3 vortNorm = vp[i].vorticity;
		Real strength = normalize(vortNorm) * scale;
	
		// transform in cylinder coordinate system
		const Real rlen = sqrt(rlen2);
		const Real z = dot(r, vortNorm);
		const Vec3 ePhi = cross(r, vortNorm) / rlen;
		const Real rho2 = rlen2 - z*z;
	
		Real vortex = 0;
		if (rho2 > 1e-10) {
			// evaluate Kernel      
			vortex = strength * sqrt(rho2) * exp (rlen2 * -0.5/sigma2);  
		}
		u += vortex * ePhi;
	}
	return u;
}

KERNEL(pts) returns(vector<Vec3> u(size))
vector<Vec3> KnVpAdvectMesh(vector<Node>& nodes, const vector<VortexParticleData>& vp, Real scale) {
	if (nodes[idx].flags & Mesh::NfFixed)
		u[idx] = 0.0;
	else
		u[idx] = VortexKernel(nodes[idx].pos, vp, scale);
}

KERNEL(pts) returns(vector<Vec3> u(size))
vector<Vec3> KnVpAdvectSelf(vector<VortexParticleData>& vp, Real scale) {
	if (vp[idx].flag & ParticleBase::PDELETE) 
		u[idx] = 0.0;
	else
		u[idx] = VortexKernel(vp[idx].pos, vp, scale);
}
	
VortexParticleSystem::VortexParticleSystem(FluidSolver* parent) :
	ParticleSystem<VortexParticleData>(parent)
{ 
}

void VortexParticleSystem::advectSelf(Real scale, int integrationMode) {
	KnVpAdvectSelf kernel(mData, scale* getParent()->getDt());
	integratePointSet( kernel, integrationMode);    
}

void VortexParticleSystem::applyToMesh(Mesh& mesh, Real scale, int integrationMode) {
	KnVpAdvectMesh kernel(mesh.getNodeData(), mData, scale* getParent()->getDt());
	integratePointSet( kernel, integrationMode);    
}

ParticleBase* VortexParticleSystem::clone() {
	VortexParticleSystem* nm = new VortexParticleSystem(getParent());
	compress();
	
	nm->mData = mData;
	nm->setName(getName());
	return nm;
}

	

} // namespace
