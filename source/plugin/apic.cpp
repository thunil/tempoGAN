// ----------------------------------------------------------------------------
//
// MantaFlow fluid solver framework
// Copyright 2016-2017 Kiwon Um, Nils Thuerey
//
// This program is free software, distributed under the terms of the
// Apache License, Version 2.0
// http://www.apache.org/licenses/LICENSE-2.0
//
// Affine Particle-In-Cell
//
// ----------------------------------------------------------------------------

#include "particle.h"
#include "grid.h"

namespace Manta {

KERNEL(pts, single)
void knApicMapLinearVec3ToMACGrid(
	const BasicParticleSystem &p, MACGrid &mg, MACGrid &vg, const ParticleDataImpl<Vec3> &vp,
	const ParticleDataImpl<Vec3> &cpx, const ParticleDataImpl<Vec3> &cpy, const ParticleDataImpl<Vec3> &cpz,
	const ParticleDataImpl<int> *ptype, const int exclude)
{
	if (!p.isActive(idx) || (ptype && ((*ptype)[idx] & exclude))) return;
	const IndexInt dX[2] = { 0, vg.getStrideX() };
	const IndexInt dY[2] = { 0, vg.getStrideY() };
	const IndexInt dZ[2] = { 0, vg.getStrideZ() };

	const Vec3 &pos = p[idx].pos, &vel = vp[idx];
	const IndexInt fi = static_cast<IndexInt>(pos.x	   ), fj = static_cast<IndexInt>(pos.y	  ), fk = static_cast<IndexInt>(pos.z	  );
	const IndexInt ci = static_cast<IndexInt>(pos.x-0.5), cj = static_cast<IndexInt>(pos.y-0.5), ck = static_cast<IndexInt>(pos.z-0.5);
	const Real wfi = clamp(pos.x-fi, Real(0), Real(1)), wfj = clamp(pos.y-fj, Real(0), Real(1)), wfk = clamp(pos.z-fk, Real(0), Real(1));
	const Real wci = clamp(Real(pos.x-ci-0.5), Real(0), Real(1)), wcj = clamp(Real(pos.y-cj-0.5), Real(0), Real(1)), wck = clamp(Real(pos.z-ck-0.5), Real(0), Real(1));
	// TODO: check index for safety
	{			// u-face
		const IndexInt gidx = fi*dX[1] + cj*dY[1] + ck*dZ[1];
		const Vec3 gpos(fi, cj+0.5, ck+0.5);
		const Real wi[2] = { Real(1)-wfi, wfi };
		const Real wj[2] = { Real(1)-wcj, wcj };
		const Real wk[2] = { Real(1)-wck, wck };
		for(int i=0; i<2; ++i)
			for(int j=0; j<2; ++j)
				for(int k=0; k<2; ++k) {
					const Real w = wi[i]*wj[j]*wk[k];
					mg[gidx+dX[i]+dY[j]+dZ[k]].x += w;
					vg[gidx+dX[i]+dY[j]+dZ[k]].x += w*vel.x;
					vg[gidx+dX[i]+dY[j]+dZ[k]].x += w*dot(cpx[idx], gpos + Vec3(i, j, k) - pos);
				}
	}
	{			// v-face
		const IndexInt gidx = ci*dX[1] + fj*dY[1] + ck*dZ[1];
		const Vec3 gpos(ci+0.5, fj, ck+0.5);
		const Real wi[2] = { Real(1)-wci, wci };
		const Real wj[2] = { Real(1)-wfj, wfj };
		const Real wk[2] = { Real(1)-wck, wck };
		for(int i=0; i<2; ++i)
			for(int j=0; j<2; ++j)
				for(int k=0; k<2; ++k) {
					const Real w = wi[i]*wj[j]*wk[k];
					mg[gidx+dX[i]+dY[j]+dZ[k]].y += w;
					vg[gidx+dX[i]+dY[j]+dZ[k]].y += w*vel.y;
					vg[gidx+dX[i]+dY[j]+dZ[k]].y += w*dot(cpy[idx], gpos + Vec3(i, j, k) - pos);
				}
	}
	if(!vg.is3D()) return;
	{			// w-face
		const IndexInt gidx = ci*dX[1] + cj*dY[1] + fk*dZ[1];
		const Vec3 gpos(ci+0.5, cj+0.5, fk);
		const Real wi[2] = { Real(1)-wci, wci };
		const Real wj[2] = { Real(1)-wcj, wcj };
		const Real wk[2] = { Real(1)-wfk, wfk };
		for(int i=0; i<2; ++i)
			for(int j=0; j<2; ++j)
				for(int k=0; k<2; ++k) {
					const Real w = wi[i]*wj[j]*wk[k];
					mg[gidx+dX[i]+dY[j]+dZ[k]].z += w;
					vg[gidx+dX[i]+dY[j]+dZ[k]].z += w*vel.z;
					vg[gidx+dX[i]+dY[j]+dZ[k]].z += w*dot(cpz[idx], gpos + Vec3(i, j, k) - pos);
				}
	}
}

PYTHON()
void apicMapPartsToMAC(
	const FlagGrid &flags, MACGrid &vel,
	const BasicParticleSystem &parts, const ParticleDataImpl<Vec3> &partVel,
	const ParticleDataImpl<Vec3> &cpx, const ParticleDataImpl<Vec3> &cpy, const ParticleDataImpl<Vec3> &cpz,
	MACGrid *mass=NULL, const ParticleDataImpl<int> *ptype=NULL, const int exclude=0)
{
	// affine map
	// let's assume that the particle mass is constant, 1.0
	const bool freeMass = !mass;
	if(!mass) mass = new MACGrid(flags.getParent());
	else mass->clear();

	vel.clear();
	knApicMapLinearVec3ToMACGrid(parts, *mass, vel, partVel, cpx, cpy, cpz, ptype, exclude);
	mass->stomp(VECTOR_EPSILON);
	vel.safeDivide(*mass);

	if(freeMass) delete mass;
}

KERNEL(pts)
void knApicMapLinearMACGridToVec3(
	ParticleDataImpl<Vec3> &vp, ParticleDataImpl<Vec3> &cpx, ParticleDataImpl<Vec3> &cpy, ParticleDataImpl<Vec3> &cpz,
	const BasicParticleSystem &p, const MACGrid &vg, const FlagGrid &flags,
	const ParticleDataImpl<int> *ptype, const int exclude)
{
	if (!p.isActive(idx) || (ptype && ((*ptype)[idx] & exclude))) return;

	vp[idx] = cpx[idx] = cpy[idx] = cpz[idx] = Vec3(Real(0));
	const IndexInt dX[2] = {0, vg.getStrideX()}, dY[2] = {0, vg.getStrideY()}, dZ[2] = {0, vg.getStrideZ()};
	const Real gw[2] = {-Real(1), Real(1)};

	const Vec3 &pos = p[idx].pos;
	const IndexInt fi = static_cast<IndexInt>(pos.x	   ), fj = static_cast<IndexInt>(pos.y	  ), fk = static_cast<IndexInt>(pos.z	 );
	const IndexInt ci = static_cast<IndexInt>(pos.x-0.5), cj = static_cast<IndexInt>(pos.y-0.5), ck = static_cast<IndexInt>(pos.z-0.5);
	const Real wfi = clamp(pos.x-fi, Real(0), Real(1)), wfj = clamp(pos.y-fj, Real(0), Real(1)), wfk = clamp(pos.z-fk, Real(0), Real(1));
	const Real wci = clamp(Real(pos.x-ci-0.5), Real(0), Real(1)), wcj = clamp(Real(pos.y-cj-0.5), Real(0), Real(1)), wck = clamp(Real(pos.z-ck-0.5), Real(0), Real(1));
	// TODO: check index for safety
	{			// u
		const IndexInt gidx = fi*dX[1] + cj*dY[1] + ck*dZ[1];
		const Real wx[2] = { Real(1)-wfi, wfi };
		const Real wy[2] = { Real(1)-wcj, wcj };
		const Real wz[2] = { Real(1)-wck, wck };
		for(int i=0; i<2; ++i)
			for(int j=0; j<2; ++j)
				for(int k=0; k<2; ++k) {
					const IndexInt vidx = gidx+dX[i]+dY[j]+dZ[k];
					Real vgx = vg[vidx].x;
					vp[idx].x  += wx[i]*wy[j]*wz[k]*vgx;
					cpx[idx].x += gw[i]*wy[j]*wz[k]*vgx;
					cpx[idx].y += wx[i]*gw[j]*wz[k]*vgx;
					cpx[idx].z += wx[i]*wy[j]*gw[k]*vgx;
				}
	}
	{			// v
		const IndexInt gidx = ci*dX[1] + fj*dY[1] + ck*dZ[1];
		const Real wx[2] = { Real(1)-wci, wci };
		const Real wy[2] = { Real(1)-wfj, wfj };
		const Real wz[2] = { Real(1)-wck, wck };
		for(int i=0; i<2; ++i)
			for(int j=0; j<2; ++j)
				for(int k=0; k<2; ++k) {
					const IndexInt vidx = gidx+dX[i]+dY[j]+dZ[k];
					Real vgy = vg[vidx].y;
					vp[idx].y  += wx[i]*wy[j]*wz[k]*vgy;
					cpy[idx].x += gw[i]*wy[j]*wz[k]*vgy;
					cpy[idx].y += wx[i]*gw[j]*wz[k]*vgy;
					cpy[idx].z += wx[i]*wy[j]*gw[k]*vgy;
				}
	}
	if(!vg.is3D()) return;
	{			// w
		const IndexInt gidx = ci*dX[1] + cj*dY[1] + fk*dZ[1];
		const Real wx[2] = { Real(1)-wci, wci };
		const Real wy[2] = { Real(1)-wcj, wcj };
		const Real wz[2] = { Real(1)-wfk, wfk };
		for(int i=0; i<2; ++i)
			for(int j=0; j<2; ++j)
				for(int k=0; k<2; ++k) {
					const IndexInt vidx = gidx+dX[i]+dY[j]+dZ[k];
					Real vgz = vg[vidx].z;
					vp[idx].z  += wx[i]*wy[j]*wz[k]*vgz;
					cpz[idx].x += gw[i]*wy[j]*wz[k]*vgz;
					cpz[idx].y += wx[i]*gw[j]*wz[k]*vgz;
					cpz[idx].z += wx[i]*wy[j]*gw[k]*vgz;
				}
	}
}

PYTHON()
void apicMapMACGridToParts(
	ParticleDataImpl<Vec3> &partVel, ParticleDataImpl<Vec3> &cpx, ParticleDataImpl<Vec3> &cpy, ParticleDataImpl<Vec3> &cpz,
	const BasicParticleSystem &parts, const MACGrid &vel, const FlagGrid &flags,
	const ParticleDataImpl<int> *ptype=NULL, const int exclude=0)
{
	knApicMapLinearMACGridToVec3(partVel, cpx, cpy, cpz, parts, vel, flags, ptype, exclude);
}

} // namespace
