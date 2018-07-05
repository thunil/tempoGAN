/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2017 Kiwon Um, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * tensorflor/numpy plugins, mostly for MLFLIP for now  (only compiled if NUMPY is enabled)
 *
 ******************************************************************************/

#include "levelset.h"
#include "commonkernels.h"
#include "particle.h"
#include <cmath>

using namespace std;

namespace Manta {




//! simple test kernel and kernel with numpy array 
KERNEL(bnd=0)
void knSimpleNumpyTest(Grid<Real>& grid, PyArrayContainer npAr, Real scale)
{
	const float* p = reinterpret_cast<float*>(npAr.pData);
	grid(i,j,k) += scale * (Real)p[j*grid.getSizeX()+i]; // calc access into numpy array, no size check here!
}

//! simple test function and kernel with numpy array 
PYTHON() void simpleNumpyTest( Grid<Real>& grid, PyArrayContainer npAr, Real scale) {
	knSimpleNumpyTest(grid, npAr, scale);
}



//! extract feature vectors

KERNEL(pts)
void knExtractFeatureVel(
	const BasicParticleSystem &p, Real *fv, const IndexInt N_row, const IndexInt off_begin,
	const MACGrid &vel, const Real scale, const ParticleDataImpl<int> *ptype, const int exclude, const int window) {
	if(!p.isActive(idx) || (ptype && ((*ptype)[idx] & exclude))) return;

	const int _k = (vel.is3D()) ? -window : 0, K = (vel.is3D()) ? window : 0;
	const IndexInt D = (vel.is3D()) ? 3 : 2;
	const IndexInt off_idx = idx*N_row;

	IndexInt off_stencil = 0;
	for(int i=-window; i<=window; ++i) {
		for(int j=-window; j<=window; ++j) {
			for(int k=_k; k<=K; ++k) {
				const Vec3 off_pos(static_cast<Real>(i), static_cast<Real>(j), static_cast<Real>(k));
				const Vec3 pos_s = p[idx].pos + off_pos;

				const Vec3 vel_s = vel.getInterpolated(pos_s)*scale;
				const IndexInt off_vel = off_idx + off_begin + off_stencil*D;
				fv[off_vel + 0] = vel_s[0];
				fv[off_vel + 1] = vel_s[1];
				if(vel.is3D()) fv[off_vel + 2] = vel_s[2];

				++off_stencil;
			}
		}
	}
}

KERNEL(pts)
void knExtractFeaturePhi(
	const BasicParticleSystem &p, Real *fv, const IndexInt N_row, const IndexInt off_begin,
	const Grid<Real> &phi, const Real scale, const ParticleDataImpl<int> *ptype, const int exclude, const int window) {
	if(!p.isActive(idx) || (ptype && ((*ptype)[idx] & exclude))) return;

	const int _k = (phi.is3D()) ? -window : 0, K = (phi.is3D()) ? window : 0;
	const IndexInt off_idx = idx*N_row;

	IndexInt off_stencil = 0;
	for(int i=-window; i<=window; ++i) {
		for(int j=-window; j<=window; ++j) {
			for(int k=_k; k<=K; ++k) {
				const Vec3 off_pos(static_cast<Real>(i), static_cast<Real>(j), static_cast<Real>(k));
				const Vec3 pos_s = p[idx].pos + off_pos;

				const Real phi_s = phi.getInterpolated(pos_s)*scale;
				const IndexInt off_phi = off_idx + off_begin + off_stencil;
				fv[off_phi] = phi_s;

				++off_stencil;
			}
		}
	}
}

KERNEL(pts)
void knExtractFeatureGeo(
	const BasicParticleSystem &p, Real *fv, const IndexInt N_row, const IndexInt off_begin,
	const FlagGrid &geo, const Real scale, const ParticleDataImpl<int> *ptype, const int exclude, const int window) {
	if(!p.isActive(idx) || (ptype && ((*ptype)[idx] & exclude))) return;

	const int _k = (geo.is3D()) ? -window : 0, K = (geo.is3D()) ? window : 0;
	const IndexInt off_idx = idx*N_row;

	IndexInt off_stencil = 0;
	for(int i=-window; i<=window; ++i) {
		for(int j=-window; j<=window; ++j) {
			for(int k=_k; k<=K; ++k) {
				const Vec3 off_pos(static_cast<Real>(i), static_cast<Real>(j), static_cast<Real>(k));
				const Vec3 pos_s = p[idx].pos + off_pos;

				const Real geo_s = static_cast<Real>(geo.getAt(pos_s))*scale;
				const IndexInt off_geo = off_idx + off_begin + off_stencil;
				fv[off_geo] = geo_s;

				++off_stencil;
			}
		}
	}
}

PYTHON()
void extractFeatureVel(PyArrayContainer fv, const int N_row, const int off_begin,
		       const BasicParticleSystem &p, const MACGrid &vel,
		       const Real scale=1.0, const ParticleDataImpl<int> *ptype=NULL, const int exclude=0,
		       const int window=1) {
	knExtractFeatureVel(p, reinterpret_cast<Real*>(fv.pData), N_row, off_begin, vel, scale, ptype, exclude, window);
}
PYTHON()
void extractFeaturePhi(PyArrayContainer fv, const int N_row, const int off_begin,
		       const BasicParticleSystem &p, const Grid<Real> &phi,
		       const Real scale=1.0, const ParticleDataImpl<int> *ptype=NULL, const int exclude=0,
		       const int window=1) {
	knExtractFeaturePhi(p, reinterpret_cast<Real*>(fv.pData), N_row, off_begin, phi, scale, ptype, exclude, window);
}
PYTHON()
void extractFeatureGeo(PyArrayContainer fv, const int N_row, const int off_begin,
		       const BasicParticleSystem &p, const FlagGrid &flag,
		       const Real scale=1.0, const ParticleDataImpl<int> *ptype=NULL, const int exclude=0,
		       const int window=1) {
	knExtractFeatureGeo(p, reinterpret_cast<Real*>(fv.pData), N_row, off_begin, flag, scale, ptype, exclude, window);
}




// non-numpy related helpers

//! region detection functions

void floodFillRegion(Grid<int> &r, const FlagGrid &flags, const IndexInt idx, const int c, const int type) {
	r(idx) = c;
	if((flags(idx-flags.getStrideX()) & type) && !r[idx-flags.getStrideX()]) floodFillRegion(r, flags, idx-flags.getStrideX(), c, type);
	if((flags(idx+flags.getStrideX()) & type) && !r[idx+flags.getStrideX()]) floodFillRegion(r, flags, idx+flags.getStrideX(), c, type);
	if((flags(idx-flags.getStrideY()) & type) && !r[idx-flags.getStrideY()]) floodFillRegion(r, flags, idx-flags.getStrideY(), c, type);
	if((flags(idx+flags.getStrideY()) & type) && !r[idx+flags.getStrideY()]) floodFillRegion(r, flags, idx+flags.getStrideY(), c, type);
	if(!flags.is3D()) return;
	if((flags(idx-flags.getStrideZ()) & type) && !r[idx-flags.getStrideZ()]) floodFillRegion(r, flags, idx-flags.getStrideZ(), c, type);
	if((flags(idx+flags.getStrideZ()) & type) && !r[idx+flags.getStrideZ()]) floodFillRegion(r, flags, idx+flags.getStrideZ(), c, type);
}

PYTHON() int getRegions(Grid<int> &r, const FlagGrid &flags, const int ctype) {
	r.clear();
	int n_regions = 0;

	FOR_IDX(flags) {
		if((flags(idx) & ctype) && !r(idx)) floodFillRegion(r, flags, idx, ++n_regions, ctype);
	}
	return n_regions;
}

PYTHON() void getRegionalCounts(Grid<int> &r, const FlagGrid &flags, const int ctype) {
	const int n_regions = getRegions(r, flags, ctype);
	std::vector<int> cnt(n_regions+1, 0);
	FOR_IDX(flags) {
		if(r[idx]>0) ++(cnt[r[idx]]);
	}
	FOR_IDX(flags) {
		r[idx] = cnt[r[idx]];
	}
}

PYTHON() void extendRegion(FlagGrid &flags, const int region, const int exclude, const int depth) {
	const int I=flags.getSizeX()-1, J=flags.getSizeY()-1, K=flags.getSizeZ()-1;
	for(int i_depth=0; i_depth<depth; ++i_depth) {
		std::vector<int> update;
		FOR_IJK(flags) {
			if(flags.get(i, j, k) & exclude) continue;
			if((i>0 && (flags.get(i-1, j, k)&region)) || (i<I && (flags.get(i+1, j, k)&region)) ||
			   (j>0 && (flags.get(i, j-1, k)&region)) || (j<J && (flags.get(i, j+1, k)&region)) ||
			   (flags.is3D() && ((k>0 && (flags.get(i, j, k-1)&region)) || (k<K && (flags.get(i, j, k+1)&region)))))
				update.push_back(flags.index(i, j, k));
		}

		for(std::vector<int>::const_iterator it=update.begin(); it!=update.end(); ++it) {
			flags[*it] = region;
		}
	}
}

bool isIsolatedFluidCell(const IndexInt idx, const FlagGrid &flags) {
	if(!flags.isFluid(idx)) return false;
	if(flags.isFluid(idx-flags.getStrideX())) return false;
	if(flags.isFluid(idx+flags.getStrideX())) return false;
	if(flags.isFluid(idx-flags.getStrideY())) return false;
	if(flags.isFluid(idx+flags.getStrideY())) return false;
	if(!flags.is3D()) return true;
	if(flags.isFluid(idx-flags.getStrideZ())) return false;
	if(flags.isFluid(idx+flags.getStrideZ())) return false;
	return true;
}

KERNEL(idx)
void knMarkIsolatedFluidCell(FlagGrid &flags, const int mark) {
	if(isIsolatedFluidCell(idx, flags)) flags[idx] = mark;
}

PYTHON()
void markIsolatedFluidCell(FlagGrid &flags, const int mark) {
	knMarkIsolatedFluidCell(flags, mark);
}

KERNEL(idx)
void knMarkSmallRegions(FlagGrid &flags, const Grid<int> &rcnt, const int mark, const int exclude, const int th) {
	if(flags[idx] & exclude) return;
	if(rcnt[idx] <= th) flags[idx] = mark;
}

PYTHON()
void markSmallRegions(FlagGrid &flags, const Grid<int> &rcnt, const int mark, const int exclude, const int th=1) {
	knMarkSmallRegions(flags, rcnt, mark, exclude, th);
}

// particle helpers

KERNEL(pts)
void KnAddForcePvel(ParticleDataImpl<Vec3> &v, const Vec3 &da, const ParticleDataImpl<int> *ptype, const int exclude) {
	if(ptype && ((*ptype)[idx] & exclude)) return;
	v[idx] += da;
} 
//! add force to vec3 particle data (ie, a velocity)
PYTHON() void addForcePvel(ParticleDataImpl<Vec3> &vel, const Vec3 &a, const Real dt, const ParticleDataImpl<int> *ptype, const int exclude) {
	KnAddForcePvel(vel, a*dt, ptype, exclude);
}

KERNEL(pts) 
void KnUpdateVelocityFromDeltaPos(const BasicParticleSystem &p, ParticleDataImpl<Vec3> &v, const ParticleDataImpl<Vec3> &x_prev, const Real over_dt, const ParticleDataImpl<int> *ptype, const int exclude) {
	if(ptype && ((*ptype)[idx] & exclude)) return;
	v[idx] = (p[idx].pos - x_prev[idx])*over_dt;
} 
//! retrieve velocity from position change
PYTHON() void updateVelocityFromDeltaPos(BasicParticleSystem& parts, ParticleDataImpl<Vec3> &vel, const ParticleDataImpl<Vec3> &x_prev, const Real dt, const ParticleDataImpl<int> *ptype, const int exclude) {
	KnUpdateVelocityFromDeltaPos(parts, vel, x_prev, 1.0/dt, ptype, exclude);
}

//! simple foward Euler integration for particle system
KERNEL(pts) 
void KnStepEuler(BasicParticleSystem &p, const ParticleDataImpl<Vec3> &v, const Real dt, const ParticleDataImpl<int> *ptype, const int exclude) {
	if(ptype && ((*ptype)[idx] & exclude)) return;
	p[idx].pos += v[idx]*dt;
}
PYTHON() void eulerStep(BasicParticleSystem& parts, const ParticleDataImpl<Vec3> &vel, const ParticleDataImpl<int> *ptype, const int exclude) {
	KnStepEuler(parts, vel, parts.getParent()->getDt(), ptype, exclude);
}


KERNEL(pts) 
void KnSetType(ParticleDataImpl<int> &ptype, BasicParticleSystem &part, const int mark, const int stype, const FlagGrid &flags, const int cflag) {
	if(flags.isInBounds(part.getPos(idx), 0) && (flags.getAt(part.getPos(idx))&cflag) && (ptype[idx]&stype)) ptype[idx] = mark;
}

PYTHON() void setPartType(BasicParticleSystem &parts, ParticleDataImpl<int> &ptype, const int mark, const int stype, const FlagGrid &flags, const int cflag) {
	KnSetType(ptype, parts, mark, stype, flags, cflag);
}


} //namespace
