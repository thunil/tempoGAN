/******************************************************************************
 *
 * MantaFlow fluid solver framework 
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * FLIP (fluid implicit particles)
 * for use with particle data fields
 *
 ******************************************************************************/

#include "particle.h"
#include "grid.h"
#include "commonkernels.h"
#include "randomstream.h"
#include "levelset.h"
#include "shapes.h"

using namespace std;
namespace Manta {



// init

//! note - this is a simplified version , sampleLevelsetWithParticles has more functionality
PYTHON() void sampleFlagsWithParticles(const FlagGrid& flags, BasicParticleSystem& parts,
				       const int discretization, const Real randomness)
{
	const bool is3D = flags.is3D();
	const Real jlen = randomness / discretization;
	const Vec3 disp (1.0 / discretization, 1.0 / discretization, 1.0/discretization);
	RandomStream mRand(9832);
 
	FOR_IJK_BND(flags, 0) {
		if ( flags.isObstacle(i,j,k) ) continue;
		if ( flags.isFluid(i,j,k) ) {
			const Vec3 pos (i,j,k);
			for (int dk=0; dk<(is3D ? discretization : 1); dk++)
			for (int dj=0; dj<discretization; dj++)
			for (int di=0; di<discretization; di++) {
				Vec3 subpos = pos + disp * Vec3(0.5+di, 0.5+dj, 0.5+dk);
				subpos += jlen * (Vec3(1,1,1) - 2.0 * mRand.getVec3());
				if(!is3D) subpos[2] = 0.5; 
				parts.addBuffered(subpos);
			}
		}
	}
	parts.insertBufferedParticles();
}

//! sample a level set with particles, use reset to clear the particle buffer,
//! and skipEmpty for a continuous inflow (in the latter case, only empty cells will
//! be re-filled once they empty when calling sampleLevelsetWithParticles during 
//! the main loop).
PYTHON() void sampleLevelsetWithParticles(const LevelsetGrid& phi, const FlagGrid& flags, BasicParticleSystem& parts,
					  const int discretization, const Real randomness, const bool reset=false, const bool refillEmpty=false)
{
	const bool is3D = phi.is3D();
	const Real jlen = randomness / discretization;
	const Vec3 disp (1.0 / discretization, 1.0 / discretization, 1.0/discretization);
	RandomStream mRand(9832);
 
	if(reset) {
		parts.clear(); 
		parts.doCompress();
	}

	FOR_IJK_BND(phi, 0) {
		if ( flags.isObstacle(i,j,k) ) continue;
		if ( refillEmpty && flags.isFluid(i,j,k) ) continue;
		if ( phi(i,j,k) < 1.733 ) {
			const Vec3 pos (i,j,k);
			for (int dk=0; dk<(is3D ? discretization : 1); dk++)
			for (int dj=0; dj<discretization; dj++)
			for (int di=0; di<discretization; di++) {
				Vec3 subpos = pos + disp * Vec3(0.5+di, 0.5+dj, 0.5+dk);
				subpos += jlen * (Vec3(1,1,1) - 2.0 * mRand.getVec3());
				if(!is3D) subpos[2] = 0.5; 
				if( phi.getInterpolated(subpos) > 0. ) continue; 
				parts.addBuffered(subpos);
			}
		}
	}

	parts.insertBufferedParticles();
}

//! sample a shape with particles, use reset to clear the particle buffer,
//! and skipEmpty for a continuous inflow (in the latter case, only empty cells will
//! be re-filled once they empty when calling sampleShapeWithParticles during
//! the main loop).
PYTHON() void sampleShapeWithParticles(const Shape& shape, const FlagGrid& flags, BasicParticleSystem& parts,
				       const int discretization, const Real randomness, const bool reset=false, const bool refillEmpty=false,
				       const LevelsetGrid *exclude=NULL)
{
	const bool is3D = flags.is3D();
	const Real jlen = randomness / discretization;
	const Vec3 disp (1.0 / discretization, 1.0 / discretization, 1.0/discretization);
	RandomStream mRand(9832);

	if(reset) {
		parts.clear();
		parts.doCompress();
	}

	FOR_IJK_BND(flags, 0) {
		if ( flags.isObstacle(i,j,k) ) continue;
		if ( refillEmpty && flags.isFluid(i,j,k) ) continue;
		const Vec3 pos (i,j,k);
		for (int dk=0; dk<(is3D ? discretization : 1); dk++)
		for (int dj=0; dj<discretization; dj++)
		for (int di=0; di<discretization; di++) {
			Vec3 subpos = pos + disp * Vec3(0.5+di, 0.5+dj, 0.5+dk);
			subpos += jlen * (Vec3(1,1,1) - 2.0 * mRand.getVec3());
			if(!is3D) subpos[2] = 0.5;
			if(exclude && exclude->getInterpolated(subpos) <= 0.) continue;
			if(!shape.isInside(subpos)) continue;
			parts.addBuffered(subpos);
		}
	}

	parts.insertBufferedParticles();
}

//! mark fluid cells and helpers
KERNEL() void knClearFluidFlags(FlagGrid& flags, int dummy=0) {
	if (flags.isFluid(i,j,k)) {
		flags(i,j,k) = (flags(i,j,k) | FlagGrid::TypeEmpty) & ~FlagGrid::TypeFluid;
	}
}
KERNEL(bnd=1) 
void knSetNbObstacle(FlagGrid& nflags, const FlagGrid& flags, const Grid<Real>& phiObs) {
	if ( phiObs(i,j,k)>0. ) return;
	if (flags.isEmpty(i,j,k)) {
		bool set=false;
		if( (flags.isFluid(i-1,j,k)) && (phiObs(i+1,j,k)<=0.) ) set=true;
		if( (flags.isFluid(i+1,j,k)) && (phiObs(i-1,j,k)<=0.) ) set=true;
		if( (flags.isFluid(i,j-1,k)) && (phiObs(i,j+1,k)<=0.) ) set=true;
		if( (flags.isFluid(i,j+1,k)) && (phiObs(i,j-1,k)<=0.) ) set=true;
		if(flags.is3D()) {
		if( (flags.isFluid(i,j,k-1)) && (phiObs(i,j,k+1)<=0.) ) set=true;
		if( (flags.isFluid(i,j,k+1)) && (phiObs(i,j,k-1)<=0.) ) set=true;
		}
		if(set) nflags(i,j,k) = (flags(i,j,k) | FlagGrid::TypeFluid) & ~FlagGrid::TypeEmpty;
	}
}
PYTHON() void markFluidCells(const BasicParticleSystem& parts, FlagGrid& flags, const Grid<Real>* phiObs=NULL, const ParticleDataImpl<int>* ptype=NULL, const int exclude=0) {
	// remove all fluid cells
	knClearFluidFlags(flags, 0);
	
	// mark all particles in flaggrid as fluid
	for(IndexInt idx=0; idx<parts.size(); idx++) {
		if (!parts.isActive(idx) || (ptype && ((*ptype)[idx] & exclude))) continue;
		Vec3i p = toVec3i( parts.getPos(idx) );
		if (flags.isInBounds(p) && flags.isEmpty(p))
			flags(p) = (flags(p) | FlagGrid::TypeFluid) & ~FlagGrid::TypeEmpty;
	}

	// special for second order obstacle BCs, check empty cells in boundary region
	if(phiObs) {
		FlagGrid tmp(flags);
		knSetNbObstacle(tmp, flags, *phiObs);
		flags.swap(tmp);
	}
}

// for testing purposes only...
PYTHON() void testInitGridWithPos(Grid<Real>& grid) {
	FOR_IJK(grid) { grid(i,j,k) = norm( Vec3(i,j,k) ); }
}



//! helper to calculate particle radius factor to cover the diagonal of a cell in 2d/3d
inline Real calculateRadiusFactor(const Grid<Real>& grid, Real factor) {
	return (grid.is3D() ? sqrt(3.) : sqrt(2.) ) * (factor+.01); // note, a 1% safety factor is added here
} 

//! re-sample particles based on an input levelset 
// optionally skip seeding new particles in "exclude" SDF
PYTHON() void adjustNumber(BasicParticleSystem& parts, const MACGrid& vel, const FlagGrid& flags,
                int minParticles, int maxParticles, const LevelsetGrid& phi, Real radiusFactor=1. , Real narrowBand=-1. ,
                const Grid<Real>* exclude=NULL )
{
	// which levelset to use as threshold
	const Real SURFACE_LS = -1.0 * calculateRadiusFactor(phi, radiusFactor);
	Grid<int> tmp( vel.getParent() );
	std::ostringstream out;

	// count particles in cells, and delete excess particles
	for (IndexInt idx=0; idx<(int)parts.size(); idx++) {
		if (parts.isActive(idx)) {
			Vec3i p = toVec3i( parts.getPos(idx) );
			if (!tmp.isInBounds(p) ) {
				parts.kill(idx); // out of domain, remove
				continue;
			}

			Real phiv = phi.getInterpolated( parts.getPos(idx) );
			if( phiv > 0 ) { parts.kill(idx); continue; }
			if( narrowBand>0. && phiv < -narrowBand) { parts.kill(idx); continue; }

			bool atSurface = false;
			if (phiv > SURFACE_LS) atSurface = true;
			int num = tmp(p);
			
			// dont delete particles in non fluid cells here, the particles are "always right"
			if ( num > maxParticles && (!atSurface) ) {
				parts.kill(idx);
			} else {
				tmp(p) = num+1;
			}
		}
	}

	// seed new particles
	RandomStream mRand(9832);
	FOR_IJK(tmp) {
		int cnt = tmp(i,j,k);
		
		// skip cells near surface
		if (phi(i,j,k) > SURFACE_LS) continue;
		if( narrowBand>0. && phi(i,j,k) < -narrowBand ) { continue; }
		if( exclude && ( (*exclude)(i,j,k) < 0.) ) { continue; }

		if (flags.isFluid(i,j,k) && cnt < minParticles) {
			for (int m=cnt; m < minParticles; m++) { 
				Vec3 pos = Vec3(i,j,k) + mRand.getVec3();
				//Vec3 pos (i + 0.5, j + 0.5, k + 0.5); // cell center
				parts.addBuffered( pos ); 
			}
		}
	}

	parts.doCompress();
	parts.insertBufferedParticles();
}

// simple and slow helper conversion to show contents of int grids like a real grid in the ui
// (use eg to quickly display contents of the particle-index grid)
PYTHON() void debugIntToReal(const Grid<int>& source, Grid<Real>& dest, Real factor=1. )
{
	FOR_IJK( source ) { dest(i,j,k) = (Real)source(i,j,k) * factor; }
}

// build a grid that contains indices for a particle system
// the particles in a cell i,j,k are particles[index(i,j,k)] to particles[index(i+1,j,k)-1]
// (ie,  particles[index(i+1,j,k)] already belongs to cell i+1,j,k)
PYTHON() void gridParticleIndex(const BasicParticleSystem& parts, ParticleIndexSystem& indexSys,
                                 const FlagGrid& flags, Grid<int>& index, Grid<int>* counter=NULL )
{
	bool delCounter = false;
	if(!counter) { counter = new Grid<int>(  flags.getParent() ); delCounter=true; }
	else         { counter->clear(); }
	
	// count particles in cells, and delete excess particles
	index.clear();
	int inactive = 0;
	for (IndexInt idx=0; idx<(IndexInt)parts.size(); idx++) {
		if (parts.isActive(idx)) {
			// check index for validity...
			Vec3i p = toVec3i( parts.getPos(idx) );
			if (! index.isInBounds(p)) { inactive++; continue; }

			index(p)++;
		} else {
			inactive++;
		}
	}

	// note - this one might be smaller...
	indexSys.resize( parts.size()-inactive );

	// convert per cell number to continuous index
	IndexInt idx=0;
	FOR_IJK( index ) {
		int num = index(i,j,k);
		index(i,j,k) = idx;
		idx += num;
	}

	// add particles to indexed array, we still need a per cell particle counter
	for (IndexInt idx=0; idx<(IndexInt)parts.size(); idx++) {
		if (!parts.isActive(idx)) continue;
		Vec3i p = toVec3i( parts.getPos(idx) );
		if (! index.isInBounds(p)) { continue; }

		// initialize position and index into original array
		//indexSys[ index(p)+(*counter)(p) ].pos        = parts[idx].pos;
		indexSys[ index(p)+(*counter)(p) ].sourceIndex = idx;
		(*counter)(p)++;
	}

	if(delCounter) delete counter;
}

KERNEL()
void ComputeUnionLevelsetPindex(const Grid<int>& index, const BasicParticleSystem& parts, const ParticleIndexSystem& indexSys,
				LevelsetGrid& phi, const Real radius,
				const ParticleDataImpl<int> *ptype, const int exclude)
{
	const Vec3 gridPos = Vec3(i,j,k) + Vec3(0.5); // shifted by half cell
	Real phiv = radius * 1.0;  // outside

	int r  = int(radius) + 1;
	int rZ = phi.is3D() ? r : 0;
	for(int zj=k-rZ; zj<=k+rZ; zj++) 
	for(int yj=j-r ; yj<=j+r ; yj++) 
	for(int xj=i-r ; xj<=i+r ; xj++) {
		if (!phi.isInBounds(Vec3i(xj,yj,zj))) continue;

		// note, for the particle indices in indexSys the access is periodic (ie, dont skip for eg inBounds(sx,10,10)
		IndexInt isysIdxS = index.index(xj,yj,zj);
		IndexInt pStart = index(isysIdxS), pEnd=0;
		if(phi.isInBounds(isysIdxS+1)) pEnd = index(isysIdxS+1);
		else                           pEnd = indexSys.size();

		// now loop over particles in cell
		for(IndexInt p=pStart; p<pEnd; ++p) {
			const int psrc = indexSys[p].sourceIndex;
			if(ptype && ((*ptype)[psrc] & exclude)) continue;
			const Vec3 pos = parts[psrc].pos;
			phiv = std::min( phiv , fabs( norm(gridPos-pos) )-radius );
		}
	}
	phi(i,j,k) = phiv;
}
 
PYTHON() void unionParticleLevelset(const BasicParticleSystem& parts, const ParticleIndexSystem& indexSys,
				    const FlagGrid& flags, const Grid<int>& index, LevelsetGrid& phi, const Real radiusFactor=1.,
				    const ParticleDataImpl<int> *ptype=NULL, const int exclude=0)
{
	// use half a cell diagonal as base radius
	const Real radius = 0.5 * calculateRadiusFactor(phi, radiusFactor);
	// no reset of phi necessary here 
	ComputeUnionLevelsetPindex(index, parts, indexSys, phi, radius, ptype, exclude);

	phi.setBound(0.5, 0);
}


KERNEL()
void ComputeAveragedLevelsetWeight(const BasicParticleSystem& parts,
				   const Grid<int>& index, const ParticleIndexSystem& indexSys,
				   LevelsetGrid& phi, const Real radius,
				   const ParticleDataImpl<int>* ptype, const int exclude)
{
	const Vec3 gridPos = Vec3(i,j,k) + Vec3(0.5); // shifted by half cell
	Real phiv = radius * 1.0; // outside 

	// loop over neighborhood, similar to ComputeUnionLevelsetPindex
	const Real sradiusInv = 1. / (4. * radius * radius) ;
	int   r = int(1. * radius) + 1;
	int   rZ = phi.is3D() ? r : 0;
	// accumulators
	Real  wacc = 0.;
	Vec3  pacc = Vec3(0.);
	Real  racc = 0.;

	for(int zj=k-rZ; zj<=k+rZ; zj++) 
	for(int yj=j-r ; yj<=j+r ; yj++) 
	for(int xj=i-r ; xj<=i+r ; xj++) {
		if (! phi.isInBounds(Vec3i(xj,yj,zj)) ) continue;

		IndexInt isysIdxS = index.index(xj,yj,zj);
		IndexInt pStart = index(isysIdxS), pEnd=0;
		if(phi.isInBounds(isysIdxS+1)) pEnd = index(isysIdxS+1);
		else                           pEnd = indexSys.size();
		for(IndexInt p=pStart; p<pEnd; ++p) {
			IndexInt   psrc = indexSys[p].sourceIndex;
			if(ptype && ((*ptype)[psrc] & exclude)) continue;

			Vec3  pos  = parts[psrc].pos; 
			Real  s    = normSquare(gridPos-pos) * sradiusInv;
			//Real  w = std::max(0., cubed(1.-s) );
			Real  w = std::max(0., (1.-s)); // a bit smoother
			wacc += w;
			racc += radius * w;
			pacc += pos    * w;
		} 
	}

	if(wacc > VECTOR_EPSILON) {
		racc /= wacc;
		pacc /= wacc;
		phiv = fabs( norm(gridPos-pacc) )-racc;
	}
	phi(i,j,k) = phiv;
}

template<class T> T smoothingValue(const Grid<T> val, int i, int j, int k, T center) {
	return val(i,j,k);
}

// smoothing, and  
KERNEL(bnd=1) template<class T> 
void knSmoothGrid(const Grid<T>& me, Grid<T>& tmp, Real factor) {
	T val = me(i,j,k) + 
			me(i+1,j,k) + me(i-1,j,k) + 
			me(i,j+1,k) + me(i,j-1,k) ;
	if(me.is3D()) {
		val += me(i,j,k+1) + me(i,j,k-1);
	}
	tmp(i,j,k) = val * factor;
}

KERNEL(bnd=1) template<class T> 
void knSmoothGridNeg(const Grid<T>& me, Grid<T>& tmp, Real factor) {
	T val = me(i,j,k) + 
			me(i+1,j,k) + me(i-1,j,k) + 
			me(i,j+1,k) + me(i,j-1,k) ;
	if(me.is3D()) {
		val += me(i,j,k+1) + me(i,j,k-1);
	}
	val *= factor;
	if(val<tmp(i,j,k)) tmp(i,j,k) = val;
	else               tmp(i,j,k) = me(i,j,k);
}

 
PYTHON() void averagedParticleLevelset(const BasicParticleSystem& parts, const ParticleIndexSystem& indexSys,
				       const FlagGrid& flags, const Grid<int>& index, LevelsetGrid& phi, const Real radiusFactor=1.,
				       const int smoothen=1, const int smoothenNeg=1,
				       const ParticleDataImpl<int>* ptype=NULL, const int exclude=0)
{
	// use half a cell diagonal as base radius
	const Real radius = 0.5 * calculateRadiusFactor(phi, radiusFactor); 
	ComputeAveragedLevelsetWeight(parts, index, indexSys, phi, radius, ptype, exclude);

	// post-process level-set
	for(int i=0; i<std::max(smoothen,smoothenNeg); ++i) {
		LevelsetGrid tmp(flags.getParent());
		if(i<smoothen) { 
			knSmoothGrid    <Real> (phi,tmp, 1./(phi.is3D() ? 7. : 5.) );
			phi.swap(tmp);
		}
		if(i<smoothenNeg) { 
			knSmoothGridNeg <Real> (phi,tmp, 1./(phi.is3D() ? 7. : 5.) );
			phi.swap(tmp);
		}
	} 
	phi.setBound(0.5, 0);
}

KERNEL(pts)
void knPushOutofObs(BasicParticleSystem& parts, const FlagGrid& flags, const Grid<Real>& phiObs, const Real shift, const Real thresh,
		    const ParticleDataImpl<int>* ptype, const int exclude) {
	if (!parts.isActive(idx) || (ptype && ((*ptype)[idx] & exclude))) return;
	Vec3i p = toVec3i( parts.getPos(idx) );

	if (!flags.isInBounds(p)) return;
	Real v = phiObs.getInterpolated(parts.getPos(idx));
	if(v < thresh) {
		Vec3 grad = getGradient( phiObs, p.x,p.y,p.z );
		if( normalize(grad) < VECTOR_EPSILON ) return;
		parts.setPos(idx, parts.getPos(idx) + grad*(thresh - v + shift));
	}
}
//! push particles out of obstacle levelset
PYTHON() void pushOutofObs(BasicParticleSystem& parts, const FlagGrid& flags, const Grid<Real>& phiObs, const Real shift=0, const Real thresh=0,
			   const ParticleDataImpl<int>* ptype=NULL, const int exclude=0) {
	knPushOutofObs(parts, flags, phiObs, shift, thresh, ptype, exclude);
}

//******************************************************************************
// grid interpolation functions

KERNEL(idx) template<class T> 
void knSafeDivReal(Grid<T>& me, const Grid<Real>& other, Real cutoff=VECTOR_EPSILON) { 
	if(other[idx]<cutoff) {
		me[idx] = 0.;
	} else {
		T div( other[idx] );
		me[idx] = safeDivide(me[idx], div ); 
	}
}

// Set velocities on the grid from the particle system

KERNEL(pts, single)
void knMapLinearVec3ToMACGrid(const BasicParticleSystem& p, const FlagGrid& flags, const MACGrid& vel, Grid<Vec3>& tmp,
			      const ParticleDataImpl<Vec3>& pvel, const ParticleDataImpl<int>* ptype, const int exclude)
{
	unusedParameter(flags);
	if (!p.isActive(idx) || (ptype && ((*ptype)[idx] & exclude))) return;
	vel.setInterpolated( p[idx].pos, pvel[idx], &tmp[0] );
}

// optionally , this function can use an existing vec3 grid to store the weights
// this is useful in combination with the simple extrapolation function
PYTHON() void mapPartsToMAC(const FlagGrid& flags, MACGrid& vel, MACGrid& velOld,
			    const BasicParticleSystem& parts, const ParticleDataImpl<Vec3>& partVel, Grid<Vec3>* weight=NULL,
			    const ParticleDataImpl<int>* ptype=NULL, const int exclude=0)
{
	// interpol -> grid. tmpgrid for particle contribution weights
	bool freeTmp = false;
	if(!weight) {
		weight = new Grid<Vec3>(flags.getParent());
		freeTmp = true;
	} else {
		weight->clear(); // make sure we start with a zero grid!
	}
	vel.clear();
	knMapLinearVec3ToMACGrid( parts, flags, vel, *weight, partVel, ptype, exclude );

	// stomp small values in weight to zero to prevent roundoff errors
	weight->stomp(Vec3(VECTOR_EPSILON));
	vel.safeDivide(*weight);
	
	// store original state
	velOld.copyFrom( vel );
	if(freeTmp) delete weight;
}

KERNEL(pts, single) template<class T>
void knMapLinear(const BasicParticleSystem& p, const FlagGrid& flags, const Grid<T>& target, Grid<Real>& gtmp,
        const ParticleDataImpl<T>& psource )
{
	unusedParameter(flags);
	if (!p.isActive(idx)) return;
	target.setInterpolated( p[idx].pos, psource[idx], gtmp );
} 

template<class T>
void mapLinearRealHelper(const FlagGrid& flags, Grid<T>& target,
                const BasicParticleSystem& parts, const ParticleDataImpl<T>& source )
{
	Grid<Real> tmp(flags.getParent());
	target.clear();
	knMapLinear<T>( parts, flags, target, tmp, source ); 
	knSafeDivReal<T>( target, tmp );
}

PYTHON() void mapPartsToGrid    (const FlagGrid& flags, Grid<Real>& target , const BasicParticleSystem& parts , const ParticleDataImpl<Real>& source ) {
	mapLinearRealHelper<Real>(flags,target,parts,source);
}
PYTHON() void mapPartsToGridVec3(const FlagGrid& flags, Grid<Vec3>& target , const BasicParticleSystem& parts , const ParticleDataImpl<Vec3>& source ) {
	mapLinearRealHelper<Vec3>(flags,target,parts,source);
}
// integers need "max" mode, not yet implemented
//PYTHON() void mapPartsToGridInt ( FlagGrid& flags, Grid<int >& target , BasicParticleSystem& parts , ParticleDataImpl<int >& source ) {
//	mapLinearRealHelper<int >(flags,target,parts,source);
//}

KERNEL(pts) template<class T>
void knMapFromGrid(const BasicParticleSystem& p, const Grid<T>& gsrc, ParticleDataImpl<T>& target )
{
	if (!p.isActive(idx)) return;
	target[idx] = gsrc.getInterpolated( p[idx].pos );
} 
PYTHON() void mapGridToParts    (const Grid<Real>& source , const BasicParticleSystem& parts , ParticleDataImpl<Real>& target ) {
	knMapFromGrid<Real>(parts, source, target);
}
PYTHON() void mapGridToPartsVec3(const Grid<Vec3>& source , const BasicParticleSystem& parts , ParticleDataImpl<Vec3>& target ) {
	knMapFromGrid<Vec3>(parts, source, target);
}


// Get velocities from grid

KERNEL(pts) 
void knMapLinearMACGridToVec3_PIC(const BasicParticleSystem& p, const FlagGrid& flags, const MACGrid& vel, ParticleDataImpl<Vec3>& pvel,
				  const ParticleDataImpl<int>* ptype, const int exclude)
{
	if (!p.isActive(idx) || (ptype && ((*ptype)[idx] & exclude))) return;
	// pure PIC
	pvel[idx] = vel.getInterpolated( p[idx].pos );
}
PYTHON() void mapMACToParts(const FlagGrid& flags, const MACGrid& vel ,
                            const BasicParticleSystem& parts , ParticleDataImpl<Vec3>& partVel,
			    const ParticleDataImpl<int>* ptype=NULL, const int exclude=0) {
	knMapLinearMACGridToVec3_PIC( parts, flags, vel, partVel, ptype, exclude );
}

// with flip delta interpolation 
KERNEL(pts) 
void knMapLinearMACGridToVec3_FLIP(const BasicParticleSystem& p, const FlagGrid& flags, const MACGrid& vel, const MACGrid& oldVel,
				   ParticleDataImpl<Vec3>& pvel,
				   const Real flipRatio, const ParticleDataImpl<int>* ptype, const int exclude)
{
	if (!p.isActive(idx) || (ptype && ((*ptype)[idx] & exclude))) return;
	Vec3 v     =        vel.getInterpolated(p[idx].pos);
	Vec3 delta = v - oldVel.getInterpolated(p[idx].pos); 
	pvel[idx] = flipRatio * (pvel[idx] + delta) + (1.0 - flipRatio) * v;    
}

PYTHON() void flipVelocityUpdate(const FlagGrid& flags, const MACGrid& vel, const MACGrid& velOld,
				 const BasicParticleSystem& parts, ParticleDataImpl<Vec3>& partVel, const Real flipRatio,
				 const ParticleDataImpl<int>* ptype=NULL, const int exclude=0) {
	knMapLinearMACGridToVec3_FLIP( parts, flags, vel, velOld, partVel, flipRatio, ptype, exclude );
}


//******************************************************************************
// narrow band 

KERNEL()
void knCombineVels(MACGrid& vel, const Grid<Vec3>& w, MACGrid& combineVel, const LevelsetGrid* phi, Real narrowBand, Real thresh ) {
	int idx = vel.index(i,j,k);

	for(int c=0; c<3; ++c)
	{
			// Correct narrow-band FLIP
			Vec3 pos(i,j,k);
			pos[(c+1)%3] += Real(0.5);
			pos[(c+2)%3] += Real(0.5);
			Real p = phi->getInterpolated(pos);

			if (p < -narrowBand) { vel[idx][c] = 0; continue; }

			if (w[idx][c] > thresh) {
				combineVel[idx][c] = vel[idx][c];
				vel[idx][c] = -1;
			}
			else
			{
				vel[idx][c] = 0;
			}
	}
}

//! narrow band velocity combination
PYTHON() void combineGridVel( MACGrid& vel, const Grid<Vec3>& weight, MACGrid& combineVel, const LevelsetGrid* phi=NULL,
    Real narrowBand=0.0, Real thresh=0.0) {
	knCombineVels(vel, weight, combineVel, phi, narrowBand, thresh);
}

//! surface tension helper
PYTHON() void getLaplacian(Grid<Real> &laplacian, const Grid<Real> &grid) {
	LaplaceOp(laplacian, grid);
}


} // namespace

