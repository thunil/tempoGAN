/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Fast marching and extrapolation
 *
 ******************************************************************************/

#include "fastmarch.h"
#include "levelset.h"
#include "kernel.h"
#include <algorithm>

using namespace std;

namespace Manta {
	
template<class COMP, int TDIR>
FastMarch<COMP,TDIR>::FastMarch(const FlagGrid& flags, Grid<int>& fmFlags, Grid<Real>& levelset, Real maxTime, MACGrid* velTransport )
	: mLevelset(levelset), mFlags(flags), mFmFlags(fmFlags)
{
	if (velTransport)
		mVelTransport.initMarching(velTransport, &flags);
	
	mMaxTime = maxTime * TDIR;
}

// helper for individual components to calculateDistance
template<class COMP, int TDIR> template<int C>
Real FastMarch<COMP,TDIR>::calcWeights(int& okcnt, int& invcnt, Real* v, const Vec3i& idx) {
	Real val = 0.;
	Vec3i idxPlus(idx), idxMinus(idx);
	idxPlus[C]++;
	idxMinus[C]--;
	
	mWeights[C*2] = mWeights[C*2+1] = 0.;
	if (mFmFlags(idxPlus)==FlagInited) {
		// somewhat arbitrary - choose +1 value over -1 ...
		val = mLevelset(idxPlus);
		v[okcnt] = val; okcnt++;
		mWeights[C*2] = 1.;
	} else if (mFmFlags(idxMinus)==FlagInited) {
		val = mLevelset(idxMinus);
		v[okcnt] = val; okcnt++;
		mWeights[C*2+1] = 1.;
	} 
	else {
		invcnt++;
	}
	return val;
}

template<class COMP, int TDIR>
inline Real FastMarch<COMP,TDIR>::calculateDistance(const Vec3i& idx) {
	//int invflag = 0;
	int invcnt = 0;
	Real v[3];
	int okcnt = 0;
	
	Real aVal = calcWeights<0>(okcnt, invcnt, v, idx);
	Real bVal = calcWeights<1>(okcnt, invcnt, v, idx);
	Real cVal = 0.;
	if (mLevelset.is3D())   cVal = calcWeights<2>(okcnt, invcnt, v, idx);
	else					{ invcnt++; mWeights[4] = mWeights[5] = 0.; }

	Real ret = InvalidTime();
	switch(invcnt) {
	case 0: {
		// take all values
		const Real ca=v[0], cb=v[1], cc=v[2];
		const Real csqrt = max(0. , 
				-2.*(ca*ca+cb*cb- cb*cc + cc*cc - ca*(cb+cc)) + 3 );
		// clamp to make sure the sqrt is valid
		ret = 0.333333*( ca+cb+cc+ TDIR*sqrt(csqrt) );

		// weights needed for transport (transpTouch)
		mWeights[0] *= fabs(ret-ca);
		mWeights[1] *= fabs(ret-ca);
		mWeights[2] *= fabs(ret-cb);
		mWeights[3] *= fabs(ret-cb);
		mWeights[4] *= fabs(ret-cc);
		mWeights[5] *= fabs(ret-cc);

		Real norm = 0.0; // try to force normalization
		for(int i=0;i<6;i++) { 
			norm += mWeights[i]; 
		}   
		norm = 1.0/norm;
		for(int i=0;i<6;i++) { mWeights[i] *= norm; } 

		} break; 
	case 1: {
		// take just the 2 ok values
		// t=0.5*( a+b+ (2*g*g-(b-a)*(b-a))^0.5) 
		const Real csqrt = max(0. , 2.-(v[1]-v[0])*(v[1]-v[0]) );
		// clamp to make sure the sqrt is valid
		ret = 0.5*( v[0]+v[1]+ TDIR*sqrt(csqrt) );

		// weights needed for transport (transpTouch)
		mWeights[0] *= fabs(ret-aVal);
		mWeights[1] *= fabs(ret-aVal);
		mWeights[2] *= fabs(ret-bVal);
		mWeights[3] *= fabs(ret-bVal);
		mWeights[4] *= fabs(ret-cVal);
		mWeights[5] *= fabs(ret-cVal);

		Real norm = 0.0; // try to force normalization
		for(int i=0;i<6;i++) { 
			norm += mWeights[i]; 
		}   
		norm = 1.0/norm;
		for(int i=0;i<6;i++) { mWeights[i] *= norm; } 
		// */

		} break; 
	case 2: {
		// just use the one remaining value
		ret = v[0]+ (Real)(TDIR) ; // direction = +- 1
		} break; 
	default:
		errMsg("FastMarch :: Invalid invcnt");
		break;
	}
	return ret;
}

template<class COMP, int TDIR>
void FastMarch<COMP,TDIR>::addToList(const Vec3i& p, const Vec3i& src) {
	if (!mLevelset.isInBounds(p,1)) return;
	const IndexInt idx = mLevelset.index(p);
	
	// already known value, value alreay set to valid value? skip cell...
	if(mFmFlags[idx] == FlagInited) return;

	// discard by source time now , TODO do instead before calling all addtolists?
	Real srct = mLevelset(src);
	if(COMP::compare(srct, mMaxTime)) return;

	Real ttime = calculateDistance(p);
	
	// remove old entry if larger
	bool found=false;

	Real oldt = mLevelset[idx];
	if (mFmFlags[idx] == FlagIsOnHeap) {
		found = true;
		// is old time better?
		if(COMP::compare(ttime,oldt)) return;        
	}

	// update field
	mFmFlags[idx] = FlagIsOnHeap;
	mLevelset[idx] = ttime;
	// debug info std::cout<<"set "<< idx <<","<< ttime <<"\n";
	
	if (mVelTransport.isInitialized())
		mVelTransport.transpTouch(p.x, p.y, p.z, mWeights, ttime);

	// the following adds entries to the heap of active cells
	// current: (!found) , previous: always add, might lead to duplicate
	//     entries, but the earlier will be handled earlier, the second one will skip to the FlagInited check above
	if(!found) 
	{
		// add list entry with source value
		COMP entry;
		entry.p    = p;
		entry.time = mLevelset[idx];

		mHeap.push( entry );
		// debug info std::cout<<"push "<< entry.p <<","<< entry.time <<"\n";
	}	

}

//! Enforce delta_phi = 0 on boundaries
KERNEL(single)
void SetLevelsetBoundaries (Grid<Real>& phi) {
	if (i==0)      phi(i,j,k) = phi(1,j,k);
	if (i==maxX-1) phi(i,j,k) = phi(i-1,j,k);

	if (j==0)      phi(i,j,k) = phi(i,1,k);
	if (j==maxY-1) phi(i,j,k) = phi(i,j-1,k);

	if(phi.is3D()) {
		if (k==0)      phi(i,j,k) = phi(i,j,1);
		if (k==maxZ-1) phi(i,j,k) = phi(i,j,k-1);
	}
}

/*****************************************************************************/
//! Walk...
template<class COMP, int TDIR>
void FastMarch<COMP,TDIR>::performMarching() {
	mReheapVal = 0.0;
	while(mHeap.size() > 0) {
		
		const COMP& ce = mHeap.top();
		Vec3i p = ce.p; 
		mFmFlags(p) = FlagInited;
		mHeap.pop();
		// debug info std::cout<<"pop "<< ce.p <<","<< ce.time <<"\n";

		addToList(Vec3i(p.x-1,p.y,p.z), p);
		addToList(Vec3i(p.x+1,p.y,p.z), p);
		addToList(Vec3i(p.x,p.y-1,p.z), p);
		addToList(Vec3i(p.x,p.y+1,p.z), p);
		if(mLevelset.is3D()) {
			addToList(Vec3i(p.x,p.y,p.z-1), p);
			addToList(Vec3i(p.x,p.y,p.z+1), p);        
		}
	}
	
	// set boundary for plain array
	SetLevelsetBoundaries setls(mLevelset); 
	setls.getArg0(); // get rid of compiler warning...
}

// explicit instantiation
template class FastMarch<FmHeapEntryIn, -1>;
template class FastMarch<FmHeapEntryOut, +1>;


/*****************************************************************************/
// simpler extrapolation functions (primarily for FLIP)

KERNEL(bnd=1)
void knExtrapolateMACSimple (MACGrid& vel, int distance , Grid<int>& tmp , const int d , const int c ) 
{
	static const Vec3i nb[6] = { 
		Vec3i(1 ,0,0), Vec3i(-1,0,0),
		Vec3i(0,1 ,0), Vec3i(0,-1,0),
		Vec3i(0,0,1 ), Vec3i(0,0,-1) };
	const int dim = (vel.is3D() ? 3:2);

	if (tmp(i,j,k) != 0) return;

	// copy from initialized neighbors
	Vec3i p(i,j,k);
	int nbs = 0;
	Real avgVel = 0.;
	for (int n=0; n<2*dim; ++n) {
		if (tmp(p+nb[n]) == d) {
			//vel(p)[c] = (c+1.)*0.1;
			avgVel += vel(p+nb[n])[c];
			nbs++;
		}
	}

	if(nbs>0) {
		tmp(p)    = d+1;
		vel(p)[c] = avgVel / nbs;
	}
}
//! copy velocity into domain side, note - don't read & write same grid, hence velTmp copy
KERNEL(bnd=0)
void knExtrapolateIntoBnd (FlagGrid& flags, MACGrid& vel, const MACGrid& velTmp)
{
	int c=0;
	Vec3 v(0,0,0);
	if( i==0 ) { 
		v = velTmp(i+1,j,k);
		if(v[0] < 0.) v[0] = 0.;
		c++;
	}
	else if( i==(flags.getSizeX()-1) ) { 
		v = velTmp(i-1,j,k);
		if(v[0] > 0.) v[0] = 0.;
		c++;
	}
	if( j==0 ) { 
		v = velTmp(i,j+1,k);
		if(v[1] < 0.) v[1] = 0.;
		c++;
	}
	else if( j==(flags.getSizeY()-1) ) { 
		v = velTmp(i,j-1,k);
		if(v[1] > 0.) v[1] = 0.;
		c++;
	}
	if(flags.is3D()) {
	if( k==0 ) { 
		v = velTmp(i,j,k+1);
		if(v[2] < 0.) v[2] = 0.;
		c++;
	}
	else if( k==(flags.getSizeZ()-1) ) { 
		v = velTmp(i,j,k-1);
		if(v[2] > 0.) v[2] = 0.;
		c++;
	} }
	if(c>0) {
		vel(i,j,k) = v/(Real)c;
	}
}

// todo - use getGradient instead?
inline Vec3 getNormal(const Grid<Real>& data, int i, int j, int k) {
	if (i > data.getSizeX()-2) i= data.getSizeX()-2;
	if (i < 1) i = 1;
	if (j > data.getSizeY()-2) j= data.getSizeY()-2;
	if (j < 1) j = 1;

	int kd = 1;
	if(data.is3D()) {
	if (k > data.getSizeZ()-2) k= data.getSizeZ()-2;
	if (k < 1) k = 1; 
	} else { kd=0; }

	return Vec3( data(i+1,j  ,k   ) - data(i-1,j  ,k   ) ,
				 data(i  ,j+1,k   ) - data(i  ,j-1,k   ) ,
				 data(i  ,j  ,k+kd) - data(i  ,j  ,k-kd) );
}
KERNEL(bnd=1)
void knUnprojectNormalComp (FlagGrid& flags, MACGrid& vel, Grid<Real>& phi, Real maxDist)
{
	// apply inside, within range near obstacle surface
	if(phi(i,j,k)>0. || phi(i,j,k)<-maxDist) return;

	Vec3 n = getNormal(phi, i,j,k);
	Vec3 v = vel(i,j,k);
	if(dot(n,v) < 0.) { 
		normalize(n);
		Real l = dot(n,v);
		vel(i,j,k) -= n*l;
	}
}
// a simple extrapolation step , used for cases where there's no levelset
// (note, less accurate than fast marching extrapolation.)
// into obstacle is a special mode for second order obstable boundaries (extrapolating
// only fluid velocities, not those at obstacles)
PYTHON() void extrapolateMACSimple (FlagGrid& flags, MACGrid& vel, int distance = 4, 
		LevelsetGrid* phiObs=NULL , bool intoObs = false ) 
{
	Grid<int> tmp( flags.getParent() );
	int dim = (flags.is3D() ? 3:2);

	for(int c=0; c<dim; ++c) {
		Vec3i dir = 0;
		dir[c] = 1;
		tmp.clear();

		// remove all fluid cells (not touching obstacles)
		FOR_IJK_BND(flags,1) {
			Vec3i p(i,j,k);
			bool mark = false;
			if(!intoObs) {
				if( flags.isFluid(p) || flags.isFluid(p-dir) ) mark = true;
			} else {
				if( (flags.isFluid(p) || flags.isFluid(p-dir) ) && 
					(!flags.isObstacle(p)) && (!flags.isObstacle(p-dir)) ) mark = true;
			}

			if(mark) tmp(p) = 1;
		}
		
		// extrapolate for distance
		for(int d=1; d<1+distance; ++d) {
			knExtrapolateMACSimple(vel, distance, tmp, d, c);
		} // d
	}

	if(phiObs) {
		knUnprojectNormalComp( flags, vel, *phiObs, distance );
	}

	// copy tangential values into sides of domain
	MACGrid velTmp( flags.getParent() ); velTmp.copyFrom(vel);
	knExtrapolateIntoBnd(flags, vel, velTmp);
}

KERNEL(bnd=1)
void knExtrapolateMACFromWeight ( MACGrid& vel, Grid<Vec3>& weight, int distance , const int d, const int c ) 
{
	static const Vec3i nb[6] = { 
		Vec3i(1 ,0,0), Vec3i(-1,0,0),
		Vec3i(0,1 ,0), Vec3i(0,-1,0),
		Vec3i(0,0,1 ), Vec3i(0,0,-1) };
	const int dim = (vel.is3D() ? 3:2);

	if (weight(i,j,k)[c] != 0) return;

	// copy from initialized neighbors
	Vec3i p(i,j,k);
	int nbs = 0;
	Real avgVel = 0.;
	for (int n=0; n<2*dim; ++n) {
		if (weight(p+nb[n])[c] == d) {
			avgVel += vel(p+nb[n])[c];
			nbs++;
		}
	}

	if(nbs>0) {
		weight(p)[c]    = d+1;
		vel(p)[c] = avgVel / nbs;
	}
}

// same as extrapolateMACSimple, but uses weight vec3 grid instead of flags to check
// for valid values (to be used in combination with mapPartsToMAC)
// note - the weight grid values are destroyed! the function is necessary due to discrepancies
// between velocity mapping on surface-levelset / fluid-flag creation. With this
// extrapolation we make sure the fluid region is covered by initial velocities
PYTHON() void extrapolateMACFromWeight ( MACGrid& vel, Grid<Vec3>& weight, int distance = 2) 
{
	const int dim = (vel.is3D() ? 3:2);

	for(int c=0; c<dim; ++c) {
		Vec3i dir = 0;
		dir[c] = 1;

		// reset weight values to 0 (uninitialized), and 1 (initialized inner values)
		FOR_IJK_BND(vel,1) {
			Vec3i p(i,j,k);
			if(weight(p)[c]>0.) weight(p)[c] = 1.0;
		}
		
		// extrapolate for distance
		for(int d=1; d<1+distance; ++d) {
			knExtrapolateMACFromWeight(vel, weight, distance, d, c);
		} // d

	}
}

// simple extrapolation functions for levelsets

static const Vec3i nb[6] = { 
	Vec3i(1 ,0,0), Vec3i(-1,0,0),
	Vec3i(0,1 ,0), Vec3i(0,-1,0),
	Vec3i(0,0,1 ), Vec3i(0,0,-1) };

KERNEL(bnd=1) template<class S>
void knExtrapolateLsSimple (Grid<S>& val, int distance , Grid<int>& tmp , const int d , S direction )
{
	const int dim = (val.is3D() ? 3:2); 
	if (tmp(i,j,k) != 0) return;

	// copy from initialized neighbors
	Vec3i p(i,j,k);
	int   nbs = 0;
	S     avg(0.);
	for (int n=0; n<2*dim; ++n) {
		if (tmp(p+nb[n]) == d) {
			avg += val(p+nb[n]);
			nbs++;
		}
	}

	if(nbs>0) {
		tmp(p) = d+1;
		val(p) = avg / nbs + direction;
	} 
}


KERNEL(bnd=1) template<class S>
void knSetRemaining (Grid<S>& phi, Grid<int>& tmp, S distance )
{
	if (tmp(i,j,k) != 0) return;
	phi(i,j,k) = distance;
}

PYTHON() void extrapolateLsSimple (Grid<Real>& phi, int distance = 4, bool inside=false )
{
	Grid<int> tmp( phi.getParent() );
	tmp.clear();
	const int dim = (phi.is3D() ? 3:2);

	// by default, march outside
	Real direction = 1.;
	if(!inside) { 
		// mark all inside
		FOR_IJK_BND(phi,1) {
			if ( phi(i,j,k) < 0. ) { tmp(i,j,k) = 1; }
		} 
	} else {
		direction = -1.;
		FOR_IJK_BND(phi,1) {
			if ( phi(i,j,k) > 0. ) { tmp(i,j,k) = 1; }
		} 
	}
	// + first layer around
	FOR_IJK_BND(phi,1) {
		Vec3i p(i,j,k);
		if ( tmp(p) ) continue;
		for (int n=0; n<2*dim; ++n) {
			if (tmp(p+nb[n])==1) {
				tmp(i,j,k) = 2; n=2*dim;
			}
		}
	} 

	// extrapolate for distance
	for(int d=2; d<1+distance; ++d) {
		knExtrapolateLsSimple<Real>(phi, distance, tmp, d, direction );
	} 

	// set all remaining cells to max
	knSetRemaining<Real>(phi, tmp, Real(direction * (distance+2)) );
}

// extrapolate centered vec3 values from marked fluid cells
PYTHON() void extrapolateVec3Simple (Grid<Vec3>& vel, Grid<Real>& phi, int distance = 4, bool inside=false)
{
	Grid<int> tmp( vel.getParent() );
	tmp.clear();
	const int dim = (vel.is3D() ? 3:2);

 	// mark initial cells, by default, march outside
 	if(!inside) {
 		// mark all inside
 		FOR_IJK_BND(phi,1) {
 			if ( phi(i,j,k) < 0. ) { tmp(i,j,k) = 1; }
 		}
 	} else {
 		FOR_IJK_BND(phi,1) {
 			if ( phi(i,j,k) > 0. ) { tmp(i,j,k) = 1; }
 		}
  	}
	// + first layer next to initial cells
	FOR_IJK_BND(vel,1) {
		Vec3i p(i,j,k);
		if ( tmp(p) ) continue;
		for (int n=0; n<2*dim; ++n) {
			if (tmp(p+nb[n])==1) {
				tmp(i,j,k) = 2; n=2*dim;
			}
		}
	} 

	for(int d=2; d<1+distance; ++d) {
		knExtrapolateLsSimple<Vec3>(vel, distance, tmp, d, Vec3(0.) );
	} 
	knSetRemaining<Vec3>(vel, tmp, Vec3(0.) );
}





} // namespace
