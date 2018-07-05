/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Wavelet noise field
 *
 ******************************************************************************/
 
#ifndef _NOISEFIELD_H_
#define _NOISEFIELD_H_

#include "vectorbase.h"
#include "manta.h"
#include "grid.h"

namespace Manta {

#define NOISE_TILE_SIZE 128

// wrapper for a parametrized field of wavelet noise
PYTHON(name=NoiseField) 
class WaveletNoiseField : public PbClass {
	public:     
		PYTHON() WaveletNoiseField( FluidSolver* parent, int fixedSeed=-1 , int loadFromFile=false );
		~WaveletNoiseField() {
			if(mNoiseTile) { delete mNoiseTile; mNoiseTile=NULL; }
		};

		//! evaluate noise
		inline Real evaluate(Vec3 pos, int tile=0) const;
		//! evaluate noise as a vector
		inline Vec3 evaluateVec(Vec3 pos, int tile=0) const;
		//! evaluate curl noise
		inline Vec3 evaluateCurl(Vec3 pos) const;

		//! direct data access
		Real* data() { return mNoiseTile; }
		
		//! compute wavelet decomposition of an input grid (stores residual coefficients)
		static void computeCoefficients(Grid<Real>& input, Grid<Real>& tempIn1, Grid<Real>& tempIn2);

		// helper
		std::string toString();
		
		// texcoord position and scale
		PYTHON(name=posOffset) Vec3 mPosOffset;
		PYTHON(name=posScale)  Vec3 mPosScale;
		// value offset & scale
		PYTHON(name=valOffset) Real mValOffset;
		PYTHON(name=valScale)  Real mValScale;
		// clamp? (default 0-1)
		PYTHON(name=clamp)     bool mClamp;
		PYTHON(name=clampNeg)  Real mClampNeg;
		PYTHON(name=clampPos)  Real mClampPos;
		// animated over time
		PYTHON(name=timeAnim)  Real mTimeAnim;
		
	protected:
		// noise evaluation functions
		static inline Real WNoiseDx (const Vec3& p, Real *data);
		static inline Vec3  WNoiseVec(const Vec3& p, Real *data);
		static inline Real WNoise   (const Vec3& p, Real *data);

		// helpers for tile generation , for periodic 128 grids only
		static void downsample(Real *from, Real *to, int n, int stride);
		static void upsample  (Real *from, Real *to, int n, int stride);

		// for grids with arbitrary sizes, and neumann boundary conditions
		static void downsampleNeumann(const Real *from, Real *to, int n, int stride);
		static void upsampleNeumann   (const Real *from, Real *to, int n, int stride);

		static inline int modSlow(int x, int n) { int m = x % n; return (m<0) ? m+n : m; }
		// warning - noiseTileSize has to be 128^3!
		#define modFast128(x)  ((x) & 127)

		inline Real getTime() const { return mParent->getTime() * mParent->getDx() * mTimeAnim; }

		// pre-compute tile data for wavelet noise
		void generateTile( int loadFromFile );
				
		// animation over time
		// grid size normalization (inverse size)
		Real mGsInvX, mGsInvY, mGsInvZ;
		// random offset into tile to simulate different random seeds
		Vec3 mSeedOffset;
		
		static Real* mNoiseTile;
		// global random seed storage
		static int randomSeed;
};



// **************************************************************************
// Implementation

#define ADD_WEIGHTED(x,y,z)\
  weight = 1.0f;\
  xC = modFast128(midX + (x));\
  weight *= w[0][(x) + 1];\
  yC = modFast128(midY + (y));\
  weight *= w[1][(y) + 1];\
  zC = modFast128(midZ + (z));\
  weight *= w[2][(z) + 1];\
  result += weight * data[(zC * NOISE_TILE_SIZE + yC) * NOISE_TILE_SIZE + xC];

//////////////////////////////////////////////////////////////////////////////////////////
// derivatives of 3D noise - unrolled for performance
//////////////////////////////////////////////////////////////////////////////////////////
inline Real WaveletNoiseField::WNoiseDx(const Vec3& p, Real *data) {
	Real w[3][3], t, result = 0;

	// Evaluate quadratic B-spline basis functions
	int midX = (int)ceil(p[0] - 0.5f); 
	t        =   midX - (p[0] - 0.5f);
	w[0][0] = -t;
	w[0][2] = (1.f - t);
	w[0][1] = 2.0f * t - 1.0f;

	int midY = (int)ceil(p[1] - 0.5f); 
	t        =   midY - (p[1] - 0.5f);
	w[1][0] = t * t * 0.5f; 
	w[1][2] = (1.f - t) * (1.f - t) *0.5f; 
	w[1][1] = 1.f - w[1][0] - w[1][2];

	int midZ = (int)ceil(p[2] - 0.5f); 
	t        =   midZ - (p[2] - 0.5f);
	w[2][0] = t * t * 0.5f; 
	w[2][2] = (1.f - t) * (1.f - t) *0.5f; 
	w[2][1] = 1.f - w[2][0] - w[2][2];

	// Evaluate noise by weighting noise coefficients by basis function values
	int xC, yC, zC;
	Real weight = 1;

	ADD_WEIGHTED(-1,-1, -1); ADD_WEIGHTED( 0,-1, -1); ADD_WEIGHTED( 1,-1, -1);
	ADD_WEIGHTED(-1, 0, -1); ADD_WEIGHTED( 0, 0, -1); ADD_WEIGHTED( 1, 0, -1);
	ADD_WEIGHTED(-1, 1, -1); ADD_WEIGHTED( 0, 1, -1); ADD_WEIGHTED( 1, 1, -1);

	ADD_WEIGHTED(-1,-1, 0);  ADD_WEIGHTED( 0,-1, 0);  ADD_WEIGHTED( 1,-1, 0);
	ADD_WEIGHTED(-1, 0, 0);  ADD_WEIGHTED( 0, 0, 0);  ADD_WEIGHTED( 1, 0, 0);
	ADD_WEIGHTED(-1, 1, 0);  ADD_WEIGHTED( 0, 1, 0);  ADD_WEIGHTED( 1, 1, 0);

	ADD_WEIGHTED(-1,-1, 1);  ADD_WEIGHTED( 0,-1, 1);  ADD_WEIGHTED( 1,-1, 1);
	ADD_WEIGHTED(-1, 0, 1);  ADD_WEIGHTED( 0, 0, 1);  ADD_WEIGHTED( 1, 0, 1);
	ADD_WEIGHTED(-1, 1, 1);  ADD_WEIGHTED( 0, 1, 1);  ADD_WEIGHTED( 1, 1, 1);

	return result;
}

inline Real WaveletNoiseField::WNoise(const Vec3& p, Real *data) {
	Real w[3][3], t, result = 0;

	// Evaluate quadratic B-spline basis functions
	int midX = (int)ceilf(p[0] - 0.5f); 
	t        =   midX - (p[0] - 0.5f);
	w[0][0] = t * t * 0.5f; 
	w[0][2] = (1.f - t) * (1.f - t) *0.5f; 
	w[0][1] = 1.f - w[0][0] - w[0][2];

	int midY = (int)ceilf(p[1] - 0.5f); 
	t        =   midY - (p[1] - 0.5f);
	w[1][0] = t * t * 0.5f; 
	w[1][2] = (1.f - t) * (1.f - t) *0.5f; 
	w[1][1] = 1.f - w[1][0] - w[1][2];

	int midZ = (int)ceilf(p[2] - 0.5f); 
	t        =   midZ - (p[2] - 0.5f);
	w[2][0] = t * t * 0.5f; 
	w[2][2] = (1.f - t) * (1.f - t) *0.5f; 
	w[2][1] = 1.f - w[2][0] - w[2][2];

	// Evaluate noise by weighting noise coefficients by basis function values
	int xC, yC, zC;
	Real weight = 1;

	ADD_WEIGHTED(-1,-1, -1); ADD_WEIGHTED( 0,-1, -1); ADD_WEIGHTED( 1,-1, -1);
	ADD_WEIGHTED(-1, 0, -1); ADD_WEIGHTED( 0, 0, -1); ADD_WEIGHTED( 1, 0, -1);
	ADD_WEIGHTED(-1, 1, -1); ADD_WEIGHTED( 0, 1, -1); ADD_WEIGHTED( 1, 1, -1);

	ADD_WEIGHTED(-1,-1, 0);  ADD_WEIGHTED( 0,-1, 0);  ADD_WEIGHTED( 1,-1, 0);
	ADD_WEIGHTED(-1, 0, 0);  ADD_WEIGHTED( 0, 0, 0);  ADD_WEIGHTED( 1, 0, 0);
	ADD_WEIGHTED(-1, 1, 0);  ADD_WEIGHTED( 0, 1, 0);  ADD_WEIGHTED( 1, 1, 0);

	ADD_WEIGHTED(-1,-1, 1);  ADD_WEIGHTED( 0,-1, 1);  ADD_WEIGHTED( 1,-1, 1);
	ADD_WEIGHTED(-1, 0, 1);  ADD_WEIGHTED( 0, 0, 1);  ADD_WEIGHTED( 1, 0, 1);
	ADD_WEIGHTED(-1, 1, 1);  ADD_WEIGHTED( 0, 1, 1);  ADD_WEIGHTED( 1, 1, 1);

	return result;
}



#define ADD_WEIGHTEDX(x,y,z)\
  weight = dw[0][(x) + 1] * w[1][(y) + 1] * w[2][(z) + 1];\
  result += weight * neighbors[x + 1][y + 1][z + 1];

#define ADD_WEIGHTEDY(x,y,z)\
  weight = w[0][(x) + 1] * dw[1][(y) + 1] * w[2][(z) + 1];\
  result += weight * neighbors[x + 1][y + 1][z + 1];

#define ADD_WEIGHTEDZ(x,y,z)\
  weight = w[0][(x) + 1] * w[1][(y) + 1] * dw[2][(z) + 1];\
  result += weight * neighbors[x + 1][y + 1][z + 1];

//////////////////////////////////////////////////////////////////////////////////////////
// compute all derivatives in at once
//////////////////////////////////////////////////////////////////////////////////////////
inline Vec3 WaveletNoiseField::WNoiseVec(const Vec3& p, Real *data)
{
	Vec3 final(0.);
	Real w[3][3];
	Real dw[3][3];
	Real result = 0;
	int xC, yC, zC;
	Real weight;

	int midX = (int)ceil(p[0] - 0.5f); 
	int midY = (int)ceil(p[1] - 0.5f); 
	int midZ = (int)ceil(p[2] - 0.5f);

	Real t0 =   midX - (p[0] - 0.5f);
	Real t1 =   midY - (p[1] - 0.5f);
	Real t2 =   midZ - (p[2] - 0.5f);

	// precache all the neighbors for fast access
	Real neighbors[3][3][3];
	for (int z = -1; z <=1; z++)
		for (int y = -1; y <= 1; y++)
			for (int x = -1; x <= 1; x++)
			{
				xC = modFast128(midX + (x));
				yC = modFast128(midY + (y));
				zC = modFast128(midZ + (z));
				neighbors[x + 1][y + 1][z + 1] = data[zC * NOISE_TILE_SIZE * NOISE_TILE_SIZE + yC * NOISE_TILE_SIZE + xC];
			}

	///////////////////////////////////////////////////////////////////////////////////////
	// evaluate splines
	///////////////////////////////////////////////////////////////////////////////////////
	dw[0][0] = -t0;
	dw[0][2] = (1.f - t0);
	dw[0][1] = 2.0f * t0 - 1.0f;

	dw[1][0] = -t1;
	dw[1][2] = (1.0f - t1);
	dw[1][1] = 2.0f * t1 - 1.0f;

	dw[2][0] = -t2;
	dw[2][2] = (1.0f - t2);
	dw[2][1] = 2.0f * t2 - 1.0f;

	w[0][0] = t0 * t0 * 0.5f; 
	w[0][2] = (1.f - t0) * (1.f - t0) *0.5f; 
	w[0][1] = 1.f - w[0][0] - w[0][2];

	w[1][0] = t1 * t1 * 0.5f; 
	w[1][2] = (1.f - t1) * (1.f - t1) *0.5f; 
	w[1][1] = 1.f - w[1][0] - w[1][2];

	w[2][0] = t2 * t2 * 0.5f; 
	w[2][2] = (1.f - t2) * (1.f - t2) *0.5f;
	w[2][1] = 1.f - w[2][0] - w[2][2];

	///////////////////////////////////////////////////////////////////////////////////////
	// x derivative
	///////////////////////////////////////////////////////////////////////////////////////
	result = 0.0f;
	ADD_WEIGHTEDX(-1,-1, -1); ADD_WEIGHTEDX( 0,-1, -1); ADD_WEIGHTEDX( 1,-1, -1);
	ADD_WEIGHTEDX(-1, 0, -1); ADD_WEIGHTEDX( 0, 0, -1); ADD_WEIGHTEDX( 1, 0, -1);
	ADD_WEIGHTEDX(-1, 1, -1); ADD_WEIGHTEDX( 0, 1, -1); ADD_WEIGHTEDX( 1, 1, -1);

	ADD_WEIGHTEDX(-1,-1, 0);  ADD_WEIGHTEDX( 0,-1, 0);  ADD_WEIGHTEDX( 1,-1, 0);
	ADD_WEIGHTEDX(-1, 0, 0);  ADD_WEIGHTEDX( 0, 0, 0);  ADD_WEIGHTEDX( 1, 0, 0);
	ADD_WEIGHTEDX(-1, 1, 0);  ADD_WEIGHTEDX( 0, 1, 0);  ADD_WEIGHTEDX( 1, 1, 0);

	ADD_WEIGHTEDX(-1,-1, 1);  ADD_WEIGHTEDX( 0,-1, 1);  ADD_WEIGHTEDX( 1,-1, 1);
	ADD_WEIGHTEDX(-1, 0, 1);  ADD_WEIGHTEDX( 0, 0, 1);  ADD_WEIGHTEDX( 1, 0, 1);
	ADD_WEIGHTEDX(-1, 1, 1);  ADD_WEIGHTEDX( 0, 1, 1);  ADD_WEIGHTEDX( 1, 1, 1);
	final[0] = result;

	///////////////////////////////////////////////////////////////////////////////////////
	// y derivative
	///////////////////////////////////////////////////////////////////////////////////////
	result = 0.0f;
	ADD_WEIGHTEDY(-1,-1, -1); ADD_WEIGHTEDY( 0,-1, -1); ADD_WEIGHTEDY( 1,-1, -1);
	ADD_WEIGHTEDY(-1, 0, -1); ADD_WEIGHTEDY( 0, 0, -1); ADD_WEIGHTEDY( 1, 0, -1);
	ADD_WEIGHTEDY(-1, 1, -1); ADD_WEIGHTEDY( 0, 1, -1); ADD_WEIGHTEDY( 1, 1, -1);

	ADD_WEIGHTEDY(-1,-1, 0);  ADD_WEIGHTEDY( 0,-1, 0);  ADD_WEIGHTEDY( 1,-1, 0);
	ADD_WEIGHTEDY(-1, 0, 0);  ADD_WEIGHTEDY( 0, 0, 0);  ADD_WEIGHTEDY( 1, 0, 0);
	ADD_WEIGHTEDY(-1, 1, 0);  ADD_WEIGHTEDY( 0, 1, 0);  ADD_WEIGHTEDY( 1, 1, 0);

	ADD_WEIGHTEDY(-1,-1, 1);  ADD_WEIGHTEDY( 0,-1, 1);  ADD_WEIGHTEDY( 1,-1, 1);
	ADD_WEIGHTEDY(-1, 0, 1);  ADD_WEIGHTEDY( 0, 0, 1);  ADD_WEIGHTEDY( 1, 0, 1);
	ADD_WEIGHTEDY(-1, 1, 1);  ADD_WEIGHTEDY( 0, 1, 1);  ADD_WEIGHTEDY( 1, 1, 1);
	final[1] = result;

	///////////////////////////////////////////////////////////////////////////////////////
	// z derivative
	///////////////////////////////////////////////////////////////////////////////////////
	result = 0.0f;
	ADD_WEIGHTEDZ(-1,-1, -1); ADD_WEIGHTEDZ( 0,-1, -1); ADD_WEIGHTEDZ( 1,-1, -1);
	ADD_WEIGHTEDZ(-1, 0, -1); ADD_WEIGHTEDZ( 0, 0, -1); ADD_WEIGHTEDZ( 1, 0, -1);
	ADD_WEIGHTEDZ(-1, 1, -1); ADD_WEIGHTEDZ( 0, 1, -1); ADD_WEIGHTEDZ( 1, 1, -1);

	ADD_WEIGHTEDZ(-1,-1, 0);  ADD_WEIGHTEDZ( 0,-1, 0);  ADD_WEIGHTEDZ( 1,-1, 0);
	ADD_WEIGHTEDZ(-1, 0, 0);  ADD_WEIGHTEDZ( 0, 0, 0);  ADD_WEIGHTEDZ( 1, 0, 0);
	ADD_WEIGHTEDZ(-1, 1, 0);  ADD_WEIGHTEDZ( 0, 1, 0);  ADD_WEIGHTEDZ( 1, 1, 0);

	ADD_WEIGHTEDZ(-1,-1, 1);  ADD_WEIGHTEDZ( 0,-1, 1);  ADD_WEIGHTEDZ( 1,-1, 1);
	ADD_WEIGHTEDZ(-1, 0, 1);  ADD_WEIGHTEDZ( 0, 0, 1);  ADD_WEIGHTEDZ( 1, 0, 1);
	ADD_WEIGHTEDZ(-1, 1, 1);  ADD_WEIGHTEDZ( 0, 1, 1);  ADD_WEIGHTEDZ( 1, 1, 1);
	final[2] = result;

	//debMsg("FINAL","at "<<p<<" = "<<final); // DEBUG
	return final;
}
#undef ADD_WEIGHTEDX
#undef ADD_WEIGHTEDY
#undef ADD_WEIGHTEDZ

inline Real WaveletNoiseField::evaluate(Vec3 pos, int tile) const { 
	pos[0] *= mGsInvX;
	pos[1] *= mGsInvY;
	pos[2] *= mGsInvZ;
	pos += mSeedOffset;

	// time anim
	pos += Vec3(getTime());

	pos[0] *= mPosScale[0];
	pos[1] *= mPosScale[1];
	pos[2] *= mPosScale[2];
	pos += mPosOffset;

	const int n3 = square(NOISE_TILE_SIZE) * NOISE_TILE_SIZE;
	Real v = WNoise(pos, &mNoiseTile[tile*n3]);

	v += mValOffset;
	v *= mValScale;
	if (mClamp) {
		if (v< mClampNeg) v = mClampNeg;
		if (v> mClampPos) v = mClampPos;
	}
	return v;
}

inline Vec3 WaveletNoiseField::evaluateVec(Vec3 pos, int tile) const { 
	pos[0] *= mGsInvX;
	pos[1] *= mGsInvY;
	pos[2] *= mGsInvZ;
	pos += mSeedOffset;

	// time anim
	pos += Vec3(getTime());

	pos[0] *= mPosScale[0];
	pos[1] *= mPosScale[1];
	pos[2] *= mPosScale[2];
	pos += mPosOffset;

	const int n3 = square(NOISE_TILE_SIZE) * NOISE_TILE_SIZE;
	Vec3 v = WNoiseVec(pos, &mNoiseTile[tile*n3]);

	v += Vec3(mValOffset);
	v *= mValScale;
	
	if (mClamp) {
		for(int i=0; i<3; i++) {
			if (v[i]< mClampNeg) v[i] = mClampNeg;
			if (v[i]> mClampPos) v[i] = mClampPos;
		}
	}
	return v;
}

inline Vec3 WaveletNoiseField::evaluateCurl(Vec3 pos) const {
	// gradients of w0-w2
	Vec3 d0 = evaluateVec(pos,0), 
		 d1 = evaluateVec(pos,1), 
		 d2 = evaluateVec(pos,2);
	
	return Vec3(d0.y-d1.z, d2.z-d0.x, d1.x-d2.y);
}

} // namespace  

#endif
