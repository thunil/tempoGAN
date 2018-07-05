/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Tools to setup fields and inflows
 *
 ******************************************************************************/

#include "vectorbase.h"
#include "shapes.h"
#include "commonkernels.h"
#include "particle.h"
#include "noisefield.h"
#include "simpleimage.h"
#include "mesh.h"

using namespace std;

namespace Manta {
	
//! Apply noise to grid
KERNEL() 
void KnApplyNoiseInfl(const FlagGrid& flags, Grid<Real>& density, const WaveletNoiseField& noise, const Grid<Real>& sdf, Real scale, Real sigma)
{
	if (!flags.isFluid(i,j,k) || sdf(i,j,k) > sigma) return;
	Real factor = clamp(1.0-0.5/sigma * (sdf(i,j,k)+sigma), 0.0, 1.0);
	
	Real target = noise.evaluate(Vec3(i,j,k)) * scale * factor;
	if (density(i,j,k) < target)
		density(i,j,k) = target;
}

//! Init noise-modulated density inside shape
PYTHON() void densityInflow(const FlagGrid& flags, Grid<Real>& density, const WaveletNoiseField& noise, Shape* shape, Real scale=1.0, Real sigma=0)
{
	Grid<Real> sdf = shape->computeLevelset();
	KnApplyNoiseInfl(flags, density, noise, sdf, scale, sigma);
}
//! Apply noise to real grid based on an SDF
KERNEL() void KnAddNoise(const FlagGrid& flags, Grid<Real>& density, const WaveletNoiseField& noise, const Grid<Real>* sdf, Real scale) {
	if (!flags.isFluid(i,j,k) || (sdf && (*sdf)(i,j,k) > 0.) ) return;
	density(i,j,k) += noise.evaluate(Vec3(i,j,k)) * scale;
}
PYTHON() void addNoise(const FlagGrid& flags, Grid<Real>& density, const WaveletNoiseField& noise, const Grid<Real>* sdf=NULL, Real scale=1.0 ) {
	KnAddNoise(flags, density, noise, sdf, scale );
}

//! sample noise field and set pdata with its values (for convenience, scale the noise values)
KERNEL(pts) template<class T>
void knSetPdataNoise(const BasicParticleSystem& parts, ParticleDataImpl<T>& pdata, const WaveletNoiseField& noise, Real scale) {
	pdata[idx] = noise.evaluate( parts.getPos(idx) ) * scale;
}
KERNEL(pts) template<class T>
void knSetPdataNoiseVec(const BasicParticleSystem& parts, ParticleDataImpl<T>& pdata, const WaveletNoiseField& noise, Real scale) {
	pdata[idx] = noise.evaluateVec( parts.getPos(idx) ) * scale;
}
PYTHON() void setNoisePdata    (const BasicParticleSystem& parts, ParticleDataImpl<Real>& pd, const WaveletNoiseField& noise, Real scale=1.) { knSetPdataNoise<Real>(parts, pd,noise,scale); }
PYTHON() void setNoisePdataVec3(const BasicParticleSystem& parts, ParticleDataImpl<Vec3>& pd, const WaveletNoiseField& noise, Real scale=1.) { knSetPdataNoiseVec<Vec3>(parts, pd,noise,scale); }
PYTHON() void setNoisePdataInt (const BasicParticleSystem& parts, ParticleDataImpl<int >& pd, const WaveletNoiseField& noise, Real scale=1.) { knSetPdataNoise<int> (parts, pd,noise,scale); }

//! SDF gradient from obstacle flags, for turbulence.py
//  FIXME, slow, without kernel...
PYTHON() Grid<Vec3> obstacleGradient(const FlagGrid& flags) {
	LevelsetGrid levelset(flags.getParent(),false);
	Grid<Vec3> gradient(flags.getParent());
	
	// rebuild obstacle levelset
	FOR_IDX(levelset) {
		levelset[idx] = flags.isObstacle(idx) ? -0.5 : 0.5;
	}
	levelset.reinitMarching(flags, 6.0, 0, true, false, FlagGrid::TypeReserved);
	
	// build levelset gradient
	GradientOp(gradient, levelset);
	
	FOR_IDX(levelset) {
		Vec3 grad = gradient[idx];
		Real s = normalize(grad);
		if (s <= 0.1 || levelset[idx] >= 0) 
			grad=Vec3(0.);        
		gradient[idx] = grad * levelset[idx];
	}
	
	return gradient;
}

//! SDF from obstacle flags, for turbulence.py
PYTHON() LevelsetGrid obstacleLevelset(const FlagGrid& flags) {
	LevelsetGrid levelset(flags.getParent(),false);

	// rebuild obstacle levelset
	FOR_IDX(levelset) {
		levelset[idx] = flags.isObstacle(idx) ? -0.5 : 0.5;
	}
	levelset.reinitMarching(flags, 6.0, 0, true, false, FlagGrid::TypeReserved);

	return levelset;
}    


//*****************************************************************************
// blender init functions 

KERNEL() 
void KnApplyEmission(const FlagGrid& flags, Grid<Real>& density, const Grid<Real>& emission, bool isAbsolute)
{
	if (!flags.isFluid(i,j,k) || emission(i,j,k) == 0.) return;
	if (isAbsolute)
		density(i,j,k) = emission(i,j,k);
	else
		density(i,j,k) += emission(i,j,k);
}

//! Add emission values
//isAbsolute: whether to add emission values to existing, or replace
PYTHON() void applyEmission(FlagGrid& flags, Grid<Real>& density, Grid<Real>& emission, bool isAbsolute) {
	KnApplyEmission(flags, density, emission, isAbsolute);
}

// blender init functions for meshes

KERNEL() 
void KnApplyDensity(const FlagGrid& flags, Grid<Real>& density, const Grid<Real>& sdf, Real value, Real sigma)
{
	if (!flags.isFluid(i,j,k) || sdf(i,j,k) > sigma) return;
	density(i,j,k) = value;
}
//! Init noise-modulated density inside mesh
PYTHON() void densityInflowMeshNoise(const FlagGrid& flags, Grid<Real>& density, const WaveletNoiseField& noise, Mesh* mesh, Real scale=1.0, Real sigma=0)
{
	LevelsetGrid sdf(density.getParent(), false);
	mesh->computeLevelset(sdf, 1.);
	KnApplyNoiseInfl(flags, density, noise, sdf, scale, sigma);
}

//! Init constant density inside mesh
PYTHON() void densityInflowMesh(const FlagGrid& flags, Grid<Real>& density, Mesh* mesh, Real value=1., Real cutoff = 7, Real sigma=0)
{
	LevelsetGrid sdf(density.getParent(), false);
	mesh->computeLevelset(sdf, 2., cutoff);
	KnApplyDensity(flags, density, sdf, value, sigma);
}


//*****************************************************************************

//! check for symmetry , optionally enfore by copying
PYTHON() void checkSymmetry( Grid<Real>& a, Grid<Real>* err=NULL, bool symmetrize=false, int axis=0, int bound=0)
{
	const int c  = axis; 
	const int s = a.getSize()[c];
	FOR_IJK(a) { 
		Vec3i idx(i,j,k), mdx(i,j,k);
		mdx[c] = s-1-idx[c];
		if( bound>0 && ((!a.isInBounds(idx,bound)) || (!a.isInBounds(mdx,bound))) ) continue;

		if(err) (*err)(idx) = fabs( (double)(a(idx) - a(mdx) ) ); 
		if(symmetrize && (idx[c]<s/2)) {
			a(idx) = a(mdx);
		}
	}
}
//! check for symmetry , mac grid version
PYTHON() void checkSymmetryVec3( Grid<Vec3>& a, Grid<Real>* err=NULL, bool symmetrize=false , int axis=0, 
								int bound=0, int disable=0)
{
	if(err) err->setConst(0.);

	// each dimension is measured separately for flexibility (could be combined)
	const int c  = axis;
	const int o1 = (c+1)%3;
	const int o2 = (c+2)%3;

	// x
	if(! (disable&1) ) {
		const int s = a.getSize()[c]+1; 
		FOR_IJK(a) { 
			Vec3i idx(i,j,k), mdx(i,j,k);
			mdx[c] = s-1-idx[c]; 
			if(mdx[c] >= a.getSize()[c]) continue; 
			if( bound>0 && ((!a.isInBounds(idx,bound)) || (!a.isInBounds(mdx,bound))) ) continue;

			// special case: center "line" of values , should be zero!
			if(mdx[c] == idx[c] ) {
				if(err) (*err)(idx) += fabs( (double)( a(idx)[c] ) ); 
				if(symmetrize) a(idx)[c] = 0.;
				continue; 
			}

			// note - the a(mdx) component needs to be inverted here!
			if(err) (*err)(idx) += fabs( (double)( a(idx)[c]- (a(mdx)[c]*-1.) ) ); 
			if(symmetrize && (idx[c]<s/2)) {
				a(idx)[c] = a(mdx)[c] * -1.;
			}
		}
	}

	// y
	if(! (disable&2) ) {
		const int s = a.getSize()[c];
		FOR_IJK(a) { 
			Vec3i idx(i,j,k), mdx(i,j,k);
			mdx[c] = s-1-idx[c]; 
			if( bound>0 && ((!a.isInBounds(idx,bound)) || (!a.isInBounds(mdx,bound))) ) continue;

			if(err) (*err)(idx) += fabs( (double)( a(idx)[o1]-a(mdx)[o1] ) ); 
			if(symmetrize && (idx[c]<s/2)) {
				a(idx)[o1] = a(mdx)[o1];
			}
		}
	} 

	// z
	if(! (disable&4) ) {
		const int s = a.getSize()[c];
		FOR_IJK(a) { 
			Vec3i idx(i,j,k), mdx(i,j,k);
			mdx[c] = s-1-idx[c]; 
			if( bound>0 && ((!a.isInBounds(idx,bound)) || (!a.isInBounds(mdx,bound))) ) continue;

			if(err) (*err)(idx) += fabs( (double)( a(idx)[o2]-a(mdx)[o2] ) ); 
			if(symmetrize && (idx[c]<s/2)) {
				a(idx)[o2] = a(mdx)[o2];
			}
		}
	} 

}


// from simpleimage.cpp
void projectImg( SimpleImage& img, const Grid<Real>& val, int shadeMode=0, Real scale=1.);

//! output shaded (all 3 axes at once for 3D)
//! shading modes: 0 smoke, 1 surfaces
PYTHON() void projectPpmFull( const Grid<Real>& val, string name, int shadeMode=0, Real scale=1.)
{
	SimpleImage img;
	projectImg( img, val, shadeMode, scale );
	img.writePpm( name );
}

// helper functions for pdata operator tests

//! init some test particles at the origin
PYTHON() void addTestParts( BasicParticleSystem& parts, int num)
{
	for(int i=0; i<num; ++i)
		parts.addBuffered( Vec3(0,0,0) );

	parts.doCompress();
	parts.insertBufferedParticles();
}

//! calculate the difference between two pdata fields (note - slow!, not parallelized)
PYTHON() Real pdataMaxDiff ( const ParticleDataBase* a, const ParticleDataBase* b )
{    
	double maxVal = 0.;
	//debMsg(" PD "<< a->getType()<<"  as"<<a->getSizeSlow()<<"  bs"<<b->getSizeSlow() , 1);
	assertMsg(a->getType()     == b->getType()    , "pdataMaxDiff problem - different pdata types!");
	assertMsg(a->getSizeSlow() == b->getSizeSlow(), "pdataMaxDiff problem - different pdata sizes!");
	
	if (a->getType() & ParticleDataBase::TypeReal) 
	{
                const ParticleDataImpl<Real>& av = *dynamic_cast<const ParticleDataImpl<Real>*>(a);
                const ParticleDataImpl<Real>& bv = *dynamic_cast<const ParticleDataImpl<Real>*>(b);
		FOR_PARTS(av) {
			maxVal = std::max(maxVal, (double)fabs( av[idx]-bv[idx] ));
		}
	} else if (a->getType() & ParticleDataBase::TypeInt) 
	{
                const ParticleDataImpl<int>& av = *dynamic_cast<const ParticleDataImpl<int>*>(a);
                const ParticleDataImpl<int>& bv = *dynamic_cast<const ParticleDataImpl<int>*>(b);
		FOR_PARTS(av) {
			maxVal = std::max(maxVal, (double)fabs( (double)av[idx]-bv[idx] ));
		}
	} else if (a->getType() & ParticleDataBase::TypeVec3) {
                const ParticleDataImpl<Vec3>& av = *dynamic_cast<const ParticleDataImpl<Vec3>*>(a);
                const ParticleDataImpl<Vec3>& bv = *dynamic_cast<const ParticleDataImpl<Vec3>*>(b);
		FOR_PARTS(av) {
			double d = 0.;
			for(int c=0; c<3; ++c) { 
				d += fabs( (double)av[idx][c] - (double)bv[idx][c] );
			}
			maxVal = std::max(maxVal, d );
		}
	} else {
		errMsg("pdataMaxDiff: Grid Type is not supported (only Real, Vec3, int)");    
	}

	return maxVal;
}


//! calculate center of mass given density grid, for re-centering
PYTHON() Vec3 calcCenterOfMass(const Grid<Real>& density)
{
	Vec3 p(0.0f);
	Real w = 0.0f;
	FOR_IJK(density){
		p += density(i, j, k) * Vec3(i + 0.5f, j + 0.5f, k + 0.5f);
		w += density(i, j, k);
	}
	if (w > 1e-6f)
		p /= w;
	return p;
}


//*****************************************************************************
// helper functions for volume fractions (which are needed for second order obstacle boundaries)



inline static Real calcFraction(Real phi1, Real phi2)
{
	if(phi1>0. && phi2>0.) return 1.;
	if(phi1<0. && phi2<0.) return 0.;

	// make sure phi1 < phi2
	if (phi2<phi1) { Real t = phi1; phi1= phi2; phi2 = t; }
	Real denom = phi1-phi2;
	if (denom > -1e-04) return 0.5; 

	Real frac = 1. - phi1/denom;
	if(frac<0.01) frac = 0.; // stomp small values , dont mark as fluid
	return std::min(Real(1), frac );
}

KERNEL (bnd=1) 
void KnUpdateFractions(const FlagGrid& flags, const Grid<Real>& phiObs, MACGrid& fractions, const int &boundaryWidth) {

	// walls at domain bounds and inner objects
	fractions(i,j,k).x = calcFraction( phiObs(i,j,k) , phiObs(i-1,j,k));
	fractions(i,j,k).y = calcFraction( phiObs(i,j,k) , phiObs(i,j-1,k));
    if(phiObs.is3D()) {
	fractions(i,j,k).z = calcFraction( phiObs(i,j,k) , phiObs(i,j,k-1));
	}

	// remaining BCs at the domain boundaries 
	const int w = boundaryWidth;
	// only set if not in obstacle
 	if(phiObs(i,j,k)<0.) return;

	// x-direction boundaries
	if(i <= w+1) {                     //min x
		if( (flags.isInflow(i-1,j,k)) ||
			(flags.isOutflow(i-1,j,k)) ||
			(flags.isOpen(i-1,j,k)) ) {
				fractions(i,j,k).x = fractions(i,j,k).y = 1.; if(flags.is3D()) fractions(i,j,k).z = 1.;
		}
	}
	if(i >= flags.getSizeX()-w-2) {    //max x
		if(	(flags.isInflow(i+1,j,k)) ||
			(flags.isOutflow(i+1,j,k)) ||
			(flags.isOpen(i+1,j,k)) ) {
			fractions(i+1,j,k).x = fractions(i+1,j,k).y = 1.; if(flags.is3D()) fractions(i+1,j,k).z = 1.;
		}
	}
	// y-direction boundaries
 	if(j <= w+1) {                     //min y
		if(	(flags.isInflow(i,j-1,k)) ||
			(flags.isOutflow(i,j-1,k)) ||
			(flags.isOpen(i,j-1,k)) ) {
			fractions(i,j,k).x = fractions(i,j,k).y = 1.; if(flags.is3D()) fractions(i,j,k).z = 1.;
		}
 	}
 	if(j >= flags.getSizeY()-w-2) {      //max y
		if(	(flags.isInflow(i,j+1,k)) ||
			(flags.isOutflow(i,j+1,k)) ||
			(flags.isOpen(i,j+1,k)) ) {
			fractions(i,j+1,k).x = fractions(i,j+1,k).y = 1.; if(flags.is3D()) fractions(i,j+1,k).z = 1.;
		}
 	}
	// z-direction boundaries
	if(flags.is3D()) {
	if(k <= w+1) {                 //min z
		if(	(flags.isInflow(i,j,k-1)) ||
			(flags.isOutflow(i,j,k-1)) ||
			(flags.isOpen(i,j,k-1)) ) {
			fractions(i,j,k).x = fractions(i,j,k).y = 1.; if(flags.is3D()) fractions(i,j,k).z = 1.;
		}
	}
	if(j >= flags.getSizeZ()-w-2) { //max z
		if(	(flags.isInflow(i,j,k+1)) ||
			(flags.isOutflow(i,j,k+1)) ||
			(flags.isOpen(i,j,k+1)) ) {
			fractions(i,j,k+1).x = fractions(i,j,k+1).y = 1.; if(flags.is3D()) fractions(i,j,k+1).z = 1.;
		}
	}
	}

}

//! update fill fraction values
PYTHON() void updateFractions(const FlagGrid& flags, const Grid<Real>& phiObs, MACGrid& fractions, const int &boundaryWidth=0) {
	fractions.setConst( Vec3(0.) );
	KnUpdateFractions(flags, phiObs, fractions, boundaryWidth);
}

KERNEL (bnd=1) 
void KnUpdateFlagsObs(FlagGrid& flags, const MACGrid* fractions, const Grid<Real>& phiObs, const Grid<Real>* phiOut ) {

	bool isObs = false;
	if(fractions) {
		Real f = 0.;
		f += fractions->get(i  ,j,k).x;
		f += fractions->get(i+1,j,k).x;
		f += fractions->get(i,j  ,k).y;
		f += fractions->get(i,j+1,k).y;
		if (flags.is3D()) {
		f += fractions->get(i,j,k  ).z;
		f += fractions->get(i,j,k+1).z; }
		if(f==0.) isObs = true;
	} else {
		if(phiObs(i,j,k) < 0.) isObs = true;
	}

	bool isOutflow = false;
 	if (phiOut && (*phiOut)(i,j,k) < 0.) isOutflow = true;

 	if (isObs)          flags(i,j,k) = FlagGrid::TypeObstacle;
 	else if (isOutflow) flags(i,j,k) = (FlagGrid::TypeEmpty | FlagGrid::TypeOutflow);
  	else                flags(i,j,k) = FlagGrid::TypeEmpty;
}

//! update obstacle and outflow flags from levelsets
//! optionally uses fill fractions for obstacle
PYTHON() void setObstacleFlags(FlagGrid& flags, const Grid<Real>& phiObs, const MACGrid* fractions=NULL, const Grid<Real>* phiOut=NULL ) {
	KnUpdateFlagsObs(flags, fractions, phiObs, phiOut );
}


//! small helper for test case test_1040_secOrderBnd.py
KERNEL() void kninitVortexVelocity(const Grid<Real> &phiObs, MACGrid& vel, const Vec3 &center, const Real &radius) {
	
	if(phiObs(i,j,k) >= -1.) {

		Real dx = i - center.x; if(dx>=0) dx -= .5; else dx += .5;
		Real dy = j - center.y;
		Real r = std::sqrt(dx*dx+dy*dy);
		Real alpha = atan2(dy,dx);

		vel(i,j,k).x = -std::sin(alpha)*(r/radius);

		dx = i - center.x;
		dy = j - center.y; if(dy>=0) dy -= .5; else dy += .5;
		r = std::sqrt(dx*dx+dy*dy);
		alpha = atan2(dy,dx);

		vel(i,j,k).y = std::cos(alpha)*(r/radius);

	}

}

PYTHON() void initVortexVelocity(const Grid<Real> &phiObs, MACGrid& vel, const Vec3 &center, const Real &radius) {
	kninitVortexVelocity(phiObs,  vel, center, radius);
}


//*****************************************************************************
// helper functions for blurring

//! class for Gaussian Blur
struct GaussianKernelCreator{
public:

	float  mSigma;
	int    mDim;
	float* mMat1D;

	GaussianKernelCreator() : mSigma(0.0f), mDim(0), mMat1D(NULL) {}
	GaussianKernelCreator(float sigma, int dim = 0) 
		: mSigma(0.0f), mDim(0), mMat1D(NULL) {
		setGaussianSigma(sigma, dim);
	}

	Real getWeiAtDis(float disx, float disy){
		float m = 1.0 / (sqrt(2.0 * M_PI) * mSigma);
		float v = m * exp(-(1.0*disx*disx + 1.0*disy*disy) / (2.0 * mSigma * mSigma));
		return v;
	}

	Real getWeiAtDis(float disx, float disy, float disz){
		float m = 1.0 / (sqrt(2.0 * M_PI) * mSigma);
		float v = m * exp(-(1.0*disx*disx + 1.0*disy*disy + 1.0*disz*disz) / (2.0 * mSigma * mSigma));
		return v;
	}

	void setGaussianSigma(float sigma, int dim = 0){
		mSigma = sigma;
		if (dim < 3)
			mDim = (int)(2.0 * 3.0 * sigma + 1.0f);
		else
			mDim = dim;
		if (mDim < 3) mDim = 3;
				
		if (mDim % 2 == 0) ++mDim;// make dim odd  

		float s2 = mSigma * mSigma;
		int c = mDim / 2;
		float m = 1.0 / (sqrt(2.0 * M_PI) * mSigma);

		// create 1D matrix
		if (mMat1D) delete[] mMat1D;
		mMat1D = new float[mDim];
		for (int i = 0; i < (mDim + 1) / 2; i++){
			float v = m * exp(-(1.0*i*i) / (2.0 * s2));
			mMat1D[c + i] = v;
			mMat1D[c - i] = v;
		}
	}

	~GaussianKernelCreator(){
		if (mMat1D) delete[] mMat1D;
	}

	float get1DKernelValue(int off){
		assertMsg(off >= 0 && off < mDim, "off exceeded boundary in Gaussian Kernel 1D!");
		return mMat1D[off];
	}

};

template<class T>
T convolveGrid(Grid<T>& originGrid, GaussianKernelCreator& gkSigma, Vec3 pos, int cdir){
	// pos should be the centre pos, e.g., 1.5, 4.5, 0.5 for grid pos 1,4,0
	Vec3 step(1.0, 0.0, 0.0);
	if (cdir == 1)// todo, z
		step = Vec3(0.0, 1.0, 0.0);
	else if (cdir == 2)
		step = Vec3(0.0, 0.0, 1.0);
	T pxResult(0);
	for (int i = 0; i < gkSigma.mDim; ++i){
		Vec3i curpos = toVec3i(pos - step*(i - gkSigma.mDim / 2));
		if (originGrid.isInBounds(curpos))
			pxResult += gkSigma.get1DKernelValue(i) * originGrid.get(curpos);
		else{ // TODO , improve...
			Vec3i curfitpos = curpos;
			if (curfitpos.x < 0) curfitpos.x = 0;
			else if (curfitpos.x >= originGrid.getSizeX()) curfitpos.x = originGrid.getSizeX() - 1;
			if (curfitpos.y < 0) curfitpos.y = 0;
			else if (curfitpos.y >= originGrid.getSizeY()) curfitpos.y = originGrid.getSizeY() - 1;
			if (curfitpos.z < 0) curfitpos.z = 0;
			else if (curfitpos.z >= originGrid.getSizeZ()) curfitpos.z = originGrid.getSizeZ() - 1;
			pxResult += gkSigma.get1DKernelValue(i) * originGrid.get(curfitpos);
		}
	}
	return pxResult;
}

KERNEL() template<class T>
void knBlurGrid(Grid<T>& originGrid, Grid<T>& targetGrid, GaussianKernelCreator& gkSigma, int cdir){
	targetGrid(i, j, k) = convolveGrid<T>(originGrid, gkSigma, Vec3(i, j, k), cdir);
}

template<class T>
int blurGrid(Grid<T>& originGrid, Grid<T>& targetGrid, float sigma){
	GaussianKernelCreator tmGK(sigma);
	Grid<T> tmpGrid(originGrid);
	knBlurGrid<T>(originGrid, tmpGrid, tmGK, 0); //blur x
	knBlurGrid<T>(tmpGrid, targetGrid, tmGK, 1); //blur y
	if (targetGrid.is3D()){
		tmpGrid.copyFrom(targetGrid);
		knBlurGrid<T>(tmpGrid, targetGrid, tmGK, 2);
	}
	return tmGK.mDim;
}


KERNEL()
void KnBlurMACGridGauss(MACGrid& originGrid, MACGrid& target, GaussianKernelCreator& gkSigma, int cdir){
	Vec3 pos(i, j, k);
	Vec3 step(1.0, 0.0, 0.0);
	if (cdir == 1)
		step = Vec3(0.0, 1.0, 0.0);
	else if (cdir == 2)
		step = Vec3(0.0, 0.0, 1.0);
	
	Vec3 pxResult(0.0f);
	for (int di = 0; di < gkSigma.mDim; ++di){
		Vec3i curpos = toVec3i(pos - step*(di - gkSigma.mDim / 2));
		if (!originGrid.isInBounds(curpos)){
			if (curpos.x < 0) curpos.x = 0;
			else if (curpos.x >= originGrid.getSizeX()) curpos.x = originGrid.getSizeX() - 1;
			if (curpos.y < 0) curpos.y = 0;
			else if (curpos.y >= originGrid.getSizeY()) curpos.y = originGrid.getSizeY() - 1;
			if (curpos.z < 0) curpos.z = 0;
			else if (curpos.z >= originGrid.getSizeZ()) curpos.z = originGrid.getSizeZ() - 1;
		}
		pxResult += gkSigma.get1DKernelValue(di) * originGrid.get(curpos);
	}
	target(i,j,k) = pxResult;
}

PYTHON() int blurMacGrid(MACGrid& oG, MACGrid& tG, float si){
	GaussianKernelCreator tmGK(si); 
	MACGrid tmpGrid(oG);
	KnBlurMACGridGauss(oG, tmpGrid, tmGK, 0); //blur x
	KnBlurMACGridGauss(tmpGrid, tG, tmGK, 1); //blur y
	if (tG.is3D()){
		tmpGrid.copyFrom(tG);
		KnBlurMACGridGauss(tmpGrid, tG, tmGK, 2);
	}
	return tmGK.mDim;
}

PYTHON() int blurRealGrid(Grid<Real>& oG, Grid<Real>& tG, float si){
	return blurGrid<Real> (oG, tG, si);
}

} // namespace

