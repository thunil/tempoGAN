/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Turbulence particles
 *
 ******************************************************************************/

#ifndef _TURBULENCEPART_H_
#define _TURBULENCEPART_H_

#include "particle.h"
#include "noisefield.h"

namespace Manta {
class Shape;
	

	
struct TurbulenceParticleData {
	TurbulenceParticleData() : pos(0.0),color(1.),tex0(0.0),tex1(0.0),flag(0) {}
	TurbulenceParticleData(const Vec3& p, const Vec3& color = Vec3(1.)) : pos(p),color(color),tex0(p),tex1(p),flag(0) {}
	Vec3 pos, color;
	Vec3 tex0, tex1;
	int flag;
	static ParticleBase::SystemType getType() { return ParticleBase::TURBULENCE; }
};

//! Turbulence particles
PYTHON() class TurbulenceParticleSystem : public ParticleSystem<TurbulenceParticleData> {
public:
	PYTHON() TurbulenceParticleSystem(FluidSolver* parent, WaveletNoiseField& noise);
  
	PYTHON() void resetTexCoords(int num, const Vec3& inflow);    
	PYTHON() void seed(Shape* source, int num);
	PYTHON() void synthesize(FlagGrid& flags, Grid<Real>& k, int octaves=2, Real switchLength=10.0, Real L0=0.1, Real scale=1.0, Vec3 inflowBias=0.0);
	PYTHON() void deleteInObstacle(FlagGrid& flags);
		
	virtual ParticleBase* clone();
	
private:
	WaveletNoiseField& noise;
};

} // namespace


#endif
