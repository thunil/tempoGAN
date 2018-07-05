/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * moving obstacles
 *
 ******************************************************************************/

#ifndef _MOVINGOBS_H
#define _MOVINGOBS_H

#include "shapes.h"
#include "particle.h"

namespace Manta {

//! Moving obstacle composed of basic shapes
PYTHON() class MovingObstacle : public PbClass {
public:
	PYTHON() MovingObstacle(FluidSolver* parent, int emptyType=FlagGrid::TypeEmpty);
	
	PYTHON() void add(Shape* shape);
	//! If t in [t0,t1], apply linear motion path from p0 to p1
	PYTHON() void moveLinear(Real t, Real t0, Real t1, Vec3 p0, Vec3 p1, FlagGrid& flags, MACGrid& vel, bool smooth=true);
	//! Compute levelset, and project FLIP particles outside obstacles
	PYTHON() void projectOutside(FlagGrid& flags, BasicParticleSystem& flip);
	
protected:
	std::vector<Shape*> mShapes;
	int mEmptyType;
	int mID;
	static int sIDcnt;
};
	

} //namespace
#endif
