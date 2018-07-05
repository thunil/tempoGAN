/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Levelset
 *
 ******************************************************************************/

#ifndef _LEVELSET_H_
#define _LEVELSET_H_

#include "grid.h"

namespace Manta {
class Mesh;

//! Special function for levelsets
PYTHON() class LevelsetGrid : public Grid<Real> {
public:
	PYTHON() LevelsetGrid(FluidSolver* parent, bool show = true);

        LevelsetGrid(FluidSolver* parent, Real* data, bool show = true);
	
	//! reconstruct the levelset using fast marching
	PYTHON() void reinitMarching(const FlagGrid& flags, Real maxTime=4.0, 
			MACGrid* velTransport=NULL, bool ignoreWalls=false, bool correctOuterLayer=true, 
			int obstacleType = FlagGrid::TypeObstacle );

	//! create a triangle mesh from the levelset isosurface
	PYTHON() void createMesh(Mesh& mesh);
	
	//! union with another levelset
	PYTHON() void join(const LevelsetGrid& o);
	PYTHON() void subtract(const LevelsetGrid& o);
	
	//! initialize levelset from flags (+/- 0.5 heaviside)
	PYTHON() void initFromFlags(const FlagGrid& flags, bool ignoreWalls=false);
	
	static Real invalidTimeValue();
};

} //namespace
#endif
