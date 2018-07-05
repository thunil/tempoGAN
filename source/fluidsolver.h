/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Main class for the fluid solver
 *
 ******************************************************************************/

#ifndef _FLUIDSOLVER_H
#define _FLUIDSOLVER_H

#include "manta.h"
#include "vectorbase.h"
#include "vector4d.h"
#include <vector>
#include <map>

namespace Manta { 
	
//! Encodes grid size, timstep etc.
PYTHON(name=Solver) 
class FluidSolver : public PbClass {
public:
	PYTHON() FluidSolver(Vec3i gridSize, int dim=3, int fourthDim=-1);
	virtual ~FluidSolver();
	
	// accessors
	PYTHON() Vec3i getGridSize() { return mGridSize; }
	inline Real  getDt() const      { return mDt; }
	inline Real  getDx() const      { return 1.0 / mGridSize.max(); }
	inline Real  getTime() const    { return mTimeTotal; }

	//! Check dimensionality
	inline bool is2D() const { return mDim==2; }
	//! Check dimensionality (3d or above)
	inline bool is3D() const { return mDim==3; }
	
	PYTHON() void printMemInfo();
	
	//! Advance the solver one timestep, update GUI if present
	PYTHON() void step();
	
	//! Update the timestep size based on given maximal velocity magnitude 
	PYTHON() void adaptTimestep(Real maxVel);
	
	//! create a object with the solver as its parent
	PYTHON() PbClass* create(PbType type, PbTypeVec T=PbTypeVec(),const std::string& name = "");
	
	// temp grid and plugin functions: you shouldn't call this manually
	template<class T> T*   getGridPointer();
	template<class T> void freeGridPointer(T* ptr);    

	//! expose animation time to python
	PYTHON(name=timestep)  Real mDt;  
	PYTHON(name=timeTotal) Real mTimeTotal;
	PYTHON(name=frame)     int  mFrame;
	//! parameters for adaptive time stepping
	PYTHON(name=cfl)          Real mCflCond;  
	PYTHON(name=timestepMin)  Real mDtMin;  
	PYTHON(name=timestepMax)  Real mDtMax;  
	PYTHON(name=frameLength)  Real mFrameLength;

protected:
	Vec3i     mGridSize;
	const int mDim;
	Real      mTimePerFrame;
	bool      mLockDt;
		
	//! subclass for managing grid memory
	//! stored as a stack to allow fast allocation
	template<class T> struct GridStorage {
		GridStorage() : used(0) {}
		T* get(Vec3i size);
		void free();
		void release(T* ptr);
		
		std::vector<T*> grids;
		int used;
	};
	
	//! memory for regular (3d) grids
	GridStorage<int>  mGridsInt;
	GridStorage<Real> mGridsReal;
	GridStorage<Vec3> mGridsVec;


	//! 4d data section, only required for simulations working with space-time data 

public:
	//! 4D enabled? note, there's intentionally no "is4D" function, there are only 3D solvers that also support 4D of a certain size
	inline bool supports4D() const        { return mFourthDim>0; }
	//! fourth dimension size
	inline int  getFourthDim() const { return mFourthDim; }
	//! 4d data allocation
	template<class T> T*   getGrid4dPointer();
	template<class T> void freeGrid4dPointer(T* ptr);    

protected:

	//! 4d size. Note - 4d is not treated like going from 2d to 3d! 4D grids are a separate data type. Normally all
	//! grids are forced to have the same size. In contrast, a solver can create and work with 3D as 
	//! well as 4D grids, when fourth-dim is >0.
	int       mFourthDim;  

	//! 4d grid storage
	GridStorage<Vec4> mGridsVec4; 
	GridStorage<int>  mGrids4dInt;
	GridStorage<Real> mGrids4dReal;
	GridStorage<Vec3> mGrids4dVec;
	GridStorage<Vec4> mGrids4dVec4;
};

}

#endif
