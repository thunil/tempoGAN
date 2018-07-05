/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Use this file to test new functionality
 *
 ******************************************************************************/

#include "levelset.h"
#include "commonkernels.h"
#include "particle.h"
#include <cmath>

using namespace std;

namespace Manta {

// two simple example kernels

KERNEL(idx, reduce=+) returns (double sum=0)
double reductionTest(const Grid<Real>& v)
{
	sum += v[idx];
}

KERNEL(idx, reduce=min) returns (double sum=0)
double minReduction(const Grid<Real>& v)
{
	if (sum < v[idx])
		sum = v[idx];
}



//! for data generation scenes
PYTHON() Vec3 centrePosofGrid(Grid<Real>& gridden){
	Vec3 centrePos(0.0f);
	Real centreWei = 0.0f;
	FOR_IJK(gridden){
		centrePos += gridden(i, j, k) * Vec3(i + 0.5f, j + 0.5f, k + 0.5f);
		centreWei += gridden(i, j, k);
	}
	if (centreWei > 1e-6f)
		centrePos /= centreWei;
	return centrePos;
}




// ... add more test code here if necessary ...

} //namespace

