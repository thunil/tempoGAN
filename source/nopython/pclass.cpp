/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * No-python dummy functions
 *
 ******************************************************************************/

#include "manta.h"
#include "general.h"

using namespace std;
namespace Manta {

//******************************************************************************
// Helpers

string PbTypeVec::str() const {
	if (T.empty()) return "";
	string s="<";
	for (int i=0; i<(int)T.size(); i++) {
		s += T[i].str();
		s += (i!=(int)T.size()-1) ? ',' : '>';
	}
	return s;
}
string PbType::str() const {
	if (S=="float") return "Real";
	if (S=="manta.vec3") return "Vec3";
	return S;
}

//******************************************************************************
// PbClass

PbClass::PbClass(FluidSolver* parent, const string& name)
	:  mParent(parent), mName(name), mHidden(false)
{
}

PbClass::PbClass(const PbClass& a) : mParent(a.mParent), mName("_unnamed"), mHidden(false)
{
}

PbClass::~PbClass() 
{
}

void PbClass::checkParent() {
	if (getParent() == NULL) {
		errMsg("New class " + mName + ": no parent given -- specify using parent=xxx !");
	}
}

} // namespace
