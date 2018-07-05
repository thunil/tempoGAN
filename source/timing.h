/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Plugin timing
 *
 ******************************************************************************/

#ifndef _TIMING_H
#define _TIMING_H

#include "manta.h"
#include <map>
namespace Manta { 


class TimingData {
private:
	TimingData();
public:
	static TimingData& instance() { static TimingData a; return a; }

	void print();
	void saveMean(const std::string& filename);
	void start(FluidSolver* parent, const std::string& name);
	void stop(FluidSolver* parent, const std::string& name);
protected:
	void step();
	struct TimingSet {
		TimingSet() : num(0),updated(false) { cur.clear(); total.clear(); }
		MuTime cur, total;
		int num;
		bool updated;
		std::string solver;
	};
	bool updated;

	int num;
	MuTime mPluginTimer;
	std::string mLastPlugin;
	std::map<std::string, std::vector<TimingSet> > mData;
};

// Python interface
PYTHON() class Timings : public PbClass {
public:
	PYTHON() Timings() : PbClass(0) {}
	
	PYTHON() void display() { TimingData::instance().print(); }
	PYTHON() void saveMean(std::string file) { TimingData::instance().saveMean(file); }
};

}

#endif
