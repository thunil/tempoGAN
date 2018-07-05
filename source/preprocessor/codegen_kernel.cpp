/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011-2014 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Preprocessor: Process replacement text of KERNEL keyword
 *
 ******************************************************************************/

#include "prep.h"
#include <cstdlib>
#include <set>
#include <sstream>
#include <iostream>
using namespace std;

#define STR(x) #x

//******************************************************
// Templates for code generation

// TP: why do we need getArg? just directly access the argument via its name...
const string TmpAccessor = STR(
inline $TYPE$ getArg$IDX$() { return $NAME$; }
typedef $TYPE_NOREF$ type$IDX$;
);

// Single kernel, default
const string TmpSingleKernel = STR(
$TEMPLATE$ struct $KERNEL$ : public KernelBase {
	$KERNEL$($ARGS$) : 
@IF(PTS) 
		KernelBase($BASE$.size()) $INIT$ $LOCALSET$
@ELSE
		KernelBase($BASE$,$BND$) $INIT$ $LOCALSET$
@END
	{
		runMessage();
		run();
	}
@IF(IJK)
	inline void op(int i, int j, int k, $ARGS$ $LOCALARG$) $CONST$ $CODE$
@ELSE
@IF(FOURD)
	inline void op(int i, int j, int k, int t, $ARGS$ $LOCALARG$) $CONST$ $CODE$
@ELSE
	inline void op(IndexInt idx, $ARGS$ $LOCALARG$) $CONST$ $CODE$
@END
@END

@IF(RET_NAME)
	inline operator $RET_TYPE$() { return $RET_NAME$; }
	inline $RET_TYPE$ & getRet() { return $RET_NAME$; }
@END
	$ACCESSORS$
	$RUNMSG_FUNC$

	$RUN$
	$MEMBERS$
	$LOCALS$
};
);

// Necesary for TBB with nontrivial return values 
const string TmpDoubleKernel = STR(
// inner kernel
$TEMPLATE$ struct _$KERNEL$ : public KernelBase {
	_$KERNEL$(const KernelBase& base, $ARGS$ $LOCALARG$) : 
		KernelBase(base) $INIT$ $LOCALINIT${}

@IF(IJK)
	inline void op(int i, int j, int k, $ARGS$ $LOCALARG$) $CONST$ $CODE$
@ELSE
@IF(FOURD)
	inline void op(int i, int j, int k, int t, $ARGS$ $LOCALARG$) $CONST$ $CODE$
@ELSE
	inline void op(IndexInt idx, $ARGS$ $LOCALARG$) $CONST$ $CODE$
@END
@END
	$RUN$
	$MEMBERS$
	$LOCALS_REF$
};

// outer kernel with accessors
$TEMPLATE$ struct $KERNEL$ : public KernelBase {
	$KERNEL$($ARGS$) :
@IF(PTS) 
		KernelBase($BASE$.size()) $COMMA$ _inner(KernelBase($BASE$.size()),$CALL$)
@ELSE
		KernelBase($BASE$,$BND$) $COMMA$ _inner(KernelBase($BASE$,$BND$),$CALL$)
@END
		$INIT$ $LOCALSET$
	{
		runMessage();
		run();
	}

	void run() { _inner.run(); }

@IF(RET_NAME)
	inline operator $RET_TYPE$() { return $RET_NAME$; }
	inline $RET_TYPE$ & getRet() { return $RET_NAME$; }
@END
	$ACCESSORS$
	$RUNMSG_FUNC$
	_$KERNEL$$TPL$ _inner;
	$MEMBERS$
	$LOCALS$
};
);

const string TmpRunSimple = STR(
void run() {
@IF(IJK)
	const int _maxX = maxX; 
	const int _maxY = maxY;
	for (int k=minZ; k< maxZ; k++)
	for (int j=$BND$; j< _maxY; j++)
	for (int i=$BND$; i< _maxX; i++)
		op(i,j,k, $CALL$);
@ELSE
@IF(FOURD)
	for (int t=minT ; t< maxT; t++)
	for (int k=minZ ; k< maxZ; k++)
	for (int j=$BND$; j< maxY; j++)
	for (int i=$BND$; i< maxX; i++)
		op(i,j,k,t, $CALL$);
@ELSE
	const IndexInt _sz = size;
	for (IndexInt i = 0; i < _sz; i++)
		op(i, $CALL$);
@END
@END
}
);

const string TmpRunTBB = STR(
void operator() (const tbb::blocked_range<IndexInt>& __r) $CONST$ {
@IF(IJK)
	const int _maxX = maxX;
	const int _maxY = maxY;
	if (maxZ>1) {
		for (int k=__r.begin(); k!=(int)__r.end(); k++)
		for (int j=$BND$; j<_maxY; j++)
		for (int i=$BND$; i<_maxX; i++)
			op(i,j,k,$CALL$);
	} else {
		const int k=0;
		for (int j=__r.begin(); j!=(int)__r.end(); j++)
		for (int i=$BND$; i<_maxX; i++)
			op(i,j,k,$CALL$);
	}
@ELSE
@IF(FOURD)
	if (maxT>1) {
		for (int t=__r.begin(); t!=(int)__r.end(); t++)
		for (int k=$BND$; k<maxZ; k++)
		for (int j=$BND$; j<maxY; j++)
		for (int i=$BND$; i<maxX; i++)
			op(i,j,k,t,$CALL$);
	} else if (maxZ>1) {
		const int t=0;
		for (int k=__r.begin(); k!=(int)__r.end(); k++)
		for (int j=$BND$; j<maxY; j++)
		for (int i=$BND$; i<maxX; i++)
			op(i,j,k,t,$CALL$);
	} else {
		const int t=0;
		const int k=0;
		for (int j=__r.begin(); j!=(int)__r.end(); j++)
		for (int i=$BND$; i<maxX; i++)
			op(i,j,k,t,$CALL$);
	}
@ELSE
	for (IndexInt idx=__r.begin(); idx!=(IndexInt)__r.end(); idx++)
		op(idx, $CALL$);
@END
@END
}
void run() {
@IF(IJK)
	if (maxZ>1)
		tbb::parallel_$METHOD$ (tbb::blocked_range<IndexInt>(minZ, maxZ), *this);
	else
		tbb::parallel_$METHOD$ (tbb::blocked_range<IndexInt>($BND$, maxY), *this);
@ELSE
@IF(FOURD)
	if (maxT>1) {
		tbb::parallel_$METHOD$ (tbb::blocked_range<IndexInt>(minT, maxT), *this);
	} else if (maxZ>1) {
		tbb::parallel_$METHOD$ (tbb::blocked_range<IndexInt>(minZ, maxZ), *this);
	} else {
		tbb::parallel_$METHOD$ (tbb::blocked_range<IndexInt>($BND$, maxY), *this); }
@ELSE
	tbb::parallel_$METHOD$ (tbb::blocked_range<IndexInt>(0, size), *this);
@END
@END
}
@IF(REDUCE)
	$IKERNEL$ ($IKERNEL$& o, tbb::split) : KernelBase(o) $COPY$ $LOCALSET$ {}
	
	void join(const $IKERNEL$ & o) {
		$JOINER$
	}
@END
);

const string TmpRunOMP = STR(
void run() {
@IF(IJK)
	const int _maxX = maxX; 
	const int _maxY = maxY;
	if (maxZ > 1) {
		$PRAGMA$ omp parallel $NL$
		{
			$OMP_DIRECTIVE$
			for (int k=minZ; k < maxZ; k++)
			for (int j=$BND$; j < _maxY; j++)
			for (int i=$BND$; i < _maxX; i++)
			   op(i,j,k,$CALL$);
		   $OMP_POST$
		}
	} else {
		const int k=0;
		$PRAGMA$ omp parallel $NL$
		{
			$OMP_DIRECTIVE$
			for (int j=$BND$; j < _maxY; j++)
			for (int i=$BND$; i < _maxX; i++)
				op(i,j,k,$CALL$);
			$OMP_POST$
		}
	}
@ELSE
@IF(FOURD)
	const int _maxX = maxX; 
	const int _maxY = maxY;
	if (maxT > 1) {
		const int _maxZ = maxZ;
		$PRAGMA$ omp parallel $NL$
		{
			$OMP_DIRECTIVE$
			for (int t=$BND$; t < maxT; t++)
			for (int k=$BND$; k < _maxZ; k++)
			for (int j=$BND$; j < _maxY; j++)
			for (int i=$BND$; i < _maxX; i++)
			   op(i,j,k,t,$CALL$);
		   $OMP_POST$
		}
	} else if (maxZ > 1) {
		const int t=0;
		$PRAGMA$ omp parallel $NL$
		{
			$OMP_DIRECTIVE$
			for (int k=minZ; k < maxZ; k++)
			for (int j=$BND$; j < _maxY; j++)
			for (int i=$BND$; i < _maxX; i++)
			   op(i,j,k,t,$CALL$);
		   $OMP_POST$
		}
	} else {
		const int t=0;
		const int k=0;
		$PRAGMA$ omp parallel $NL$
		{
			$OMP_DIRECTIVE$
			for (int j=$BND$; j < _maxY; j++)
			for (int i=$BND$; i < _maxX; i++)
				op(i,j,k,t,$CALL$);
			$OMP_POST$
		}
	}
@ELSE
	const IndexInt _sz = size;
	$PRAGMA$ omp parallel $NL$
	{ 
		$OMP_DIRECTIVE$ 
		for (IndexInt i = 0; i < _sz; i++)
			op(i,$CALL$);
		$OMP_POST$
	}
@END
@END
}
);

const string TmpOMPDirective = STR (
@IF(REDUCE)
	$OMP_PRE$
	$PRAGMA$ omp for nowait $OMP_FOR_OPT$ $NL$
@ELSE
	$PRAGMA$ omp for $OMP_FOR_OPT$ $NL$
@END
);

// hard coded type names for now
const std::string gGridNames[] = {
	std::string("Grid"),
	std::string("MACGrid"),
	std::string("LevelsetGrid"),
	std::string("FlagGrid"),
	std::string("Grid4d") };
const std::string gGrid4dNames[] = {
	std::string("Grid4d") };
static bool isGridCheck(const std::string& type, const std::string* gridNames, int num) { 
	for(int i=0; i<num; ++i) {
		if(type.compare( gridNames[i] )==0) return true;
	}
	return false;
}
bool isGridType(const std::string& type) { 
	return isGridCheck(type, gGridNames, 5);
}
bool isGrid4dType(const std::string& type) { 
	return isGridCheck(type, gGrid4dNames, 1);
}

#define kernelAssert(x,msg) if(!(x)){errMsg(block.line0,string("KERNEL: ") + msg);}

void processKernel(const Block& block, const string& code, Sink& sink) {
	const Function& kernel = block.func;
	
	if (gDocMode) {
		sink.inplace << "//! \\ingroup Kernels\n" << block.func.minimal << "{}\n";
		return;
	}

	// process options
	bool idxMode = false, reduce = false, pts = false, fourdMode = false;
	bool hasLocals = !block.locals.empty(), hasRetType = kernel.returnType.name != "void";
	string bnd = "0", reduceOp="", ompForOpt="";

	MType mtType = gMTType;
	for (size_t i=0; i<block.options.size(); i++) {
		const string& opt = block.options[i].name;
		if (opt == "ijk") 
			idxMode = false;
		else if (opt == "index" || opt == "idx")
			idxMode = true;
		else if (opt == "st" || opt == "single")
			mtType = MTNone;
		else if (opt == "pts" || opt == "particle" || opt == "points")
			pts = true;
		else if (opt == "bnd")
			bnd = block.options[i].value;
		else if (opt == "fourd" )
			fourdMode = true;
		else if (opt == "reduce") {
			reduce = true;
			reduceOp = block.options[i].value;
			if (!(reduceOp == "+" || reduceOp == "-" || reduceOp == "*" ||
				  reduceOp == "/" || reduceOp == "min" || reduceOp == "max"))
				errMsg(block.line0, "invalid 'reduce' operator. Expected reduce= +|-|*|/|min|max");
		} else if (opt == "imbalanced") {
			// The kernels' workload is imbalanced and we need "intelligent" scheduling: 
			// - OpenMP: use chunksize 1 to distribute threads more randomly/evenly
			// - TBB: default (auto_partitioner) is sufficient, do nothing
			ompForOpt.append(" schedule(static,1)"); 
		} else
			errMsg(block.line0, "illegal kernel option '"+ opt +
								"' Supported options are: 'ijk', 'idx', 'bnd=x', 'reduce=x', 'st', 'pts'");
	}
	
	// point out illegal paramter combinations
	kernelAssert (bnd == "0" || !idxMode, "can't combine index mode with bounds iteration.");    
	kernelAssert (!pts || (!idxMode && bnd == "0" ), 
		"KERNEL(opt): Modes 'ijk', 'idx' and 'bnd' can't be applied to particle kernels.");

	// check type consistency of first 'returns' with return type
	if (hasRetType) {
		kernelAssert(hasLocals, "for kernels not returning void a 'returns' statement is required");
		kernelAssert(block.locals.size() == 1, "multiple 'returns' statements only work for 'void' kernels");
		const Type& rt = block.locals[0].type;
		kernelAssert(rt == kernel.returnType, "return type does not match type in first 'returns' statement");
	} else {
		// for void kernels, any number of returns is fine...
		//kernelAssert( block.locals.size() >= 1, "return type specified without a suitable 'returns' statement");
	}
	
	// figure out basegrid
	string baseGrid;
	for (int i=0; i<(int)kernel.arguments.size(); i++) {
		const string& type = kernel.arguments[i].type.name;
		bool isGrid = isGridType(type);
		if (isGrid || pts) { 
			baseGrid = kernel.arguments[i].name;
			if (isGrid && !kernel.arguments[i].type.isPointer)
				baseGrid = "&"+baseGrid;
			break;
		}
	}
	// prevents grids being passed by value
	for (int i=0; i<(int)kernel.arguments.size(); i++) {
		const string& type = kernel.arguments[i].type.name;
		if( isGridType(type) && !(kernel.arguments[i].type.isPointer || kernel.arguments[i].type.isRef) ) {
		errMsg(block.line0, "don't pass grid objects by value!");
		}
	}
	// first arg 4d?
	if(kernel.arguments.size()>0) {
		const string& type = kernel.arguments[0].type.name;
		bool is4d = isGrid4dType(type); 
		if (is4d && (!fourdMode && !idxMode)) {
			errMsg(block.line0, "enable 4d mode to loop over 4d grids!");
		}
	}
	// TODO, potentially add auto check for 4d grid (or data type of first arg in general)

	kernelAssert(!baseGrid.empty(), "use at least one grid to call the kernel.");

	// build accesors
	stringstream accessors;
	for (int i=0; i<(int)kernel.arguments.size(); i++) {
		stringstream num; num << i;
		Type noref = kernel.arguments[i].type;
		noref.isPointer = noref.isRef = noref.isConst = false;
		const string table[] = { "TYPE", kernel.arguments[i].type.build(true),
								 "TYPE_NOREF", noref.build(),
								 "NAME", kernel.arguments[i].name, 
								 "IDX", num.str(), 
								 "" };
		accessors << replaceSet(TmpAccessor, table);
	}

	// optional - print cmd line msg for each kernel call
	stringstream runMsgFunc;
	if(true) {
		runMsgFunc << "void runMessage() { debMsg(\"Executing kernel "+ kernel.name +" \", 3); ";
		runMsgFunc << "debMsg(\"Kernel range\" << ";
		if(!pts) {
			runMsgFunc << " \" x \"<<  maxX  << \" y \"<< maxY  << \" z \"<< minZ<<\" - \"<< maxZ  << \" \"  ";
			if(fourdMode) runMsgFunc << " \" t \"<< minT<<\" - \"<< maxT ";
		} else {
			runMsgFunc << " \" size \"<<  size  << \" \"  ";
		}
		runMsgFunc << " , 4); };";
	}  else { 
		runMsgFunc << "void runMessage() { };"; // disable run msgs
	}

	// build locals, and reduce joiners
	stringstream joiner, preReduce, postReduce;
	for (int i=0; i<(int)block.locals.size(); i++) {
		const string& name = block.locals[i].name;
		const string type = block.locals[i].type.build();
		const string& value = block.locals[i].value;

		preReduce << type << " " << name << " = " << value << ";";
		if (reduceOp == "min" || reduceOp == "max") {
			joiner << name << " = " << reduceOp << "(" << name << ",o." << name << "); ";
			postReduce << "this->" << name << " = " << reduceOp << "(" << name << ", this->" << name << "); ";
		} else {
			joiner << name << " " << reduceOp << "= o." << name << "; ";
			postReduce << "this->" << name << " " << reduceOp << "= " << name << "; ";
		}         
	}
	const string ompPost = reduce ? "\n#pragma omp critical\n{"+postReduce.str()+"}":"";
	bool doubleKernel = mtType == MTTBB && hasRetType && !reduce;
	
	const string table[] = { "IDX", idxMode ? "Y":"",
							 "PTS", pts ? "Y":"",
							 "IJK", (!pts && !idxMode && !fourdMode) ? "Y":"",
							 "FOURD", (fourdMode) ? "Y":"",
							 "REDUCE", reduce ? "Y":"",
							 "TEMPLATE", kernel.isTemplated() ? "template "+kernel.templateTypes.minimal : "",
							 "TPL", kernel.isTemplated() ? "<"+kernel.templateTypes.names()+">" : "",
							 "KERNEL", kernel.name,
							 "IKERNEL", (doubleKernel ? "_":"") + kernel.name,
							 "ARGS", kernel.arguments.listText,
							 "LOCALARG", block.locals.full(true),
							 "BASE", baseGrid,
							 "INIT", kernel.arguments.copier("", false),
							 "LOCALINIT", block.locals.copier("", false),
							 "LOCALSET", block.locals.copier("", true),
							 "COPY", kernel.arguments.copier("o.", false),
							 "MEMBERS", kernel.arguments.createMembers(),
							 "LOCALS", block.locals.createMembers(false),
							 "LOCALS_REF", block.locals.createMembers(true),
							 "ACCESSORS", accessors.str(),
							 "RUNMSG_FUNC", runMsgFunc.str(),
							 "CONST", (!reduce && mtType==MTTBB) ? "const" : "",
							 "CODE", code,
							 "RET_TYPE", hasRetType ? block.locals[0].type.minimal : "",
							 "RET_NAME", hasRetType ? block.locals[0].name : "",
							 "BND", bnd,
							 "CALL", kernel.callString() + (hasLocals ? ","+block.locals.names() : ""),
							 "METHOD", reduce ? "reduce" : "for",
							 "PRAGMA", "\n#pragma",
							 "NL", "\n",
							 "COMMA", ",",
							 "JOINER", joiner.str(),
							 "OMP_PRE", preReduce.str(),
							 "OMP_FOR_OPT", ompForOpt,
							 "OMP_POST", ompPost,
							 "" };

	// generate kernel
	string templ = doubleKernel ? TmpDoubleKernel : TmpSingleKernel;
	if (mtType == MTNone)
		replaceAll(templ, "$RUN$", TmpRunSimple);
	else if (mtType == MTTBB)
		replaceAll(templ, "$RUN$", TmpRunTBB);
	else if (mtType == MTOpenMP) {
		string ompTempl = TmpRunOMP;
		replaceAll(ompTempl, "$OMP_DIRECTIVE$", TmpOMPDirective);
		replaceAll(templ, "$RUN$", ompTempl);
	}

	// synthesize code
	sink.inplace << block.linebreaks() << replaceSet(templ, table);

	// adjust lines after OMP block
	if ( (mtType == MTOpenMP) && (!gDebugMode) )
		sink.inplace << "\n#line " << block.line1 << " \"" << sink.infile << "\"\n" << endl;
}
