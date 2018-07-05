/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Preprocessor Main
 *
 ******************************************************************************/

#include <string>
#include <sstream>
#include <fstream>
#include <streambuf>
#include <iostream>
#include <cstdlib>
#include <cstring>

#include "prep.h"

using namespace std;

string gFilename;
bool gDebugMode = false;
bool gDocMode;
bool gIsHeader;
MType gMTType = MTNone;


void usage() {
	cerr << "preprocessor error: Unknown parameters." << endl;
	cerr << "  Usage : prep generate <dbg_mode> <mt_type> <inputdir> <inputfile> <outputfile>" << endl;
	cerr << "     or : prep docgen <dbg_mode> <mt_type> <inputdir> <inputfile> <outputfile>" << endl;
	cerr << "     or : prep link <regfiles...>" << endl;
	exit(1);
}

void doMerge(int argc, char* argv[]) {    
	if (argc < 3) usage();
	
	generateMerge(argc-2, &argv[2]);
}

void doGenerate(int argc, char* argv[], bool docs) {
	gDocMode = docs;
	gDebugMode = false;
	gMTType    = MTNone;
	if (argc != 7) usage();
	
	// set constants
	const string indir(argv[4]), infile(argv[5]), outfile(argv[6]);
	bool isPython = infile.size() > 3 && !infile.compare(infile.size()-3, 3, ".py");

	// only enable in cmake's PREP_DEBUG mode (passed via cmd line option dbg_mode)
	gDebugMode = atoi(argv[2]) != 0;
	if (!strcmp(argv[3],"TBB")) gMTType = MTTBB;
	if (!strcmp(argv[3],"OPENMP")) gMTType = MTOpenMP;
	
	// load complete file into buffer    
	gFilename = indir+infile;
	string text = readFile(gFilename);
	if (text.empty()) {
		cerr << "preprocessor error: Can't read file '" << infile << "'" << endl;
		exit(1);
	}
	// pad text for easier lexing lookups
	text += "\n\n\n";
	
	Sink sink(infile,outfile);
	if (gDocMode) {
		sink.inplace << "/*! \\file " + infile + " */\n";
	} else {
		sink.inplace << "\n\n\n\n\n// DO NOT EDIT !\n";
		sink.inplace << "// This file is generated using the MantaFlow preprocessor (prep generate).";
		sink.inplace << "\n\n\n\n\n";
	}
	
	if (isPython) {
		// python file, only registering
		replaceAll(text, "\n", "\\n");
		replaceAll(text, "\r", "");
		replaceAll(text, "\t", "\\t");
		replaceAll(text, "\"", "<qtm>"); // split into two ops to avoid interference
		replaceAll(text, "<qtm>", "\\\"");
		sink.link << "#include \"registry.h\"\n";
		sink.link << "static const Pb::Register _reg(\"" + infile + "\", \"" + text + "\");\n";
	} else {
		if (!gDocMode) {
			sink.link << "#include \"" + infile + "\"\n";
			if (!gDebugMode)
				sink.inplace << "#line 1 \"" << indir << infile << "\"\n";
		}
		std::vector<Instantiation> inst;
		processText(text, 1, sink, 0, inst);
		postProcessInstantiations(sink, inst);
	}
	sink.write();
}

void doRegister(int argc, char* argv[]) {
	if(argc<3) { 
		cerr << "Wrong call for prep register, not enough arguments" << endl; 
		exit(1); 
	}
	std::ofstream output(argv[argc-1]); // full path from call
	std::string newl("\n"); 
	
	std::stringstream RegistrationsDefs;
	std::stringstream Registrations;
	for (int i = 2; i < argc-1; i++) {
		std::ifstream input(argv[i]);

		for (std::string line; getline(input, line); )	{
			int pos = line.find("void PbRegister_");
			if (pos != (int)std::string::npos) {
				std::string lineRegEnd = line.substr( pos, line.length() );
				int endpos = lineRegEnd.find("{");
				std::string lineFunc = lineRegEnd.substr( 0, endpos );

				RegistrationsDefs << "\t\textern " << lineFunc << ";" << newl;
				replaceAll(lineFunc, "void ", "");
				Registrations << "\t\t" << lineFunc << ";" << newl;
			}
		}
		input.close();
	}

	// add external declarations
	output << "extern \"C\" {\n";
	output << RegistrationsDefs.str();
	output << "}\n\n";

	// add dummy function calling all of them
	output << "namespace Pb {\n";
	output << "\tvoid MantaEnsureRegistration()\n\t{\n";
	output << Registrations.str();
	output << "\t}\n";
	output << "}\n";
	output.close();
}


int main(int argc, char* argv[]) {
	// command line options
	if (argc < 2) usage();

	// prep modes:
	// - preprocessing of header / cpp files, "main" mode
	// - generating doxygen info
	// - merge / link generated functions
	// - create dummy registry call to prevent linkers from discarding python bindings
	if (!strcmp(argv[1],"generate")) 
		doGenerate(argc, argv, false);
	else if (!strcmp(argv[1],"docgen"))
		doGenerate(argc, argv, true);
	else if (!strcmp(argv[1],"link")) 
		doMerge(argc, argv);
	else if (!strcmp(argv[1], "register"))
		doRegister(argc, argv);
	else 
		usage();
	
	return 0;
}
 
