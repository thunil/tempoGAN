/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Preprocessor merge file gen
 *
 ******************************************************************************/

#include <map>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include "prep.h"
using namespace std;

struct Chain {
	Chain() {};
	Chain(const string& tpl, const string& target, const string& targetTpl) : 
		tpl(tpl), target(target), targetTpl(targetTpl) {}
	string tpl, target, targetTpl;
};

struct Request {
	Request(const string& c, const string& t) : cls(c), tpl(t), base("") {}
	string cls, tpl, base;
};

struct RegFile {
	RegFile(const string& name, int index) : filename(name), idx(index) {}

	string filename;
	ostringstream out, header, footer;
	vector<Request> req;
	int idx;
};

struct ClassInfo {
	map<string,bool> tplDone;
	vector<string> snippets;
};

static map<string, ClassInfo> classes;
static map<string, Chain> chains;
static vector<RegFile *> regFiles;

// compare template arguments
string mapArgs(const string& inst, const string& match, const string& target) {
	vector<string> inArg = split(inst,',');
	vector<string> maArg = split(match,',');
	vector<string> taArg = split(target,',');
	vector<string> destArg = taArg;
	
	for (int i=0; i<(int)maArg.size(); i++) {
		for (int j=0; j<(int)taArg.size(); j++) {
			if (maArg[i] == taArg[j]) {
				destArg[j] = inArg[i];
			}
		}
	}
	stringstream s;
	for (int i=0; i<(int)destArg.size(); i++) {
		s << destArg[i];
		if (i != (int)destArg.size()-1) s << ',';
	}
	return s.str();
}

void resolveChains(RegFile& file) {
	for (int i=0; i<(int)file.req.size(); i++) {
		Request& req = file.req[i];
		map<string, Chain>::iterator it = chains.find(req.cls);
		if (it != chains.end()) {
			Chain& chain = it->second;
			string tpl = mapArgs(req.tpl, chain.tpl, chain.targetTpl);
			req.base = tpl;
			file.req.push_back(Request(chain.target, tpl));
		}
	}
}

void resolveRequests(RegFile& file) {
	static int FileID = 0;
	// sort request by class
	map<string, vector<Request*> > sortedReqs;
	for (int i=0; i<(int)file.req.size(); i++) {
		Request& req = file.req[i];
		ClassInfo& info = classes[req.cls];
		if (!info.tplDone[req.tpl]) {
			info.tplDone[req.tpl] = true;
			sortedReqs[req.cls].push_back(&req);
		}
	}
	
	bool IsPython = file.filename.find(".py") != std::string::npos;
	stringstream FileidxStr;
	FileidxStr << FileID++;
	file.footer << "extern \"C\" {\n";
	file.footer << "void PbRegister_file_" << FileidxStr.str() << "()\n{\n";

	// process requests
	for(map<string,vector<Request*> >::iterator it = sortedReqs.begin(); it != sortedReqs.end(); ++it) {
		ClassInfo& info = classes[it->first];
		file.out << "#ifdef _C_" << it->first << '\n';
		for (int i=0; i<(int)it->second.size(); i++) {
			Request& req = *(it->second[i]);
			for (int j=0; j<(int)info.snippets.size(); j++) {
				stringstream idxStr;
				idxStr << file.idx++;
				const string table[] = {"CT", req.tpl, "BT", req.base, "IDX", idxStr.str(), ""};
				file.out << replaceSet(info.snippets[j], table) << '\n';
				file.footer << "\tKEEP_UNUSED(_R_" << idxStr.str() << ");\n";

			}
		}
		file.out << "#endif\n";
	}
	if (IsPython) 
		file.footer << "\tKEEP_UNUSED(_reg);\n";
	file.footer << "}\n";
	file.footer << "}";

}

// create data structure from regfiles
void parseLine(const string& line, RegFile& file) {
	if (line.empty()) return;

	vector<string> parts = split(line,'^');
	string cls = parts[0].substr(1);
	
	if (line[0] == '+')
		classes[cls].snippets.push_back(parts[1]);
	else if (line[0] == '>') 
		file.req.push_back(Request(cls,parts[1]));
	else if (line[0] == '@')
		chains[cls] = Chain(parts[1],parts[2],parts[3]);
	else if (line[0] == '#')
		file.header << line << '\n';
	else if (line[0] == '&') {
		string txt = line.substr(1);
		stringstream num; num << file.idx++;
		replaceAll(txt, "$IDX$", num.str());
		file.footer << txt << '\n';
	} else {
		file.out << line << '\n';
	}
}

void generateMerge(int num, char* files[]) {
	// parse files
	for (int i=0; i<num; i++) {
		regFiles.push_back(new RegFile(files[i],i));

		string text = readFile(files[i]);
		replaceAll(text,"\r","");
		vector<string> lines = split(text,'\n');

		for (int j=0; j<(int)lines.size(); j++)
			parseLine(lines[j], *regFiles.back());
	}

	// process and save files
	for (int i=0; i<num; i++) {
		resolveChains(*regFiles[i]);
		resolveRequests(*regFiles[i]);

		string text = "", fn = regFiles[i]->filename;
		bool isPython = fn.compare(fn.size()-7, 7, ".py.reg") == 0;

		if (regFiles[i]->idx > 0) {
			text  = "\n\n\n\n\n// DO NOT EDIT !\n";
			text += "// This file is generated using the MantaFlow preprocessor (prep link).";
			text += "\n\n\n\n\n";
			text += regFiles[i]->header.str();
			text += "namespace Manta {\n";
			text += regFiles[i]->out.str();
			text += regFiles[i]->footer.str();
			text += "}";
		} else if (isPython) {
			text = regFiles[i]->header.str() + regFiles[i]->out.str() + regFiles[i]->footer.str();
		}
		string filename = fn + ".cpp";
		// only write if content is different
		if (!fileExists(filename) || readFile(filename) != text)
			writeFile(filename, text);
		delete regFiles[i];
	}
}
