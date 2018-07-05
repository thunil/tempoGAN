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

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <stack>
#include "prep.h"
using namespace std;

void errMsg(int line, const string& text) {
	cerr << gFilename << ":" << line << ": error: " << text << endl;
	exit(1);
}

void debMsgHelper(int line, const string& text) {
	cout << gFilename << ":" << line << ": info: " << text << endl;
}

bool fileExists(const string& name) {
	ifstream t(name.c_str());
	if (!t.good()) return false;
	t.close();
	return true;
}

string readFile(const string& name) {
	ifstream t(name.c_str());
	if (!t.good()) return "";
	
	string str;
	
	t.seekg(0, ios::end);   
	str.reserve((unsigned)t.tellg());
	t.seekg(0, ios::beg);

	str.assign((istreambuf_iterator<char>(t)),
				istreambuf_iterator<char>());
	
	t.close();
	return str;
}

void writeFile(const string& name, const string& text) {
	ofstream ofs(name.c_str(), ios::binary | ios::out);
	if (!ofs.good()) {
		cerr << "preprocessor error: Can't write to file '" << name << "'" << endl;
		exit(1);
	}
	ofs << text;
	ofs.close();
}

Sink::Sink(const string& infile, const string& outfile):
	infile(infile), filename(outfile)
{
	isHeader = outfile.compare(outfile.size()-3, 3, ".py") == 0 ||
			   outfile.compare(outfile.size()-2, 2, ".h") == 0;
}

void Sink::write() {
	writeFile(filename, inplace.str());
	if (isHeader && !gDocMode) {
		writeFile(filename + ".reg", link.str());
	}
}

vector<string> split(const string& text, char sep) {
	vector<string> bins;
	string cur;
	for (int i=0; i<(int)text.size(); i++) {
		if (text[i]==sep) {
			bins.push_back(cur);
			cur = "";
		} else
			cur += text[i];
	}
	bins.push_back(cur);
	return bins;
}

string strip(const string& s0) {
	string s = s0;
	while (s[0] == ' ') s=s.substr(1);
	while (s[s.size()-1] == ' ') s=s.substr(0,s.size()-1);
	return s;
}

void replaceAll(string& source, string const& find, string const& replace)
{
	for(string::size_type i = 0; (i = source.find(find, i)) != std::string::npos;)
	{
		source.replace(i, find.length(), replace);
		i += replace.length() - find.length() + 1;
	}
}

string makeSafe(const string& s) {
	string t="_X_";
	string source = "+=-<>!()";
	string trans  = "12345678";
	for (int i=0; i<(int)s.size(); i++) {
		size_t idx = source.find(s[i]);
		t += (idx == string::npos) ? s[i] : trans[idx];
	}
	return t;
}

void stealLinebreaks(string& code, int num) {
	// list all line breaks
	vector<int> lb;
	lb.push_back(-1);
	for (int i=0; i<(int)code.size(); i++)
		if (code[i] == '\n')
			lb.push_back(i);
	lb.push_back(code.size());

	for (int i=lb.size()-2; i>0; i--) {
		// make sure we don't mess with comments and defines
		string curLine = code.substr(lb[i-1]+1,lb[i]-lb[i-1]-1);
		string nextLine = code.substr(lb[i]+1, lb[i+1]-lb[i]-1);
		if (curLine.find("//") == string::npos &&
			curLine.find('#') == string::npos &&
			nextLine.find('#') == string::npos) {
			code[lb[i]] = ' ';
			if (--num == 0)
				return;
		}
	}
	return;
}

// Helpers for replaceSet
static string getBracketArg(const string& a, int &pos) {
	string ret="";
	pos++;
	for (;pos<(int)a.size(); pos++) {
		if (a[pos]!='(' && a[pos]!=' ' && a[pos]!='$') break;
	}
	for (; pos<(int)a.size(); pos++) {
		if (a[pos]==')' || a[pos]=='$') return ret;
		ret += a[pos];
	}
	return "";
}

inline bool compareKW(const string& a, int& pos, const string& kw) {
	if (a.compare(pos+1,kw.size(),kw) == 0) {
		pos += kw.size() ;
		return true;
	}
	return false;
}

string replaceSet(const string& templ, const string table[]) {
	vector<string> key, value;
	for (int i=0;!table[i].empty();i+=2) {
		key.push_back(table[i]);
		value.push_back(table[i+1]);
	}
	vector<bool> conditionStack;
	conditionStack.push_back(true);
	stringstream s;
	int elifs = 0;
	for (int i=0; i<(int)templ.size(); i++) {
		char c = templ[i];
		if (c=='@') {
			if (compareKW(templ,i,"IF")) {
				string cond = getBracketArg(templ,i);
				vector<string>::iterator it = find(key.begin(),key.end(),cond);
				bool res = (it != key.end()) && (!value[it-key.begin()].empty());
				conditionStack.push_back(res);
			} else if (compareKW(templ,i,"ELIF")) {
				conditionStack.back() = !conditionStack.back();
				string cond = getBracketArg(templ,i);
				vector<string>::iterator it = find(key.begin(),key.end(),cond);
				bool res = (it != key.end()) && (!value[it-key.begin()].empty());
				conditionStack.push_back(res);
				elifs++;
			} else if (compareKW(templ,i,"END")) {
				for (int k=0; k<elifs+1; k++)
					conditionStack.pop_back();
			} else if (compareKW(templ,i,"ELSE")) {
				conditionStack.back() = !conditionStack.back();
			}
			continue;
		}
		// check condition
		bool valid = true;
		for (int k=0; k<(int)conditionStack.size(); k++)
			if (!conditionStack[k]) valid = false;
		if (!valid) continue;

		if (c=='$') {
			string kw = getBracketArg(templ,i);
			vector<string>::iterator it = find(key.begin(),key.end(),kw);
			if (it == key.end())
				s << '$' << kw << '$';
			else
				s << value[it - key.begin()];
			continue;
		}
		s << templ[i];
		// format output slightly nicer
		if ( gDebugMode && (c=='{' || c=='}')) { s << "\n"; } 
	}
	return s.str();
}
