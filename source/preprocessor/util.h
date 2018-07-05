/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011-2014 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Preprocessor Declarations
 *
 ******************************************************************************/

#ifndef _P_UTIL_H
#define _P_UTIL_H

#include <string>
#include <sstream>
#include <vector>
 
// warnings, errors
void errMsg(int line, const std::string& text);

// text tools
void replaceAll(std::string& text, const std::string& pattern, const std::string& repl);
std::string replaceSet(const std::string& templ, const std::string table[]);
std::vector<std::string> split(const std::string& text, char sep);
void stealLinebreaks(std::string& code, int num);
std::string strip(const std::string& s);
std::string makeSafe(const std::string& s);
inline bool isNameChar(char c) {
	return (c>='A' && c<='Z') || (c>='a' && c<='z') || (c>='0' && c<='9') || c=='_';
}

// file tools
std::string readFile(const std::string&);
bool fileExists(const std::string& name);
void writeFile(const std::string& name, const std::string& text);

struct Sink {
	Sink(const std::string& infile,const std::string& outfile);
	void write();

	std::ostringstream inplace;
	std::ostringstream link;
	bool isHeader;
	std::string infile;
private:
	std::string filename;
};

// simple string-based stack
struct BracketStack {
	bool empty() { return stack.empty(); }
	char top() { return empty() ? '\0' : *(stack.rbegin()); }
	void push_back(char c) { stack += c; }
	char pop() { if (empty()) return '\0'; char c = *(stack.rbegin()); stack.erase(stack.end()-1); return c; }
	
	std::string stack;
};


#endif