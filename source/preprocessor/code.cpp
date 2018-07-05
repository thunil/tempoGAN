/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011-2014 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Preprocessor Code Structs
 *
 ******************************************************************************/

#include <iostream>
#include <cstdlib>
#include <algorithm>
#include "prep.h"

using namespace std;

//*************************************************************
// Helpers

bool Type::operator==(const Type& a) const {
	if (a.templateTypes.size() != templateTypes.size()) return false;
	for (size_t i=0; i<templateTypes.size(); i++)
		if (templateTypes[i].name != a.templateTypes[i].name) return false;
	return a.name==name && a.isConst==isConst && a.isRef==isRef && a.isPointer==isPointer;
}

string Type::build(bool refify) const {
	string s="";
	if (isConst) s+= "const ";
	s += name;
	if (!templateTypes.empty())
		s+= templateTypes.minimal;
	if (isRef || (refify && !isPointer)) s+= "&";
	if (isPointer) s+= "*";            
	return s;
}

string Text::linebreaks() const {
	int num = count(original.begin(), original.end(), '\n') - 
		count(minimal.begin(),minimal.end(), '\n');
	string s="";
	for (int i=0; i< num; i++)
		s += '\n';
	return s;
}

string Function::signature() const {
	string s;
	if (isTemplated())
		s = "template " + templateTypes.minimal;
	if (isVirtual)
		s += "virtual ";
	if (isInline)
		s += "inline ";
	s+= returnType.minimal + name + arguments.minimal;
	if (isConst)
		s += "const ";
	return s;    
}

template<> std::string List<Argument>::full(bool refify) const {
	std::stringstream s;
	for (int i=0; i<(int)_data.size(); i++) {
		s << "," << _data[i].type.build(refify) << " " << _data[i].name;
	}
	return s.str();
}

template<> std::string List<Argument>::createMembers(bool refify) const {
	std::stringstream s;
	for (int i=0; i<(int)_data.size(); i++)
		s << _data[i].type.build(refify) << ' ' << _data[i].name << "; ";
	return s.str();
}

template<> std::string List<Argument>::copier(const std::string& prefix, bool useVal) const {
	std::stringstream s;
	for (int i=0; i<(int)_data.size(); i++) {
		s << ',' << _data[i].name << '(';
		s << (useVal ? _data[i].value : prefix+_data[i].name);
		s << ')';
	}
	return s.str();
}

template struct List<Argument>;
template struct List<Type>;

