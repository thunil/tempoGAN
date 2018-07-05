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

#ifndef _CODE_H
#define _CODE_H

#include <string>
#include <vector>

struct Text {
	int line0;
	std::string minimal, original;
	void reset() { minimal = original = ""; line0=0; }
	std::string linebreaks() const; 
	virtual std::string dynamicClass() { return ""; }
	void prequel(const Text* a) { minimal = a->minimal + minimal; original = a->original + original; }
	void post(const Text* a) { minimal += a->minimal; original += a->original; }
};

template <class T>
struct List : Text {
	std::vector<T> _data;
	std::string listText;

	// exposed vector members
	inline size_t size() const { return _data.size(); }
	inline T& operator[](int i) { return _data[i]; }
	inline const T& operator[](int i) const { return _data[i]; }
	inline bool empty() const { return _data.empty(); }
	inline T& back() { return _data.back(); }
	inline void push_back(const T& a) { _data.push_back(a); }

	std::string names() const;
	std::string full(bool refify=false) const { return ""; };
	std::string copier(const std::string& prefix, bool useValue) const { return ""; };
	std::string createMembers(bool refify = false) const { return ""; };
	virtual std::string dynamicClass() { return "List"; }
};

struct Type : Text {
	Type() : isConst(false), isRef(false), isPointer(false) {};
	
	std::string name;
	bool isConst, isRef, isPointer;
	List<Type> templateTypes;

	inline bool isTemplated() const { return !templateTypes.empty(); }
	bool operator==(const Type& a) const;
	std::string build(bool refify = false) const;
	virtual std::string dynamicClass() { return "Type"; }
};

struct Argument : Text {
	Argument() : type(),index(-1) {};
	
	Type type;
	int index;
	std::string name, value;
	std::string completeText, minimalText;
	virtual std::string dynamicClass() { return "Argument"; }
};

struct Function : Text {
	Function() : returnType(),isInline(false),isVirtual(false),isConst(false),noParentheses(false),isOperator(false) {}

	std::string name;
	Type returnType;
	bool isInline, isVirtual, isConst, noParentheses, isOperator;
	List<Type> templateTypes;
	List<Argument> arguments;
	std::string signature() const;
	inline std::string callString() const { return arguments.names(); }
	inline bool isTemplated() const { return !templateTypes.empty(); }
	virtual std::string dynamicClass() { return "Function"; }
};

struct Instantiation {
	Instantiation(const std::string& c, const std::string& n) : 
		cls(c), name(n) {}
	std::string cls, name;
	std::vector<std::string> wrapName;
	std::vector<Function> func;
	std::vector<std::string> templates;    
};

struct Class : Text {
	Class() {};

	std::string name;
	Type baseClass;
	List<Type> templateTypes;  
	std::string fullName() const { return isTemplated() ? (name+"<"+templateTypes.names()+">") : name; }
	inline bool isTemplated() const { return !templateTypes.empty(); }
	virtual std::string dynamicClass() { return "Class"; }
};

struct Block : Text {
	Block() {};

	int line1;
	std::string initList;
	Class cls;
	const Class *parent;
	Function func;
	List<Argument> options;
	List<Argument> locals;
	virtual std::string dynamicClass() { return "Block"; }
};

// list functions , notify compiler of specializations in code.cpp

template<class T> std::string List<T>::names() const {
	std::stringstream s;
	for (int i=0; i<(int)_data.size(); i++) {
		s << _data[i].name;
		if (i != (int)_data.size()-1) s << ',';
	}
	return s.str();
} 
template<> std::string List<Argument>::full(bool refify) const;
template<> std::string List<Argument>::createMembers(bool refify) const;
template<> std::string List<Argument>::copier(const std::string& prefix, bool useVal) const;

#endif
