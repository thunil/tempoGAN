/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011-2014 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Preprocessor Parsing
 *
 ******************************************************************************/

#include <iostream>
#include <cstdlib>
#include <algorithm>
#include "prep.h"

using namespace std;

List<Type> parseTypeList(TokenPointer& parentPtr);

//*************************************************************
// parsers

#define tkAssert(x,msg) {if(!(x)){tk.errorMsg(msg);}}
#define typeAssert(x) tkAssert(x," is not a valid type.")
#define argAssert(x) tkAssert(x," is not a valid argument.")

string parseRunaway(TokenPointer& parentPtr) {
	Text text;
	TokenPointer tk(parentPtr, &text);

	// don't validate, just track bracket level
	BracketStack stack;
	for(;!tk.done();tk.next()) {
		if (stack.empty() && (tk.curType() == TkComma || 
			tk.curType() == TkBracketR || tk.curType() == TkTBracketR)) {
				break;
		}
		if (tk.curType() == TkBracketL || tk.curType() == TkTBracketL) {
			stack.push_back(tk.cur().text[0]);
		} else if (tk.curType() == TkBracketR) {
			argAssert(stack.pop() == '(');
		} else if (tk.curType() == TkTBracketR) {
			argAssert(stack.pop() == '<');
		}
	}
	argAssert(stack.empty());
	return text.minimal;
}

void parsePointer(TokenPointer& tk, Type& cur) {
	if (tk.done()) return;
	if (tk.cur().type == TkOperator && tk.cur().text == "*") {
		cur.isPointer = true;
		tk.next();
	} else if (tk.cur().type == TkRef) {
		cur.isRef = true;
		tk.next();
	}
}

Type parseType(TokenPointer& parentPtr) {
	Type cur = Type();
	TokenPointer tk(parentPtr, &cur);

	// constness
	if (tk.curType() == TkConst) {
		cur.isConst = true;
		tk.next();
	}
	typeAssert(!tk.done());
	
	// signed / unsigned
	if (tk.curType() == TkTypeQualifier) {
		cur.name = tk.cur().text;
		tk.next();
		if (tk.done())
			return cur;
		typeAssert(tk.curType() == TkSimpleType);
		cur.name += " " + tk.cur().text;
		tk.next();
		parsePointer(tk, cur);
		return cur;
	}

	// template argument
	if (tk.curType() == TkClass) {
		tk.next();
	}

	typeAssert(tk.curType() == TkDescriptor || tk.curType() == TkSimpleType);
	cur.name = tk.cur().text;
	tk.next();
	if (cur.name == "operator") {
		while (tk.curType() == TkOperator) {
			cur.name += tk.cur().text;
			tk.next();
		}
	}
	
	// namespace
	if (tk.curType() == TkDoubleColon) {
		cur.name += "::";
		tk.next();
		typeAssert(tk.curType() == TkDescriptor || tk.curType() == TkSimpleType);
		cur.name += tk.cur().text;
		tk.next();
	}

	// template
	if (tk.curType() == TkTBracketL) {
		cur.templateTypes = parseTypeList(tk);
	}

	parsePointer(tk, cur);
	return cur;
}

Argument parseArgument(TokenPointer& parentPtr, bool requireName, bool requireType) {
	Argument cur = Argument();
	TokenPointer tk(parentPtr, &cur);

	if (requireType)
		cur.type = parseType(tk);

	if (tk.curType() != TkDescriptor && !requireName) 
		return cur;

	argAssert(tk.curType() == TkDescriptor);
	cur.name = tk.cur().text;
	tk.next();

	// default value ?
	if (tk.curType() == TkBracketL || (tk.curType() == TkOperator && tk.cur().text == "=")) {
		if (tk.curType() == TkOperator)
			tk.next();
		cur.value = parseRunaway(tk);
	}
	return cur;
}

List<Type> parseTypeList(TokenPointer& parentPtr) {
	List<Type> list;
	TokenPointer tk(parentPtr, &list);
	
	tkAssert(tk.curType() == TkTBracketL, "expected template opening bracket");
	tk.next();
	if (tk.curType() != TkTBracketR) {
		for(;;) {
			list.push_back(parseType(tk));
			if (tk.curType() == TkTBracketR) 
				break;
			tkAssert(tk.curType() == TkComma, "expected comma or closing bracket");
			tk.next();
		}
	}
	list.listText = list.minimal.substr(1);
	tkAssert(tk.curType() == TkTBracketR, "expected template closing bracket");
	tk.next();
	return list;
}

List<Argument> parseArgumentList(TokenPointer& parentPtr, bool requireName, bool requireType) {
	List<Argument> list;
	TokenPointer tk(parentPtr, &list);
	
	tkAssert(tk.curType() == TkBracketL, "expected opening bracket");
	tk.next();
	if (tk.curType() != TkBracketR) {
		for(int idx=0;;idx++) {
			list.push_back(parseArgument(tk, requireName, requireType));
			list.back().index = idx;
			if (tk.curType() == TkBracketR) 
				break;
			tkAssert(tk.curType() == TkComma, "expected comma or closing bracket");
			tk.next();
		}
	}
	list.listText = list.minimal.substr(1);
	tkAssert(tk.curType() == TkBracketR, "expected closing bracket");
	tk.next();
	return list;
}

Function parseFunction(TokenPointer& parentPtr, bool requireNames, bool requireType, bool requireArgs) {
	Function cur;
	TokenPointer tk(parentPtr, &cur);

	// templated
	if (tk.curType() == TkTemplate) {
		tk.next();
		cur.templateTypes = parseTypeList(tk);            
	}

	for (;tk.curType() == TkInline || tk.curType() == TkVirtual; tk.next()) {
		if (tk.curType() == TkInline) cur.isInline = true;
		if (tk.curType() == TkVirtual) cur.isVirtual = true;
	}

	if (requireType)
		cur.returnType = parseType(tk);
	tkAssert(tk.curType() == TkDescriptor, "malformed function/kernel");
	cur.name = tk.cur().text;
	tk.next();

	if (cur.name == "operator") {
		cur.isOperator = true;
		while(tk.curType() == TkOperator) {
			cur.name += tk.cur().text;
			tk.next();
		}
	}

	if (requireArgs || tk.curType() == TkBracketL)
		cur.arguments = parseArgumentList(tk, requireNames, true);
	else
		cur.noParentheses = true;

	if (tk.curType() == TkConst) {
		cur.isConst = true;
		tk.next();
	}
	return cur;
}

Class parseClass(TokenPointer& parentPtr) {
	Class cur;
	TokenPointer tk(parentPtr, &cur);

	tkAssert(tk.curType() == TkClass, "");
	tk.next();
	tkAssert(tk.curType() == TkDescriptor, "malformed preprocessor keyword block. Expected 'PYTHON class name : public X {}'");
	cur.name = tk.cur().text;
	tk.next();
	tkAssert(tk.curType() == TkColon, "PYTHON class must publicly derive from PbClass (or a subclass)");
	tk.next();
	tkAssert(tk.curType() == TkPublic, "PYTHON class must publicly derive from PbClass (or a subclass)");
	tk.next();
	tkAssert(tk.curType() == TkDescriptor, "PYTHON class must publicly derive from PbClass (or a subclass)");
	cur.baseClass = parseType(tk);

	return cur;
}

// Parse syntax KEYWORD(opt1, opt2, ...) STATEMENTS [ {} or ; ]    
void parseBlock(const string& kw, const vector<Token>& tokens, const Class* parent, Sink& sink, vector<Instantiation>& inst) {
	Block block = Block();
	block.parent = parent;
	TokenPointer tk(tokens, &block);

	// parse keyword options
	if (tk.curType() == TkBracketL)
		block.options = parseArgumentList(tk, true, false);

	if (kw == "KERNEL") {
		List<Type> templTypes;

		// templated kernel
		if (tk.curType() == TkTemplate) {
			tk.next();
			templTypes = parseTypeList(tk);            
		}
		
		// return values
		while (tk.curType() == TkDescriptor && tk.cur().text == "returns") {
			tk.next();
			tkAssert(tk.curType() == TkBracketL, "expext opening bracket");
			tk.next();
			block.locals.push_back(parseArgument(tk, true, true));
			tkAssert(tk.curType() == TkBracketR, "expected closing bracket");
			tk.next();            
		}

		block.func = parseFunction(tk, true, true, true);
		if (!templTypes.empty())
			block.func.templateTypes = templTypes;

		tkAssert(tk.curType() == TkCodeBlock && tk.isLast(), 
			"Malformed KERNEL, expected KERNEL(opts...) ret_type name(args...) { code }");

		block.line1 = tk.cur().line;
		processKernel(block, tk.cur().text, sink);
	}
	else if (kw == "PYTHON")
	{
		// template instantiation / alias
		if (tk.curType() == TkDescriptor && (tk.cur().text == "alias" || tk.cur().text == "instantiate"))  {
			string kw = tk.cur().text;
			Type aliasType;
			do {
				tk.next();
				aliasType = parseType(tk);
				processPythonInstantiation(block, aliasType, sink, inst);            
			} while (tk.curType() == TkComma);

			if (kw == "alias") {
				tkAssert(tk.curType() == TkDescriptor, "malformed preprocessor block. Expected 'PYTHON alias cname pyname;'");
				string aliasName = tk.cur().text;
				processPythonAlias(block, aliasType, aliasName, sink);
				tk.next();
			}
			tkAssert(tk.curType() == TkSemicolon && tk.isLast(), "malformed preprocessor block. Expected 'PYTHON alias/instantiate cname [pyname];'");
			sink.inplace << block.linebreaks();
				
			return;
		}
		List<Type> templTypes;

		// resolve template class
		Text templText;
		if (tk.curType() == TkTemplate) {
			TokenPointer t2(tk, &templText);
			t2.next();
			templTypes = parseTypeList(t2);            
		}

		// python class
		if (tk.curType() == TkClass && tk.cur().text != "typename") {
			block.cls = parseClass(tk);
			block.cls.templateTypes = templTypes;
			block.cls.prequel(&templText);
			tkAssert(tk.curType() == TkCodeBlock && tk.isLast(), "malformed preprocessor keyword block. Expected 'PYTHON class name : public X {}'");
			processPythonClass(block, tk.cur().text, sink, inst);
		}
		else // function or member
		{
			bool isConstructor = parent && tk.curType() == TkDescriptor && 
				 parent->name == tk.cur().text && tk.previewType() == TkBracketL;
			block.func = parseFunction(tk, false, !isConstructor, false);
			block.func.templateTypes = templTypes;
			block.func.prequel(&templText);
			
			if (isConstructor && tk.curType() == TkColon) {
				// read till end
				while(!tk.done() && tk.curType() != TkSemicolon && tk.curType() != TkCodeBlock) {
					block.initList += tk.cur().text;
					tk.next();
				}
				tkAssert(!tk.done(), "Constructor initializer list not limited");
			}

			if (tk.curType() == TkSemicolon && block.func.noParentheses) {
				tkAssert(tk.curType() == TkSemicolon && tk.isLast(), 
					"malformed preprocessor keyword block. Expected 'PYTHON type varname;'");
				processPythonVariable(block, sink);
			} else {
				tkAssert((tk.curType() == TkCodeBlock || tk.curType() == TkSemicolon) && tk.isLast(), 
					"malformed preprocessor keyword block. Expected 'PYTHON type funcname(args) [{}|;]'");
				processPythonFunction(block, tk.cur().text, sink, inst);
			}
		}

	}
}
