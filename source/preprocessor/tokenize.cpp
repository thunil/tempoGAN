/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011-2014 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Preprocessor Tokenizer
 *
 ******************************************************************************/

#include <string>
#include <iostream>
#include <vector>
#include <algorithm>

#include "prep.h"

using namespace std;

//*************************************************************************************
// TokenPointer class members

TokenType TokenPointer::previewType() {
	for (int i=ptr+1; i < (int)queue.size(); i++) {
		if (queue[i].type != TkWhitespace && queue[i].type != TkComment)
			return queue[i].type;
	}
	return TkNone;
}

void TokenPointer::consumeWhitespace() {
	if (done()) return;
	while (cur().type == TkWhitespace || cur().type == TkComment) {
		string minimal = (!txt->minimal.empty() && *(txt->minimal).rbegin() == ' ') ? "" : " ";
		forward(minimal,cur().text,1);
		if (ptr >= (int)queue.size())
			errMsg(-1, "Preprocessor ran out of tokens. This shouldn't happen.");
	}
}

void TokenPointer::next() {
	forward(cur().text, cur().text, 1);
	consumeWhitespace();
}

void TokenPointer::forward(const string& minimal, const string& original, int offset) {
	txt->minimal += minimal; 
	txt->original += original;
	ptr += offset;
	if (parent)
		parent->forward(minimal, original, offset);
}

string TokenPointer::backtrace() {
	string bt = "[" + txt->dynamicClass() + "]: " + txt->original + "\n";
	if (parent)
		bt += parent->backtrace();
	return bt;
}

void TokenPointer::errorMsg(const string& msg) {
	string emsg = "'" + txt->original;
	if (!done()) emsg += cur().text;
	emsg += "' " + msg + "\n";
	emsg += " Preprocessor backtrace:\n" + backtrace(); 

	int line = -1;
	if (!queue.empty()) 
		line = done() ? queue.back().line : cur().line; 
	
	errMsg(line, emsg);
}

//*************************************************************************************
// Lexing functions

// tokenize and parse until keyword section ends
void tokenizeBlock(vector<Token>& tokens, const string& kw, const string& text, size_t& i, int& line, bool& hasBrackets) {
	tokens.push_back(Token(TkWhitespace, line));
	BracketStack brackets;
		
	// tokenize loop
	bool comment=false, slComment=false, define=false, extendLine=false;
	bool isString=false;
	bool foundChars=false;
	int codeblockLevel = 0;
	for (; i<text.size(); i++) {
		char c = text[i];
		if (c=='\r') continue;

		string& curstr = tokens.back().text;
		bool lastchar = (i == text.size()-1) || (i == text.size()-2 && text[i+1]=='\r');
		
		// track lines
		bool isEOL = !extendLine && c=='\n';
		if (c=='\\') 
			extendLine = true;
		else if (c!='\n') 
			extendLine = false;
		if (c=='\n') line++;
		
		// track comments and defines
		if (comment) {
			if (!lastchar && c=='*' && text[i+1]=='/') comment = false;            
			curstr += c;
		}
		else if (slComment) {
			if (isEOL) slComment = false;
			curstr += c;
		} 
		else if (define) {
			if (isEOL) define = false;
			curstr += c;
		}
		else if (!lastchar && c=='/' && text[i+1]=='*') {
			comment = true;
			if (codeblockLevel==0) 
				tokens.push_back(Token(TkComment, line, c));
			else
				curstr += c;
		}
		else if (!lastchar && c=='/' && text[i+1]=='/') {
			slComment = true;
			if (codeblockLevel==0) 
				tokens.push_back(Token(TkComment, line, c));
			else
				curstr += c;
		}
		else if (!isString && c=='\"' && codeblockLevel == 0) {
			isString = true;
			tokens.push_back(Token(TkString, line, c));
		}
		else if (isString) {
			if (c=='\"') isString = false;
			curstr += c;
		}        
		else if (c=='#') {
			define = true;
			if (codeblockLevel==0) 
				tokens.push_back(Token(TkComment, line, c));
			else
				curstr += c;
		}   
		
		// track codeblock
		else if (codeblockLevel > 0) {
			curstr += c;
			if (c == '}') codeblockLevel--;
			if (c == '{') codeblockLevel++;
			
			if (codeblockLevel == 0) {
				// block finished, return to parseText
				return;
			}
		}
		else if (c=='{') {
			codeblockLevel++;
			if (!brackets.empty())
				errMsg(line, "codeblock {} is not allowed inside brackets.");                
			tokens.push_back(Token(TkCodeBlock, line, c));
		}
				
		// track brackets
		else if (c=='(') {
			tokens.push_back(Token(TkBracketL, line, c));
			brackets.push_back(c);
			// detect first parentheses are keyword
			if(!foundChars) hasBrackets = true;
		}
		else if (c==')') {
			if (brackets.pop() != '(') 
				errMsg(line, "bracket mismatch, can't close ')'");
			tokens.push_back(Token(TkBracketR, line, c));
		}
		else if (c=='<') {
			tokens.push_back(Token(TkTBracketL, line, c));
			brackets.push_back(c);
		}
		else if (c=='>') {
			if (brackets.pop() != '<')
				errMsg(line, "bracket mismatch, can't close '>'");
			tokens.push_back(Token(TkTBracketR, line, c));
		}
		
		// track symbol tokens
		else if (c==' ' || c=='\t' || c=='\n') {
			tokens.push_back(Token(TkWhitespace, line, c));
		}
		else if (c==',') {
			tokens.push_back(Token(TkComma, line, c));
		}
		else if (c=='*' || c=='-' || c=='+' || c=='/' || c=='=') {
			tokens.push_back(Token(TkOperator, line, c));
		}
		else if (c=='&') {
			tokens.push_back(Token(TkRef, line, c));
		}
		else if (c == ':') {
			if (!lastchar && text[i+1]==':') {
				tokens.push_back(Token(TkDoubleColon, line));
				tokens.back().text = "::";
				i++;
			} else {
				tokens.push_back(Token(TkColon, line, c));
			}
		}
		else if (c==';') {
			if (!brackets.empty())
				errMsg(line, "Semicolon is not allowed inside brackets.");                
			
			tokens.push_back(Token(TkSemicolon, line, c));
						
			// block finished, return to parseText
			return;
		}
		
		// track descriptors
		else if (isNameChar(c)) {
			foundChars=true;
			if (tokens.back().type != TkDescriptor)
				tokens.push_back(Token(TkDescriptor, line, c));
			else
				curstr += c;
		}
		// misc stuff: don't bother with more complex tokens
		else {
			if (tokens.back().type != TkNone) {
				tokens.push_back(Token(TkNone, line, c));
			} else
				curstr += c;
		}        
	}
	
	// EOF before block end ?
	errMsg(line, "EOF before block for keyword '" + kw + "' was closed.");
	
	return;
}

void convertKeywords(vector<Token>& tokens) {
	for (size_t i=0; i<tokens.size(); i++) {
		if (tokens[i].type == TkDescriptor) {
			const string& t = tokens[i].text;
			if (t == "const")
				tokens[i].type = TkConst;
			else if (t == "unsigned" || t == "signed")
				tokens[i].type = TkTypeQualifier;
			else if (t == "char" || t == "int" || t == "short" || t == "long" ||
					 t == "float" || t == "double" || t == "string" || t == "bool" ||
					 t == "Vec3" || t == "Vec3i")
				tokens[i].type = TkSimpleType;   
			else if (t == "PYTHON" || t == "KERNEL")
				tokens[i].type = TkManta;
			else if (t == "struct" || t == "class" || t == "typename")
				tokens[i].type = TkClass;
			else if (t == "inline")
				tokens[i].type = TkInline;
			else if (t == "public")
				tokens[i].type = TkPublic;
			else if (t == "virtual")
				tokens[i].type = TkVirtual;
			else if (t == "static")
				tokens[i].type = TkStatic;
			else if (t == "template")
				tokens[i].type = TkTemplate;
			else {
				for (int i=0; !unsupportedKeywords[i].empty(); i++) {
					if (t == unsupportedKeywords[i])
						tokens[i].type = TkUnsupportedKW;
				}
			}
		}
	}
}

// parse complete file.
// Defer parts to be processed to processTextBlock, and directly copy the rest
void processText(const string& text, int baseline, Sink& sink, const Class* parent, vector<Instantiation>& inst) {
	ostream& newText = sink.inplace;
	
	// no real lexing yet, only track define and comment blocks
	string word = "";
	int line=baseline;
	bool comment=false, slComment=false, define=false, extendLine=false, isString=false;
	for (size_t i=0; i<text.size()-1; i++) {
		char c = text[i];
		if (c == '\r') continue;
		
		// track lines
		bool isEOL = !extendLine && c=='\n';
		if (c=='\\') 
			extendLine = true;
		else if (c!='\n') 
			extendLine = false;
		if (c=='\n') line++;
		
		// track comments and defines
		if (comment) {
			if (c=='*' && text[i+1]=='/') comment = false;            
			newText << c;
		} 
		else if (slComment) {
			if (isEOL) slComment = false;
			newText << c;
		} 
		else if (define) {
			if (isEOL) define = false;
			newText << c;
		}
		else if (c=='/' && text[i+1]=='*') {
			comment = true;
			newText << word;
			word = "";
			newText << c;
		}
		else if (c=='/' && text[i+1]=='/') {
			slComment = true;
			newText << word;
			word = "";
			newText << c;
		}
		else if (!isString && c=='\"') {
			isString = true;
			newText << word;
			word = "";
			newText << c;
		}
		else if (isString) {
			if (c=='\"') isString = false;
			newText << c;
		}
		else if (c=='#') {
			define = true;
			newText << word;
			word = "";
			newText << c;            
		} 
		else {
			// tokenize keywords
			if (isNameChar(c))
				word += c;
			else {
				if (word == "KERNEL" || word == "PYTHON") {
					vector<Token> tokens;
					bool brackets = false;
					tokenizeBlock(tokens, word, text, i, line, brackets);
					if(!brackets) {
						errMsg(line, "KERNEL and PYTHON keywords must have \"()\" ");
					}
					convertKeywords(tokens);
					parseBlock(word, tokens, parent, sink, inst); 
				} else {
					newText << word;                    
					newText << c; 
				}
				word = "";         
			}
		}
	}
	newText << word;
}


