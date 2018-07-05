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

#ifndef _P_TOKENIZE_H
#define _P_TOKENIZE_H

#include <string>
#include <vector>
#include "code.h"

enum TokenType { TkNone = 0, TkComment, TkWhitespace, TkCodeBlock, TkDescriptor, TkComma, TkBracketL, 
				 TkBracketR, TkTBracketL, TkTBracketR, TkPointer, TkRef, TkDoubleColon, TkSemicolon, 
				 TkSimpleType, TkTypeQualifier, TkConst, TkEnd, TkManta, TkUnsupportedKW, TkClass,
				 TkInline, TkTemplate, TkStatic, TkVirtual, TkString, TkPublic, TkColon, TkOperator };

const std::string unsupportedKeywords[] = {"and", "and", "and_eq", "auto", "bitand", "bitor", "break", 
	"catch", "const_cast", "continue", "default", "delete", "do", "dynamic_cast", "else", "enum", 
	"explicit", "export", "extern", "for", "friend", "goto", "if", "mutable", "namespace", "new", 
	"not", "not_eq", "or", "or_eq", "private", "protected", "register", 
	"reinterpret_cast", "return", "sizeof", "static_cast", "switch", "this", "throw", "try", "typedef", 
	"union", "using", "volatile", "while", "xor", "xor_eq", "" };

inline bool isIntegral(const std::string& t) {
	return t=="int" || t=="char" || t=="unsigned" || t=="bool" || t=="float" || t=="long" || 
		   t=="double" || t=="Real" || t=="Vec3" || t=="Vec3i" || t=="string" || t=="std::string" ||
		   t=="PbType" || t=="PbTypeVec";
}

struct Token {
	Token(TokenType t, int l) : type(t), text(""), line(l) {}
	Token(TokenType t, int l, char c) : type(t), text(""), line(l) { text += c; }
	TokenType type;
	std::string text;
	int line;
};

// tracks a set of tokens, and the current position in this list
struct TokenPointer {
	TokenPointer(const std::vector<Token>& t, Text *track) : parent(0),queue(t),ptr(0),txt(track) { reset(); }
	TokenPointer(TokenPointer& t, Text *track) : parent(&t),queue(t.queue),ptr(t.ptr),txt(track) { reset(); }
	TokenPointer *parent;
	const std::vector<Token>& queue;
	int ptr;
	Text *txt;

	inline void reset() { txt->reset(); consumeWhitespace(); if(!done()) txt->line0 = cur().line;  }
	inline TokenType curType() { return done() ? TkEnd : cur().type; }
	TokenType previewType();
	inline const Token& cur() { return queue[ptr]; }
	inline bool done() { return ptr >= (int)queue.size(); }
	inline bool isLast() { return ptr == (int)queue.size()-1;}
	void forward(const std::string& minimal, const std::string& original, int offset);
	void next();
	void consumeWhitespace();
	void errorMsg(const std::string& msg);
	std::string backtrace();
};

#endif