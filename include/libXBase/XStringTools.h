//Francois Becispahic/IGN/2000-2003/Projet Camera Numerique
// J.P.Papelard - SR MATIS

#ifndef __XSTRINGTOOLS_H
#define __XSTRINGTOOLS_H

#include <cstdio>
#include "libXBase/XBase.h"
#include <string>
#include <vector>

class XStringTools 
{
public:
	std::string Format(uint16 sizeOut, char defaultChar, uint32 valToFormat);
	std::string Fusion(std::string str1, std::string str2);
	std::string Uppercase(std::string word);
	std::string lowercase(std::string word);
	bool equal(const std::string& s1, const std::string& s2);
	uint32 StringToUint32(std::string& s);
	uint32 ExtractInteger(std::string& s ,bool front = true);
	uint32 ExtractIntegerBack(std::string& s);
	uint32 ExtractIntegerBack(std::string& s, std::string::size_type* posIntegerBegin, std::string::size_type* posIntegerEnd);
	std::string RemoveBlanks(std::string s);
	std::string RemoveAllBlanks(std::string s);
	bool AsBlankInside(std::string s);//renvoie true si la chaine contient un espace
	bool DecodeStringToBool(std::string strvalIn, bool* valOut);
	std::vector<std::string> Tokenize(const std::string &str, char sep, bool multiSepAsOne);
	void Tokenize(const std::string& str,std::vector<std::string>& tokens,const std::string& delimiters = " ");
	std::string Substr(char* ini, char* fin);
	std::string ExtractStringBetween(std::string& s,std::string& beginAfter, std::string& StopBefore);
	bool IdentifyLineAsComment(std::string& ligne);
	uint16 numberOfDigit(uint32 input);
	std::string itoa(int i);
	std::string CompatibleFilename(std::string s);
};

#endif //__XSTRINGTOOLS_H
