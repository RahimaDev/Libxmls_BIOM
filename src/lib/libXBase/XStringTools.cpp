#include "libXBase/XStringTools.h"
#include <sstream>

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
std::vector<std::string> XStringTools::Tokenize(const std::string &str, char sep, bool multiSepAsOne)
{
	std::vector<std::string> vec;
	std::string token;
	std::istringstream is(str);
	while (getline(is, token, sep)){
		if(multiSepAsOne){//on considere plusieurs separateurs consecutifs comme un seul
			if(!token.length())
				continue;
		}
		vec.push_back(token);
	}

return vec;
}
//-----------------------------------------------------------------------------
void XStringTools::Tokenize(const std::string& str,std::vector<std::string>& tokens,const std::string& delimiters)
{
    // Skip delimiters at beginning.
    std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    std::string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (std::string::npos != pos || std::string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}


//-----------------------------------------------------------------------------
std::string XStringTools::Substr(char* ini, char* fin)
{
	std::string extract;
	uint32 count = (uint32)(fin-ini);
	char* temp = new char[count+1];
	memcpy(temp,ini ,count);
	temp[count]='\0';
	extract = temp;
	delete temp;
	return extract;
}

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
std::string XStringTools::Format(uint16 size, char defaultChar, uint32 valToFormat)
{
	std::string result;
	char temp[1024];
	sprintf(temp,"%d",valToFormat);
	std::string strVal = temp;
	// OTO : decalaration de la variable de boucle, sinon, probleme de scope sous linux
	int i;
	for(i =0; i< size; i++){
		result += defaultChar;
	}
	uint16 rang = 0;
	for(i = (int)(size - strVal.size()); i<(int) size; i++){
		result[i] = strVal[rang];
		rang++;
	}
	return result;
}
//-----------------------------------------------------------------------------
// str1= ABC000  str2= 12		result= ABC012
// str1= ABC000  str2= 1234567  result= 1ABC000
// str1= ABC000  str2= 123456	result= ABC000
//-----------------------------------------------------------------------------
std::string XStringTools::Fusion(std::string str1, std::string str2)
{
	std::string ref = str1;
	uint32 sizeRef = str1.size();

	std::string fus = str2;
	uint32 sizeFus = str2.size();

	if(sizeRef == sizeFus)
		return str2;

	if(sizeRef < sizeFus){
		ref = str2;
		sizeRef = str2.size();
		fus = str1;
		sizeFus = str1.size();
	}
	for(uint32 i = 0; i< sizeFus; i++)
		ref[sizeRef-sizeFus+i] = fus[i];
	
	return ref;
}
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
std::string XStringTools::Uppercase(std::string word)
{
	if(word.empty())
		return word;

	std::string ret_str = "";
	for (unsigned int index = 0; index < word.size(); index++)
		ret_str += ::toupper(word[index]);

	return ret_str;
}
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
std::string XStringTools::lowercase(std::string word)
{
	if(word.empty())
		return word;

	std::string ret_str = "";

	for (unsigned int index = 0; index < word.size(); index++)
		ret_str += ::tolower(word[index]);

  return ret_str;
}
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
bool XStringTools::equal(const std::string& s1, const std::string& s2)
{
	unsigned int count =(unsigned int)( (s2.size() < s1.size() ? s2.size(): s1.size()) );
  std::string str1 = lowercase(s1);
  std::string str2 = lowercase(s2);

  for(unsigned int i=0;i<count;i++)
    if(str1[i] != str2[i])
      return false;
  return true;
}
//-----------------------------------------------------------------------------
//Decode un booleen dans une chaine de caractere
//-----------------------------------------------------------------------------
bool XStringTools::DecodeStringToBool(std::string strvalIn,bool* valOut)
{
	//TRUE FALSE True False true false
	if(strvalIn.size() >= 4){
		if(equal(strvalIn,"true")){
			*valOut = true;
			return true;
		}
		if(equal(strvalIn,"false")){
			*valOut = false;
			return true;
		}
		return false;
	}

	//0 1 00 01 000 001
	int n;
	sscanf(strvalIn.c_str(), "%d", &n);
	if(n==1){
		*valOut = true;
		return true;
	}
	if(n==0){
		*valOut = false;
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
//Convertit une string en uint32 en concatenant tous les caracteres numeriques contenus dans la chaîne
//-----------------------------------------------------------------------------
uint32 XStringTools::StringToUint32(std::string& s)
{
	std::string strInteger = "0";
	for(uint32 i=0; i<s.size(); i++)
		if( isdigit((char)s[i]) )
			strInteger += s[i];
	
	uint32 integer = (uint32)atof(strInteger.c_str());
	return integer;
}
//-----------------------------------------------------------------------------
//Extrait une valeur entiere dans une string.a partir du debut
//-----------------------------------------------------------------------------
uint32 XStringTools::ExtractInteger(std::string& s, bool front)
{
	if(!front)
		return ExtractIntegerBack(s);

	std::string::size_type pos1 = std::string::npos;
	std::string::size_type count = std::string::npos;
	//recherche du premier digit
	std::string::size_type i;
	for(i=0; i<s.size(); i++){
		if( isdigit((char)s[i]) ){
			pos1 = i;
			break;
		}
	}
	if(pos1 == std::string::npos)
		return 0;
	count = 0;
	//recherche du dernier digit
	for(i= pos1; i<s.size(); i++){
		if(!isdigit((char)s[i]) )
			break;

		count++;
	}
	std::string strInteger = s.substr(pos1,count);
	uint32 integer = (uint32)atof(strInteger.c_str());

	return integer;
}
//-----------------------------------------------------------------------------
//Extrait une valeur entiere dans une string.a partir de la fin
//-----------------------------------------------------------------------------
uint32 XStringTools::ExtractIntegerBack(std::string& s)
{
	std::string::size_type pos2 = std::string::npos;
	std::string::size_type count = std::string::npos;

	//recherche du dernier digit
	int i;
	for(i= (int)s.size(); i>0; i--){
		if( isdigit((char)s[i-1]) ){
			pos2 = i;
			break;
		}
	}
	if(pos2 == std::string::npos)
		return 0;

	count = 0;
	//recherche du premier digit
	for(i= (int)pos2; i>0; i--){
		if(!isdigit((char)s[i-1]) )
			break;

		count++;
	}
	std::string strInteger = s.substr(pos2-count,count);
	uint32 integer = (uint32)atof(strInteger.c_str());

	return integer;
}
//-----------------------------------------------------------------------------
//Extrait une valeur entiere dans une string.a partir de la fin
//-----------------------------------------------------------------------------
uint32 XStringTools::ExtractIntegerBack(std::string& s, std::string::size_type* posIntegerBegin, std::string::size_type* posIntegerEnd)
{
	std::string::size_type posEnd = std::string::npos;
	std::string::size_type count = std::string::npos;

	//recherche du dernier digit
	int i;
	for(i= (int)s.size(); i>0; i--){
		if( isdigit((char)s[i-1]) ){
			posEnd = i;
			break;
		}
	}
	if(posEnd == std::string::npos)
		return 0;

	count = 0;
	//recherche du premier digit
	for(i= (int)posEnd; i>0; i--){
		if(!isdigit((char)s[i-1]) )
			break;

		count++;
	}
	std::string::size_type posBegin = posEnd-count;
	std::string strInteger = s.substr(posBegin,count);
	uint32 integer = (uint32)atof(strInteger.c_str());

	if(posIntegerBegin != NULL)
		*posIntegerBegin = posBegin;
	if(posIntegerEnd != NULL)
		*posIntegerEnd = posEnd;

	return integer;
}

//-----------------------------------------------------------------------------
//Retire les blancs en debut et fin de de chaine
//-----------------------------------------------------------------------------
std::string XStringTools::RemoveBlanks(std::string s)
{
  for(unsigned int i=0; i<s.size(); i++)
  {
    if(s[i] != ' ' && s[i] != '\t')//premier caractere non blanc
	{
	for(unsigned int j=(int)(s.size())-1;j>=i;j--)
	  {
		if(s[j] != ' ' && s[j] != '\t')//dernier caractere non blanc
		{
			return s.substr(i,j+1-i);
		}
      }
    }
  }
  return "";
}
//-----------------------------------------------------------------------------
bool XStringTools::AsBlankInside(std::string s)
{
  for(unsigned int i=0; i<s.size(); i++)
    if(s[i] == ' ' )
		return true;
	
  return false;
}
//-----------------------------------------------------------------------------
//Retire tous les blancs
//-----------------------------------------------------------------------------
std::string XStringTools::RemoveAllBlanks(std::string s)
{
  std::string tmp = "";

  for(unsigned int i=0;i<s.size();i++){
    if(s[i] != ' ' && s[i] != '\t' && s[i] != '\n')
      tmp += s[i];
  }
  return tmp;
}
//-----------------------------------------------------------------------------
//Extrait une chaine de caractere entre 2 chaines contenues
//-----------------------------------------------------------------------------
std::string XStringTools::ExtractStringBetween(std::string& s,std::string& beginAfter, std::string& StopBefore)
{
	std::string::size_type pos1 = s.find(beginAfter);
	if(pos1 == std::string::npos)
		return std::string();
	pos1 += beginAfter.size();
	std::string::size_type pos2 = s.rfind(StopBefore);
	if(pos2 == std::string::npos)
		return std::string();
	return s.substr(pos1 ,pos2-pos1);
}
//-----------------------------------------------------------------------------
//methode commune a differents fichiers texte permettant de considerer une ligne comme un commentaire
//une ligne de commentaire commence par une * ( comp3D par exemple)
//un commentaire comence par // (entete geoconcept par exemple)
//si le premier caractere rencontre est different de " \t*" la ligne est valide
bool XStringTools::IdentifyLineAsComment(std::string& ligne)
{
	//on cherche le premier caractere different de ' 'et '\t' 
  for(unsigned int i=0; i<ligne.size();i++)
  {
    if(ligne[i] == ' ' || ligne[i] == '\t' )
		continue;
	//test du premier caractere valide
    if(ligne[i] == '*' || ligne[i] == '/' )
		return true;
	 return false;
  }
  return true;// aucun caractere valide dans la ligne
}
//-----------------------------------------------------------------------------
//nombre de digit pour un nombre
//-----------------------------------------------------------------------------
uint16 XStringTools::numberOfDigit(uint32 input)
{
     uint16 numOfDigit = 0;
     while (input / 10 >= 1) {
           numOfDigit++;
           input /= 10;
      }
      return numOfDigit + 1;
}
//-----------------------------------------------------------------------------
std::string XStringTools::itoa(int i)
{
	char str[1024];
	sprintf(str,"%d",i);
	return std::string(str);
}
//-----------------------------------------------------------------------------
//rend un chaine compatible avec un nom de fichier
//-----------------------------------------------------------------------------
std::string XStringTools::CompatibleFilename(std::string s)
{
  std::string tmp = "";
  
  for(unsigned int i=0; i<s.size(); i++){
    if(    s[i] == '\\' 
		|| s[i] == '/' 
		|| s[i] == '\t' 
		|| s[i] == '\n'
		|| s[i] == ':'
		|| s[i] == '?'
		|| s[i] == '"'
		|| s[i] == '<'
		|| s[i] == '>'
		|| s[i] == '|'
		|| s[i] == '*')
		tmp += '_';
	else
      tmp += s[i];
  }
  return tmp;
}
//-----------------------------------------------------------------------------
/*
std::string lowercase(std::string word){
#ifdef NOCASE
  std::string ret_str = "";

  for (unsigned int index = 0; index < word.size(); index++)
    ret_str += ::tolower(word[index]);

  return ret_str;
#else
  return word;
#endif
}

std::string uppercase(std::string word){
#ifdef NOCASE
  std::string ret_str = "";

  for (unsigned int index = 0; index < word.size(); index++)
    ret_str += ::toupper(word[index]);

  return ret_str;
#else
  return word;
#endif
}

	*/
