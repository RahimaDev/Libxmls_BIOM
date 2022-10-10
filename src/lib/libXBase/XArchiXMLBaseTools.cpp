#include "libXBase/XArchiXMLBaseTools.h"

#include <iostream>

#include "tinyxml.h"

#include "libXBase/XStringTools.h"
#include "libXBase/XArchiXMLException.h"


namespace XArchiXML
{
using namespace std;
//-------------------------------------------------------------------------
// STL independent ValueStr Recherche d'un noeud
//-------------------------------------------------------------------------
std::string ValueStr(TiXmlNode* node)
{
//#ifdef TIXML_USE_STL
//    return node->ValueStr();
//#else
    return std::string( node->Value() );
//#endif // TIXML_USE_STL
}

//-------------------------------------------------------------------------
// STL independent FirstChild
//-------------------------------------------------------------------------
TiXmlNode* FirstChild(TiXmlNode* node, std::string nodename)
{
#ifdef TIXML_USE_STL
    return node->FirstChild(nodename);
#else
    return node->FirstChild(nodename.c_str());
#endif // TIXML_USE_STL
}
//-------------------------------------------------------------------------
//Recherche d'un noeud obligatoire
//-------------------------------------------------------------------------
TiXmlNode* FindAssertNode(TiXmlNode* node, std::string nodename)
{
    TiXmlNode* child = FirstChild(node, nodename);
    if(child == NULL){
        std::string ident = NodeGeneaologie(node) + '\\' + nodename;
        throw(XmlException(ident.c_str(), XmlException::eTagObligatoire));
        return NULL ;
    }
    return child;
}
//-------------------------------------------------------------------------
//Recherche d'un noeud a n'importe quel niveau de l'arborescence
//-------------------------------------------------------------------------
TiXmlNode* FindSubSubNode(TiXmlNode* node, std::string nodename)
{
    cout << nodename << "==" << ValueStr(node) << endl;
    if(nodename.compare(ValueStr(node)) == 0) return node;
    TiXmlNode* child = NULL;
    while( (child = node->IterateChildren( child )) )
    {
        TiXmlNode* mat = FindSubSubNode(child, nodename);
        if(mat != NULL) return mat;
    }
    return NULL ;
}

//-------------------------------------------------------------------------
//Genealogie d'un noeud
//-------------------------------------------------------------------------
std::string NodeGeneaologie(TiXmlNode* node)
{
    std::string genealogie = ValueStr(node);
    TiXmlNode* parent = node->Parent();
    while(parent != NULL)
    {
        genealogie = ValueStr(parent) + '\\' + genealogie;
        parent = parent->Parent();
    }
    return genealogie;
}
//-------------------------------------------------------------------------
//Nom du fichier initial
//-------------------------------------------------------------------------
std::string InitFileName(TiXmlNode* node)
{
    std::string filename = ValueStr(node);
    TiXmlNode* parent = node->Parent();
    while(parent != NULL){
        filename = ValueStr(parent) ;
        parent = parent->Parent();
    }
    return filename;
}

//-------------------------------------------------------------------------
std::string ReadAssertNodeAsString(TiXmlNode* node, std::string nodename)
{
    TiXmlNode* n = FindAssertNode(node,nodename);
    const char* val = n->ToElement()->GetText();

    if (val == NULL){
        std::string ident = NodeGeneaologie(node) + '\\' + nodename;
        throw(XmlException(ident.c_str(), XmlException::eTagVide));
        return std::string();
    }

    return std::string(val);
}
//-------------------------------------------------------------------------
std::string ReadNodeAsString(TiXmlNode* node, std::string nodename, std::string valDef)
{
    TiXmlNode* n = FirstChild(node, nodename);
    if(n == NULL)
        return valDef;

    TiXmlElement *  elt = n->ToElement ();
    const char* val = elt->GetText();
    if (val == NULL)
        return valDef;
    return std::string(val);

}

//-----------------------------------------------------------------------------
int ReadNodeAsInt(TiXmlNode* node, std::string nodename,int valDef )
{
    std::string result = ReadNodeAsString(node,nodename);
    if(result.empty())
        return valDef;
    int n;
    sscanf(result.c_str(), "%d", &n);
    return n;
}
//-----------------------------------------------------------------------------
int ReadAssertNodeAsInt(TiXmlNode* node, std::string nodename )
{
    std::string strval = ReadAssertNodeAsString(node,nodename);
    int n;
    sscanf(strval.c_str(), "%d", &n);
    return n;
}
//-----------------------------------------------------------------------------
unsigned int ReadNodeAsUint32(TiXmlNode* node, std::string nodename,uint32 valDef)
{
    std::string result = ReadNodeAsString(node,nodename);
    if(result.empty())
        return valDef;
    uint32 n;
    sscanf(result.c_str(), "%u", &n);
    return n;
}
//-----------------------------------------------------------------------------
unsigned int ReadAssertNodeAsUint32(TiXmlNode* node, std::string nodename)
{
    std::string strval = ReadAssertNodeAsString(node,nodename);
    uint32 n;
    sscanf(strval.c_str(), "%u", &n);
    return n;
}
//-------------------------------------------------------------------------
uint32 ReadAssertAttributeOrNodeAsUint32(TiXmlNode* node, std::string nodename )
{
    TiXmlNode* n = node->FirstChild(nodename.c_str());
    if(n != NULL)
        return ReadAssertNodeAsUint32(node,nodename );

    int valInt = 0;
    TiXmlElement* elt =  node->ToElement();
    const char* strval = elt->Attribute(nodename.c_str(),&valInt);
    if(strval == NULL){
        std::string ident = NodeGeneaologie(node) + '\\' + nodename;
        throw(XmlException(ident.c_str(), XmlException::eTagVide));
    }
    if(valInt < 0){
        std::string ident = NodeGeneaologie(node) + '\\' + nodename;
        throw(XmlException(ident.c_str(), XmlException::eBadData));
    }
    return (uint32)valInt;
}
//-----------------------------------------------------------------------------
_int64 ReadAssertNodeAsInt64(TiXmlNode* node, std::string nodename)
{
    std::string strval = ReadAssertNodeAsString(node,nodename);
#ifdef WIN32
    return _atoi64(strval.c_str());
#else
    return atoll(strval.c_str());
#endif
}
//-----------------------------------------------------------------------------
_int64 ReadNodeAsInt64(TiXmlNode* node, std::string nodename, _int64 valdef)
{
	std::string strval = ReadNodeAsString(node,nodename);
	if(strval.empty())
		return valdef;
	_int64 n;
	sscanf(strval.c_str(), "%lld", &n);
	return n;
}
//-----------------------------------------------------------------------------
double ReadNodeAsDouble(TiXmlNode* node, std::string nodename, double valDef)
{
    std::string result = ReadNodeAsString(node,nodename);
    if(result.empty())
        return valDef;
    double x;
    sscanf(result.c_str(), "%lf", &x);
    return x;
}

//-----------------------------------------------------------------------------
double ReadAssertNodeAsDouble(TiXmlNode* node, std::string nodename)
{
    std::string strval = ReadAssertNodeAsString(node,nodename);
    double x;
    sscanf(strval.c_str(), "%lf", &x);
    return x;
}
//-----------------------------------------------------------------------------
float ReadAssertNodeAsFloat(TiXmlNode* node, std::string nodename)
{
    std::string result = ReadAssertNodeAsString(node,nodename);
    float x;
    sscanf(result.c_str(), "%f", &x);
    return x;
}

//-------------------------------------------------------------------------
bool ReadNodeAsBool(TiXmlNode* node, std::string nodename, bool valDef)
{
    std::string result = ReadNodeAsString(node,nodename);
    if(result.empty())
        return valDef;
    XStringTools ST;
    bool val;
    if(!ST.DecodeStringToBool(result,&val))
        return valDef;

    return val;
}
//-------------------------------------------------------------------------
bool ReadAssertNodeAsBool(TiXmlNode* node, std::string nodename)
{
    std::string strval = ReadAssertNodeAsString(node,nodename);
    XStringTools ST;
    bool val;
    if(!ST.DecodeStringToBool(strval,&val)){
        std::string ident = NodeGeneaologie(node) + '\\' + nodename;
        throw(XmlException(ident.c_str(), XmlException::eBadData));
        return false;
    }
    return val;
}
//-------------------------------------------------------------------------
void AssertRoot(TiXmlNode* node, std::string rootname)
{
    std::string strRoot = ValueStr(node);
    if(strRoot != rootname)
    {
        std::string ident = NodeGeneaologie(node) + '\\' + rootname;
        throw(XmlException(ident.c_str(),XmlException::eRacine));
    }
}
//-------------------------------------------------------------------------
double ReadAssertAttributeOrNodeAsDouble(TiXmlNode* node, std::string nodename )
{
    TiXmlNode* n = node->FirstChild(nodename.c_str());
    if(n != NULL)
        return ReadAssertNodeAsDouble(node,nodename );

    double valDouble = 0;
    TiXmlElement* elt =  node->ToElement();
    const char* strval = elt->Attribute(nodename.c_str(),&valDouble);
    if(strval == NULL){
        std::string ident = NodeGeneaologie(node) + '\\' + nodename;
        throw(XmlException(ident.c_str(), XmlException::eTagVide));
    }

    return valDouble;
}
//-------------------------------------------------------------------------
double ReadAttributeOrNodeAsDouble(TiXmlNode* node, std::string nodename,double valdef )
{
    TiXmlNode* n = node->FirstChild(nodename.c_str());
    if(n != NULL)
        return ReadNodeAsDouble(node,nodename,valdef );

    double valDouble = 0;
    TiXmlElement* elt =  node->ToElement();
    const char* strval = elt->Attribute(nodename.c_str(),&valDouble);
    if(strval == NULL)
        return valdef;

    return valDouble;
}

//-------------------------------------------------------------------------
uint32 ReadAttributeOrNodeAsUint32(TiXmlNode* node, std::string nodename, uint32 valdef )
{
    TiXmlNode* n = node->FirstChild(nodename.c_str());
    if(n != NULL)
        return ReadNodeAsUint32(node,nodename,valdef );

    int valInt = 0;
    TiXmlElement* elt =  node->ToElement();
    const char* strval = elt->Attribute(nodename.c_str(),&valInt);
    if(strval == NULL)
        return valdef;

    if(valInt < 0){
        std::string ident = NodeGeneaologie(node) + '\\' + nodename;
        throw(XmlException(ident.c_str(), XmlException::eBadData));
    }
    return (uint32)valInt;
}
//-------------------------------------------------------------------------
int ReadAssertAttributeOrNodeAsInt(TiXmlNode* node, std::string nodename )
{
    TiXmlNode* n = node->FirstChild(nodename.c_str());
    if(n != NULL)
        return ReadAssertNodeAsInt(node,nodename );

    int valInt = 0;
    TiXmlElement* elt =  node->ToElement();
    const char* strval = elt->Attribute(nodename.c_str(),&valInt);
    if(strval == NULL){
        std::string ident = NodeGeneaologie(node) + '\\' + nodename;
        throw(XmlException(ident.c_str(), XmlException::eTagVide));
    }
    return valInt;
}
//-------------------------------------------------------------------------
int ReadAttributeOrNodeAsInt(TiXmlNode* node, std::string nodename,int valdef )
{
    TiXmlNode* n = node->FirstChild(nodename.c_str());
    if(n != NULL)
        return ReadNodeAsInt(node,nodename,valdef );

    int valInt = 0;
    TiXmlElement* elt =  node->ToElement();
    const char* strval = elt->Attribute(nodename.c_str(),&valInt);
    if(strval == NULL)
        return valdef;
    return valInt;
}
//-------------------------------------------------------------------------
std::string ReadAssertAttributeOrNodeAsString(TiXmlNode* node, std::string nodename )
{
    TiXmlNode* n = node->FirstChild(nodename.c_str());
    if(n != NULL)
        return ReadAssertNodeAsString(node,nodename );

	TiXmlElement* elt =  node->ToElement();
        const char* strval = elt->Attribute(nodename.c_str());
	if(strval == NULL){
		std::string ident = NodeGeneaologie(node) + '\\' + nodename;
		throw(XmlException(ident.c_str(), XmlException::eTagVide));
	}

    return std::string(strval);
}
//-----------------------------------------------------------------------------
// Lecture d'un tableau de noeuds
//-----------------------------------------------------------------------------
uint32 ReadArrayNode(TiXmlNode* node, std::string nodename, std::vector<std::string>* V)
{
    TiXmlHandle hdl(node);
    uint32 count =0;
    TiXmlElement* champ =  hdl.FirstChild(nodename.c_str()).ToElement();
    while(champ)
    {
        count++;
        V->push_back(champ->GetText());
        champ = champ->NextSiblingElement(nodename.c_str());
    }
    return count;
}

//-----------------------------------------------------------------------------
// Lecture d'un tableau de noeuds avec conversion de type
//-----------------------------------------------------------------------------
uint32 ReadArrayNodeAsInt(TiXmlNode* node, std::string nodename, std::vector<int>* V)
{
    std::vector<std::string> T;
    ReadArrayNode(node,nodename, &T);
    int x;
    for (uint32 i = 0; i < T.size(); i++) {
        sscanf(T[i].c_str(), "%d", &x);
        V->push_back(x);
    }
    return V->size();
}
//-----------------------------------------------------------------------------
uint32 ReadArrayNodeAsUint32(TiXmlNode* node, std::string nodename, std::vector<uint32>* V)
{
    std::vector<std::string> T;
    ReadArrayNode(node,nodename, &T);
    uint32 x;
    for (uint32 i = 0; i < T.size(); i++) {
        sscanf(T[i].c_str(), "%u", &x);
        V->push_back(x);
    }
    return V->size();
}
//-----------------------------------------------------------------------------
uint32 ReadArrayNodeAsDouble(TiXmlNode* node, std::string nodename, std::vector<double>* V)
{
    std::vector<std::string> T;
    ReadArrayNode(node,nodename, &T);
    double x;
    for (uint32 i = 0; i < T.size(); i++) {
        sscanf(T[i].c_str(), "%lf", &x);
        V->push_back(x);
    }
    return V->size();
}

}
