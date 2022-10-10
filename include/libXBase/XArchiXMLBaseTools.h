//Jean-Pierre Papelard/IGN/2006-2015/SR/MATIS
#ifndef _X_ARCHI_XML_BASE_TOOLS_H_
#define _X_ARCHI_XML_BASE_TOOLS_H_

#include <string>
#include "libXBase/XBase.h"
#include <tinyxml.h>

#ifndef WIN32
#ifndef _int64
typedef long long _int64;
#endif
#endif

class TiXmlNode;

namespace XArchiXML
{
	extern std::string NodeGeneaologie(TiXmlNode* node);
	extern std::string InitFileName(TiXmlNode* node);
	
	extern void AssertRoot(TiXmlNode* node, std::string rootname);

	extern TiXmlNode* FindNode(TiXmlNode* node, std::string nodename);
	extern TiXmlNode* FindAssertNode(TiXmlNode* node, std::string nodename);
	extern TiXmlNode* FindSubSubNode(TiXmlNode* node, std::string nodename);


	extern std::string ReadNodeAsString(TiXmlNode* node, std::string nodename, std::string valDef = std::string());
	extern std::string ReadAssertNodeAsString(TiXmlNode* node, std::string nodename);

	extern int ReadNodeAsInt(TiXmlNode* node, std::string nodename,int valDef = 0);
	extern int ReadAssertNodeAsInt(TiXmlNode* node, std::string nodename);


	extern bool ReadNodeAsBool(TiXmlNode* node, std::string nodename,bool valDef);
	extern bool ReadAssertNodeAsBool(TiXmlNode* node, std::string nodename);

	extern unsigned int ReadNodeAsUint32(TiXmlNode* node, std::string nodename,uint32 valDef = 0);
	extern unsigned int ReadAssertNodeAsUint32(TiXmlNode* node, std::string nodename);
	extern uint32 ReadAssertAttributeOrNodeAsUint32(TiXmlNode* node, std::string nodename );
    
    extern _int64 ReadAssertNodeAsInt64(TiXmlNode* node, std::string nodename);
    extern _int64 ReadNodeAsInt64(TiXmlNode* node, std::string nodename, _int64 valdef);

	extern double ReadNodeAsDouble(TiXmlNode* node, std::string nodename,double valDef = 0.);
	extern double ReadAssertNodeAsDouble(TiXmlNode* node, std::string nodename);

	extern float ReadAssertNodeAsFloat(TiXmlNode* node, std::string nodename);

	extern double ReadAssertAttributeOrNodeAsDouble(TiXmlNode* node, std::string nodename );
	extern double ReadAttributeOrNodeAsDouble(TiXmlNode* node, std::string nodename,double valdef );
	extern uint32 ReadAssertAttributeOrNodeAsUint32(TiXmlNode* node, std::string nodename );
	extern uint32 ReadAttributeOrNodeAsUint32(TiXmlNode* node, std::string nodename, uint32 valdef );
	extern int ReadAssertAttributeOrNodeAsInt(TiXmlNode* node, std::string nodename );
	extern int ReadAttributeOrNodeAsInt(TiXmlNode* node, std::string nodename,int valdef );
	extern std::string ReadAssertAttributeOrNodeAsString(TiXmlNode* node, std::string nodename );

	extern uint32 ReadArrayNode(TiXmlNode* node, std::string nodename, std::vector<std::string>* V);
	extern uint32 ReadArrayNodeAsInt(TiXmlNode* node, std::string nodename, std::vector<int>* V);
	extern uint32 ReadArrayNodeAsUint32(TiXmlNode* node, std::string nodename, std::vector<uint32>* V);
	extern uint32 ReadArrayNodeAsDouble(TiXmlNode* node, std::string nodename, std::vector<double>* V);

}

#endif
