//Jean-Pierre Papelard/IGN/2006-2015/SR/MATIS
#ifndef _X_ARCHI_XML_TOOLS_H_
#define _X_ARCHI_XML_TOOLS_H_

#include <string>

#include "libXBase/XQuaternion.h"

class XRect;
class XPt3D;
class XPt2D;

class XFrame;
class XDefect;
class XColorInfo;
//class XGrille;
class XMat3D;
class TiXmlNode;


namespace XArchiXML
{
	extern bool WriteXMLHeaderIso(std::ostream* out,std::string MainTag, std::string XslDtdFile);
	extern std::string OemToXml(std::string s);

	extern bool XRect_LoadFromNode(XRect* rect, TiXmlNode* node);
	extern XRect XRect_LoadFromNode(TiXmlNode* node);

	extern bool XPt3D_LoadFromNode(XPt3D* pt3d, TiXmlNode* node);
	extern XPt3D XPt3D_LoadFromNode(TiXmlNode* node);

	extern bool XPt2D_LoadFromNode(XPt2D* pt2d, TiXmlNode* node);
	extern XPt2D XPt2D_LoadFromNode(TiXmlNode* node);
	extern XPt2D XPt2D_LoadLineInNode(TiXmlNode* node);

	extern bool XQuaternion_LoadFromNode(XQuaternion* q, TiXmlNode* node);
	extern XQuaternion XQuaternion_LoadFromNode(TiXmlNode* node);

	extern XMat3D XMat3D_LoadFileXml(std::string filename);
	extern XMat3D XMat3D_LoadInFileXml(std::string filename);
	extern XMat3D XMat3D_LoadFromNode(TiXmlNode* node);
	extern XMat3D XMat3D_LoadSubNode(TiXmlNode* node);

	extern XMat3D XRotation_LoadFromNode(TiXmlNode* node);

	extern void XFrame_LoadFromNode( XFrame* frame, TiXmlNode* node);

    //extern bool XGrille_LoadFromNode(XGrille* grille, TiXmlNode* node);
}

#endif
