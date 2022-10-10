#include "libXBase/XArchiGeorefXML.h"
#include "libXBase/XArchiGeoref.h"
#include "libXBase/XArchiXMLBaseTools.h"
#include "libXBase/XArchiXMLTools.h"
#include "libXBase/XArchiXMLException.h"

#include "tinyxml.h"

namespace XArchiXML
{
//-------------------------------------------------------------------------
void XArchiGeoref_LoadFromNode(XArchiGeoref* georef, TiXmlNode* node)
{
    georef->Translation(XPt3D_LoadFromNode(node));
    georef->Rotation(XRotation_LoadFromNode(node));
    return;
}

bool XArchiGeoref_LoadFromNode(XArchiGeoref* georef, std::string filename)
{
    try
    {
        //validation du formalisme XML
        TiXmlDocument doc( filename.c_str() );
        if (!doc.LoadFile())
            std::cout << "ERROR: Invalid XML format in " << filename << " Error=" << doc.ErrorDesc() << std::endl;
        XArchiGeoref_LoadFromNode(georef, XArchiXML::FindSubSubNode(doc.RootElement(), "georef"));
    }
    catch(XArchiXML::XmlException e)
    {
        std::cout << "Reading " << filename << " raised " << e.Erreur() << std::endl;
        return false;
    }
    return true;
}

}
