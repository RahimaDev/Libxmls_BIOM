//Bertrand Cannelle/IGN/SR/MATIS
#ifndef _X_ARCHI_GEOREF_XML_H_
#define _X_ARCHI_GEOREF_XML_H_

#include <iostream>
class XArchiGeoref;
class TiXmlNode;

namespace XArchiXML
{
    void XArchiGeoref_LoadFromNode(XArchiGeoref* georef, TiXmlNode* node);

    bool XArchiGeoref_LoadFromNode(XArchiGeoref* georef, std::string filename);
}
#endif
