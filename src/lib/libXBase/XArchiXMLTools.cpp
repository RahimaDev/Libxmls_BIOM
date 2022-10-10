#include "libXBase/XArchiXMLTools.h"
#include "libXBase/XArchiXMLBaseTools.h"
#include "libXBase/XRect.h"
#include "libXBase/XPt2D.h"
#include "libXBase/XPt3D.h"
//#include "libXSensor/XGrille.h"
#include "libXBase/XFrame.h"
#include "libXBase/XStringTools.h"
#include "libXBase/XColorTools.h"
#include "libXBase/XArchiXMLException.h"

#include "tinyxml.h"

namespace XArchiXML
{
//-------------------------------------------------------------------------
bool WriteXMLHeaderIso(std::ostream* out,std::string MainTag, std::string XslDtdFile)
{
	*out << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>" << std::endl;
	if(!XslDtdFile.empty())
		*out << "<?xml-stylesheet type=\"text/xsl\" href=\"" << XslDtdFile<< ".xsl\"?>" << std::endl;
	if(!MainTag.empty())
		*out << "<!DOCTYPE " << MainTag << " SYSTEM \"" << XslDtdFile <<".dtd\">"<< std::endl;
	return out->good();
}
//-----------------------------------------------------------------------------
// Conversion d'une chaine en chaine compatible avec XML
//-----------------------------------------------------------------------------
std::string OemToXml(std::string s)
{
	std::string result = "";
	char newStr[10];
	for (uint32 i = 0; i < s.length(); i++) {
		if (((byte)s[i] < 160) && (s[i] != '<') && (s[i] != '>'))
			result += s[i];
		else {
			sprintf(newStr, "&#%d;", (byte)s[i]);
			result += newStr;
		}
	}
	return result;
}
//-------------------------------------------------------------------------
bool XRect_LoadFromNode(XRect* rect, TiXmlNode* node)
{
	TiXmlNode* r = node->FirstChild("rect");
	if (r == NULL)
		return false;

	rect->X = ReadNodeAsUint32(r,"x",rect->X);
	rect->Y = ReadNodeAsUint32(r,"y",rect->Y);
	rect->W = ReadNodeAsUint32(r,"w",rect->W);
	rect->H = ReadNodeAsUint32(r,"h",rect->H);
	return true;
}
//-------------------------------------------------------------------------
XRect XRect_LoadFromNode(TiXmlNode* node)
{
	TiXmlNode* r = FindAssertNode(node,"rect");
	XRect rect;
	rect.X = ReadAssertNodeAsUint32(r,"x");
	rect.Y = ReadAssertNodeAsUint32(r,"y");
	rect.W = ReadAssertNodeAsUint32(r,"w");
	rect.H = ReadAssertNodeAsUint32(r,"h");
	return rect;
}

//-------------------------------------------------------------------------
void XFrame_LoadFromNode( XFrame* frame, TiXmlNode* node)
{
	frame->Xmin = ReadAssertNodeAsDouble(node,"xmin");
	frame->Xmax = ReadAssertNodeAsDouble(node,"xmax");
	frame->Ymin = ReadAssertNodeAsDouble(node,"ymin");
	frame->Ymax = ReadAssertNodeAsDouble(node,"ymax");
}

//-------------------------------------------------------------------------
bool XPt3D_LoadFromNode(XPt3D* pt3d, TiXmlNode* node)
{

	pt3d->X = ReadNodeAsDouble(node,"x",pt3d->X);
	pt3d->Y = ReadNodeAsDouble(node,"y",pt3d->Y);
	pt3d->Z = ReadNodeAsDouble(node,"z",pt3d->Z);
	return true;
}
//  ----------------------------------------------------------
XPt3D XPt3D_LoadFromNode(TiXmlNode* node)
{
	TiXmlNode* pt = FindAssertNode(node,"pt3d");
	return XPt3D(ReadAssertNodeAsDouble(pt,"x"),ReadAssertNodeAsDouble(pt,"y"),ReadAssertNodeAsDouble(pt,"z"));
}

//  ----------------------------------------------------------
XPt2D XPt2D_LoadFromNode(TiXmlNode* node)
{
	TiXmlNode* pt = FindAssertNode(node,"pt2d");
	return XPt2D(ReadAssertNodeAsDouble(pt,"x"),ReadAssertNodeAsDouble(pt,"y"));

}
//  ----------------------------------------------------------
// lecture de 2 doubles dans un enregistrement du type : <CDist>1042.79129219386891 1023.34369425995885</CDist>
XPt2D XPt2D_LoadLineInNode(TiXmlNode* node)
{
	TiXmlElement *  elt = node->ToElement ();
	const char* val = elt->GetText();
	if (val == NULL)
	{
		std::string ident = NodeGeneaologie(node) + "la balise est vide " ;
		throw(XmlException(ident.c_str(), XmlException::eBadData));
	}
	std::string test =val;
	XStringTools st;
	std::vector<std::string> vec = st.Tokenize(std::string(val),' ',true);
	if (vec.size() != 2)
	{
		std::string ident = NodeGeneaologie(node) + "la balise ne contient pas 2 valeurs " ;
		throw(XmlException(ident.c_str(), XmlException::eBadData));
	}
	XPt2D p;
	sscanf(vec[0].c_str(), "%lf", &p.X);
	sscanf(vec[1].c_str(), "%lf", &p.Y);
	return p;
}
//-------------------------------------------------------------------------
bool XPt2D_LoadFromNode(XPt2D* pt2d, TiXmlNode* node)
{
	TiXmlNode* pt = node->FirstChild("pt2d");
	if (pt == NULL)
		return false;

	pt2d->X = ReadNodeAsDouble(pt,"x",pt2d->X);
	pt2d->Y = ReadNodeAsDouble(pt,"y",pt2d->Y);
	return true;
}
//avec assertion   ----------------------------------------------------------
//-------------------------------------------------------------------------
//permet de lire au choix une matrice ou un quaternion dans un tag rotation
XMat3D XRotation_LoadFromNode(TiXmlNode* node)
{
	TiXmlNode* rot = FindAssertNode(node,"rotation");
	TiXmlNode* mat = rot->FirstChild("mat3d");
	if (mat != NULL){
		XMat3D mat = XMat3D_LoadFromNode(rot);
		//on verifie qu'il sagit bien d'une rotation
		//code a deplacer dans XMat3D
		XMat3D trn = mat.Trn();
		XMat3D diff = trn*mat;
		diff = diff-XMat3D(XPt3D(1,0,0),XPt3D(0,1,0),XPt3D(0,0,1));
		double max = diff.A.X;
		max = XMax(diff.A.Y, max);
		max = XMax(diff.A.Z, max);
		max = XMax(diff.B.X, max);
		max = XMax(diff.B.Y, max);
		max = XMax(diff.B.Z, max);
		max = XMax(diff.C.X, max);
		max = XMax(diff.C.Y, max);
		max = XMax(diff.C.Z, max);
		if( max > 0.0001){
			std::string ident = NodeGeneaologie(node) + "\\mat3d  n'est pas une matrice rotation (trans_R*R != id)" ;
			throw(XmlException(ident.c_str(), XmlException::eBadData));
		}
		return mat;
	}
	XQuaternion quat = XQuaternion_LoadFromNode(rot);
	quat.Normalise();
	return quat.GetRotationMatrix();
}
//-------------------------------------------------------------------------
bool XQuaternion_LoadFromNode(XQuaternion* q, TiXmlNode* node)
{
	q->x(ReadAssertNodeAsDouble(node,"x"));
	q->y(ReadAssertNodeAsDouble(node,"y"));
	q->z(ReadAssertNodeAsDouble(node,"z"));
	q->w(ReadAssertNodeAsDouble(node,"w"));
	return  true;
}
//-------------------------------------------------------------------------
XQuaternion XQuaternion_LoadFromNode(TiXmlNode* node)
{
	TiXmlNode* q = FindAssertNode(node,"quaternion");

	double x = ReadAssertNodeAsDouble(q,"x");
	double y = ReadAssertNodeAsDouble(q,"y");
	double z = ReadAssertNodeAsDouble(q,"z");
	double w = ReadAssertNodeAsDouble(q,"w");

	return   XQuaternion(x, y, z, w);
}


//-------------------------------------------------------------------------
XMat3D XMat3D_LoadFileXml(std::string filename)
{
	char message[1024];
	//validation du formalisme XML
	TiXmlDocument doc( filename.c_str() );
	if (!doc.LoadFile() )
	{
		sprintf( message,"Format XML non conforme pour %s. Error='%s'\n",filename.c_str(), doc.ErrorDesc() );
		throw(XmlException(message,XmlException::eBadFormat));
		return XMat3D();
	}
	TiXmlElement* root = doc.RootElement();
	std::string strRoot = std::string( root->Value() );
	if(strRoot == "mat3d")
	{
		return XMat3D_LoadSubNode(root);
	}
	std::string ident = filename + '\\' + strRoot;
	throw(XmlException(ident.c_str(),XmlException::eRacine));
}
//-------------------------------------------------------------------------
XMat3D XMat3D_LoadInFileXml(std::string filename)//charge la matrice a un sous niveau quelconque de l'arborescence
{
	char message[1024];
	//validation du formalisme XML
	TiXmlDocument doc( filename.c_str() );
	if (!doc.LoadFile() )
	{
		sprintf( message,"Format XML non conforme pour %s. Error='%s'\n",filename.c_str(), doc.ErrorDesc() );
		throw(XmlException(message,XmlException::eBadFormat));
		return XMat3D();
	}
	TiXmlNode* child = FindSubSubNode(doc.RootElement(),"mat3d");
	if(child == NULL){
		throw(XmlException("mat3d", XmlException::eTagObligatoire));
        return XMat3D() ;
	}

	return XMat3D_LoadSubNode(child);
}

//-------------------------------------------------------------------------
XMat3D XMat3D_LoadFromNode(TiXmlNode* node)
{
	TiXmlNode* mat = FindAssertNode(node,"mat3d");
	return XMat3D_LoadSubNode(mat);
}
//-------------------------------------------------------------------------
XMat3D XMat3D_LoadSubNode(TiXmlNode* mat3Dnode)
{
	TiXmlNode* l1 = FindAssertNode(mat3Dnode,"l1");
	XPt3D A = XPt3D_LoadFromNode(l1);
	TiXmlNode* l2 = FindAssertNode(mat3Dnode,"l2");
	XPt3D B = XPt3D_LoadFromNode(l2);
	TiXmlNode* l3 = FindAssertNode(mat3Dnode,"l3");
	XPt3D C = XPt3D_LoadFromNode(l3);

	return XMat3D(A,B,C);
}
//-------------------------------------------------------------------------

/*bool XGrille_LoadFromNode(XGrille* grille, TiXmlNode* node)
{
	TiXmlNode* grid = FindAssertNode(node,"grid");
	TiXmlNode* size = FindAssertNode(grid,"size");
	TiXmlNode* origine = FindAssertNode(grid,"origine");
        FindAssertNode(grid,"step");
	TiXmlNode* filename = FindAssertNode(grid,"filename");
        FindAssertNode(filename,"x");
        FindAssertNode(filename,"y");

	grille->SetSize(ReadNodeAsUint32(size,"x",0),ReadNodeAsUint32(size,"y",0));
	grille->SetOrigine(ReadNodeAsDouble(origine,"x",0),ReadNodeAsDouble(origine,"y",0));
	grille->SetStep(ReadNodeAsDouble(grid,"step",0));
	grille->FilenameX(ReadNodeAsString(filename,"x"));
	grille->FilenameY(ReadNodeAsString(filename,"y"));

	return true;
}*/

}
