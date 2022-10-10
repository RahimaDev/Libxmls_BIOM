//Bertrand Cannelle/IGN/SR/MATIS
#ifndef _X_ARCHI_GEOREF_H_
#define _X_ARCHI_GEOREF_H_

#include "libXBase/XPt3D.h"
#include "libXBase/XMat3D.h"

//Donnees relatives au georeferencement
class XArchiGeoref
{
protected ://donnees membres
	XPt3D			m_Translation;	// Rattachement dans le referentiel vehicule
    //double       m_dEchelle; //inutilise pour l'instant !

	//donnees calculees
	//Les donnees sont exprimees dans le sens : repere capteur --> repere vehicule
	XMat3D m_Rotation;
public :
    XArchiGeoref(XPt3D T = XPt3D(), XMat3D M= XMat3D::Identite()):
        m_Translation(T),
        m_Rotation(M){}
    XArchiGeoref(const XArchiGeoref & g0, const XArchiGeoref & g1, double alpha)
    {
        double uma=1-alpha;
        m_Translation = uma*g0.Translation() + alpha*g1.Translation();
        m_Rotation = uma*g0.Rotation() + alpha*g1.Rotation();
        m_Rotation.A.Normalise();
        m_Rotation.B.Normalise();
        m_Rotation.C.Normalise();
    }

    XPt3D & Translation(){return m_Translation;}
    const XPt3D & Translation() const {return m_Translation;}
    void Translation(XPt3D p){m_Translation = p;}

    XMat3D & Rotation(){return m_Rotation;}
    const XMat3D & Rotation() const {return m_Rotation;}
	void Rotation(XMat3D M){m_Rotation = M;}

	XPt3D Applique_transfo(XPt3D Pt);//repere capteur --> repere vehicule
	XPt3D Applique_inverse_transfo(XPt3D Pt);//repere vehicule --> repere capteur

	XMat3D Applique_transfo(XMat3D M);
	XMat3D Applique_inverse_transfo(XMat3D M);

	XArchiGeoref Applique_transfo(XArchiGeoref G);
	XArchiGeoref Applique_inverse_transfo(XArchiGeoref G);

	std::string InfoTexte();
	bool XmlWrite(std::ostream* out);

};

#endif
