#include "libXBase/XArchiGeoref.h"
#include "libXBase/XQuaternion.h"

#include <sstream>

//-----------------------------------------------------------------------------
//resultat = eche* Rot*Pt+Trans  //repere capteur --> repere vehicule
XPt3D XArchiGeoref::Applique_transfo(XPt3D Pt)
{
    return m_Rotation*Pt + m_Translation; // m_dEchelle *
}
//-----------------------------------------------------------------------------
//resultat = transpose(Rot)*(Pt-Trans)/eche //repere vehicule --> repere capteur
XPt3D XArchiGeoref::Applique_inverse_transfo(XPt3D Pt)
{
    return (m_Rotation.Trn()*(Pt-m_Translation)); // 1/m_dEchelle *
}
//-----------------------------------------------------------------------------
XMat3D XArchiGeoref::Applique_transfo(XMat3D M)
{
    return m_Rotation * M;
}
//-----------------------------------------------------------------------------
XMat3D XArchiGeoref::Applique_inverse_transfo(XMat3D M)
{
	XMat3D R = m_Rotation * M;
    return R.Trn(); // BV: should be return m_Rotation.Trn() * M
}

//-----------------------------------------------------------------------------
XArchiGeoref XArchiGeoref::Applique_transfo(XArchiGeoref G)
{
	//XMat3D R = G.m_Rotation * m_Rotation ;
	XMat3D R = m_Rotation * G.m_Rotation  ;
	XPt3D Tori = m_Rotation * G.m_Translation;
	XPt3D P =  Tori + m_Translation;

	return XArchiGeoref(P,R);
}
//-----------------------------------------------------------------------------
XArchiGeoref XArchiGeoref::Applique_inverse_transfo(XArchiGeoref G)
{
	XMat3D R = m_Rotation * G.m_Rotation;
	XPt3D P = m_Rotation.Trn()*(G.m_Translation-m_Translation);

	return XArchiGeoref(P,R.Trn());
}
//-----------------------------------------------------------------------------
std::string XArchiGeoref::InfoTexte()
{
	std::ostringstream oss;
	oss << "Rattachement du sommet : " << m_Translation << std::endl;
	oss << "Matrice rotation : "  << std::endl;
	oss << m_Rotation << std::endl;

	return oss.str();
}
//-----------------------------------------------------------------------------
bool XArchiGeoref::XmlWrite(std::ostream* out)
{
	std::streamsize prec = out->precision(2);// Sauvegarde des parametres du flux
	std::ios::fmtflags flags = out->setf(std::ios::fixed);

	*out << "<georef>" << std::endl;
		Translation().XmlWrite(out);
		out->precision(12);
		XQuaternion quaternion = XQuaternion(Rotation());
		quaternion.XmlWrite(out);
	*out << "</georef>" << std::endl;

	out->precision(prec);		// Restauration des parametres du flux
	out->unsetf(std::ios::fixed);
	out->setf(flags);
	return out->good();
}

