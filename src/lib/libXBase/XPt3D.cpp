// XPt3D.cpp Version initiale : F.Becirspahic - Projet Camera Numerique //2000

#include <cmath>
#include "libXBase/XPt3D.h"

//-----------------------------------------------------------------------------
// Conversion en XPt2D
//-----------------------------------------------------------------------------
XPt3D::operator XPt2D() const
{
	return XPt2D(X, Y);
}

//-----------------------------------------------------------------------------
// Ecriture dans un fichier XML
//-----------------------------------------------------------------------------
bool XPt3D::XmlWrite(std::ostream* out)
{
	*out << "<pt3d> " << std::endl;
	*out << "<x> " << X << " </x>" << std::endl;
	*out << "<y> " << Y << " </y>" << std::endl;
	*out << "<z> " << Z << " </z>" << std::endl;
	*out << "</pt3d>" << std::endl;
	return out->good();
}


//-----------------------------------------------------------------------------
void XPt3D::BinaryRead(std::ifstream& in)
{
   in.read((char*)&X,sizeof(double));
   in.read((char*)&Y,sizeof(double));
   in.read((char*)&Z,sizeof(double));
}
//-----------------------------------------------------------------------------
void XPt3D::BinaryWrite(std::ofstream& out) const
{
   out.write((char*)&X,sizeof(double));
   out.write((char*)&Y,sizeof(double));
   out.write((char*)&Z,sizeof(double));
}



//-----------------------------------------------------------------------------
// Operateurs de calcul
//-----------------------------------------------------------------------------
XPt3D operator+(XPt3D A, XPt3D B)
{
	XPt3D C = A;
	return C += B;
}

XPt3D operator-(XPt3D A, XPt3D B)
{
	XPt3D C = A;
	return C -= B;
}

XPt3D operator*(XPt3D A, double k)
{
	XPt3D B = A;
	return B *= k;
}

XPt3D operator*(double k, XPt3D A)
{
	XPt3D B = A;
	return B *= k;
}

XPt3D operator/(XPt3D A, double k)
{
	XPt3D B = A;
	return B /= k;
}

XPt3D operator/(double k, XPt3D A)
{
	XPt3D B = A;
	return B /= k;
}

//-----------------------------------------------------------------------------
// Operateurs logiques
//-----------------------------------------------------------------------------
bool operator==(XPt3D A, XPt3D B)
{
	return A.X==B.X && A.Y==B.Y && A.Z ==B.Z; 
}

bool operator!=(XPt3D A, XPt3D B)
{
	return !(A==B);
}

//-----------------------------------------------------------------------------
// Fonctions de calcul de distances
//-----------------------------------------------------------------------------
double dist(XPt3D A, XPt3D B)			// Distance
{
	return sqrt(dist2(A, B));
}

double dist2(XPt3D A, XPt3D B)		// Distance au carre
{
	XPt3D C = A - B;
	return prodScal(C, C);
}

double dist_plani(XPt3D A, XPt3D B)	// Distance planimetrique
{
	return sqrt(dist_plani2(A, B));
}

double dist_plani2(XPt3D A, XPt3D B)	// Distance planimetrique au carre
{
	return (A.X-B.X)*(A.X-B.X) + ((A.Y-B.Y)*(A.Y-B.Y));
}

double dist_alti(XPt3D A, XPt3D B)	// Distance altimetrique
{
	return fabs(A.Z - B.Z);
}

double dist_polar(XPt3D A)			// Distance polaire en plani
{
	return sqrt(dist_polar2(A));
}

double dist_polar2(XPt3D A)			// Distance polaire carre en plani
{
	return A.X*A.X + A.Y*A.Y;
}

//-----------------------------------------------------------------------------
// Fonctions vectorielles
//-----------------------------------------------------------------------------
double prodScal(XPt3D A, XPt3D B)	// Produit Scalaire
{
	return A.X*B.X + A.Y*B.Y + A.Z*B.Z;
}

XPt3D prodVect(XPt3D A, XPt3D B)	// Produit Vectoriel
{
	double x = A.Y*B.Z - B.Y*A.Z;
	double y = B.X*A.Z - A.X*B.Z;
	double z = A.X*B.Y - B.X*A.Y;
	return XPt3D(x, y, z);
}

double prodMixt(XPt3D A, XPt3D B, XPt3D C)		// Produit Mixte
{
	return prodScal(A, prodVect(B, C));
}

//-----------------------------------------------------------------------------
// Fonctions d'entree / sortie
//-----------------------------------------------------------------------------
std::istream& operator>>(std::istream& s, XPt3D& M)
{
	s >> M.X >> M.Y >> M.Z;
	return s;
}

std::ostream& operator<<(std::ostream& s, XPt3D M)
{
	return s << M.X << "\t" << M.Y << "\t" << M.Z;
}

double PseudoInter2Droites(const XPt3D& S1,const XPt3D& V1,const XPt3D& S2,const XPt3D& V2,XPt3D& out){
  XPt3D B,M1,M2,D;
  B=S2-S1;
  double a1=prodScal(V1,V1);
  double a2=prodScal(V2,V2);
  double a12=prodScal(V1,V2);
  double aB1=prodScal(B,V1);
  double aB2=prodScal(B,V2);

  M1=S1+(aB2*a12-aB1*a2)/(a12*a12-a1*a2)*V1;
  M2=S2+(aB2*a1-aB1*a12)/(a12*a12-a1*a2)*V2;
  out=0.5*(M1+M2);
  D=0.5*(M2-M1);
  return D.Norme();
}
