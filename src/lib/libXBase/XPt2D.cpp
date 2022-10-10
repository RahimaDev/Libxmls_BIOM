// XPt2D.cpp Version initiale : F.Becirspahic - Projet Camera Numerique //2000

#include "libXBase/XPt2D.h"
#include "libXBase/XPt3D.h"

//-----------------------------------------------------------------------------
// Ecriture dans un fichier XML
//-----------------------------------------------------------------------------
bool XPt2D::XmlWrite(std::ostream* out)
{
	*out << "<pt2d> " << std::endl;
	*out << "<x> " << X << " </x>" << std::endl;
	*out << "<y> " << Y << " </y>" << std::endl;
	*out << "</pt2d>" << std::endl;
	return out->good();
}
//-----------------------------------------------------------------------------
// Conversion en XPt3D
//-----------------------------------------------------------------------------
XPt2D::operator XPt3D() const
{
	return XPt3D(X, Y, 0);
}

//-----------------------------------------------------------------------------
// Operateurs de calcul
//-----------------------------------------------------------------------------
XPt2D operator+(XPt2D A, XPt2D B)
{
	XPt2D C = A;
	return C += B;
}

XPt2D operator-(XPt2D A, XPt2D B)
{
	XPt2D C = A;
	return C -= B;
}

XPt2D operator*(XPt2D A, double k)
{
	XPt2D B = A;
	return B *= k;
}

XPt2D operator*(double k, XPt2D A)
{
	XPt2D B = A;
	return B *= k;
}

XPt2D operator/(XPt2D A, double k)
{
	XPt2D B = A;
	return B /= k;
}

XPt2D operator/(double k, XPt2D A)
{
	XPt2D B = A;
	return B /= k;
}

//-----------------------------------------------------------------------------
// Operateurs logiques
//-----------------------------------------------------------------------------
bool operator==(XPt2D A, XPt2D B)
{
	return A.X==B.X && A.Y==B.Y; 
}

bool operator!=(XPt2D A, XPt2D B)
{
	return !(A==B);
}

//-----------------------------------------------------------------------------
// Fonctions de calcul de distances
//-----------------------------------------------------------------------------
double dist(XPt2D A, XPt2D B)			// Distance
{
	return sqrt(dist2(A, B));
}

double dist2(XPt2D A, XPt2D B)		// Distance au carre
{
	return (B.X - A.X) * (B.X - A.X) + (B.Y - A.Y) * (B.Y - A.Y);
}

//-----------------------------------------------------------------------------
// Fonctions vectorielles
//-----------------------------------------------------------------------------
double prodScal(XPt2D A, XPt2D B)	// Produit Scalaire
{
	return A.X*B.X + A.Y*B.Y;
}

double prodCross(XPt2D A, XPt2D B, XPt2D C)	// Produit en croix
{
	return (B.X - A.X)*(C.Y - A.Y) - (C.X - A.X)*(B.Y - A.Y);
}

