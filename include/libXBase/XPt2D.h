#ifndef _XPT2D_H
#define _XPT2D_H
//Francois Becispahic/IGN/2000-2003/Projet Camera Numerique

#include "libXBase/XBase.h"

class XPt3D;

class XPt2D 
{
public:
	double	X;
	double	Y;

	// Constructeurs
	XPt2D(double x=0,double y=0): X(x), Y(y) {}

	virtual bool XmlWrite(std::ostream* out);

	// Operations
	inline XPt2D& operator+=(XPt2D M){X+=M.X;Y+=M.Y; return *this;}
	inline XPt2D& operator-=(XPt2D M){X-=M.X;Y-=M.Y; return *this;}
	inline XPt2D& operator*=(double k){X*=k; Y*=k; return *this;}
	inline XPt2D& operator/=(double k){X/=k; Y/=k; return *this;}

	// Operateur de conversion
	operator XPt3D() const;

	double Norme(){return sqrt((X*X) + (Y*Y));}
	void Normalise()
	{
		double n=Norme();
		if (n != 0){
		  X /= n;
		  Y /= n;
		}
	}
};

// Fonctions d'aide

// Operateurs de calcul
XPt2D operator+(XPt2D, XPt2D);
XPt2D operator-(XPt2D, XPt2D);

XPt2D operator*(XPt2D, double);
XPt2D operator*(double, XPt2D);
XPt2D operator/(XPt2D, double);
XPt2D operator/(double, XPt2D);

// Operateurs logiques
bool operator==(XPt2D, XPt2D);
bool operator!=(XPt2D, XPt2D);

// Fonctions de calcul de distances
double dist(XPt2D, XPt2D);			// Distance
double dist2(XPt2D, XPt2D);			// Distance au carre

// Fonctions vectorielles
double prodScal(XPt2D, XPt2D);					// Produit Scalaire
double prodCross(XPt2D, XPt2D, XPt2D);	// Produit en croix


#endif //_XPT2D_H
