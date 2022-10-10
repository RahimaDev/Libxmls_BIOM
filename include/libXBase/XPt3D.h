//Francois Becispahic/IGN/2000-2003/Projet Camera Numerique
#ifndef _XPT3D_H
#define _XPT3D_H

#include <fstream>
#include "libXBase/XPt2D.h"

class XPt2D;

class XPt3D 
{
public:
    double	X;
    double	Y;
    double	Z;

    // Constructeurs
    XPt3D(double x=0,double y=0,double z=0): X(x), Y(y), Z(z) {}
    XPt3D(const XPt3D & p0, const XPt3D & p1, double alpha)
    {
        double uma=1-alpha;
        X = uma*p0.X + alpha*p1.X;
        Y = uma*p0.Y + alpha*p1.Y;
        Z = uma*p0.Z + alpha*p1.Z;
    }

    bool XmlWrite(std::ostream* out);

    // Operations
    inline XPt3D& operator+=(XPt3D M){X+=M.X;Y+=M.Y;Z+=M.Z; return *this;}
    inline XPt3D& operator-=(XPt3D M){X-=M.X;Y-=M.Y;Z-=M.Z; return *this;}
    inline XPt3D& operator*=(double k){X*=k; Y*=k; Z*=k; return *this;}
    inline XPt3D& operator/=(double k){operator*=(1./k); return *this;}

    // Operateur de conversion
    operator XPt2D() const;

    double Norme(){return sqrt((X*X) + (Y*Y) +(Z*Z));}
    void Normalise()
    {
        double n = Norme();
        if(n != 0) operator*=(1./n);
    }
    void BinaryRead(std::ifstream& in);
    void BinaryWrite(std::ofstream& out) const;

};

// Fonctions d'aide

// Operateurs de calcul
XPt3D operator+(XPt3D, XPt3D);
XPt3D operator-(XPt3D, XPt3D);

XPt3D operator*(XPt3D, double);
XPt3D operator*(double, XPt3D);
XPt3D operator/(XPt3D, double);
XPt3D operator/(double, XPt3D);

// Operateurs logiques
bool operator==(XPt3D, XPt3D);
bool operator!=(XPt3D, XPt3D);

// Fonctions de calcul de distances
double dist(XPt3D, XPt3D);			// Distance
double dist2(XPt3D, XPt3D);			// Distance au carre
double dist_plani(XPt3D, XPt3D);	// Distance planimetrique
double dist_plani2(XPt3D, XPt3D);	// Distance planimetrique au carre
double dist_alti(XPt3D, XPt3D);		// Distance altimetrique
double dist_polar(XPt3D);			// Distance polaire en plani
double dist_polar2(XPt3D);			// Distance polaire carre en plani

// Fonctions vectorielles
double prodScal(XPt3D, XPt3D);		// Produit Scalaire
XPt3D prodVect(XPt3D, XPt3D);			// Produit Vectoriel
double prodMixt(XPt3D, XPt3D, XPt3D);	// Produit Mixte

// Fonctions d'entree / sortie
std::istream& operator>>(std::istream&, XPt3D&);	// entree
std::ostream& operator<<(std::ostream&, XPt3D);	// sortie

double PseudoInter2Droites(const XPt3D& S1,const XPt3D& V1,const XPt3D& S2,const XPt3D& V2,XPt3D& out);

#endif //_XPT3D_H
