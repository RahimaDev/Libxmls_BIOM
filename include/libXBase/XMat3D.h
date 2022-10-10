//Francois Becispahic/IGN/2000-2003/Projet Camera Numerique
#ifndef _XMAT3D_H
#define _XMAT3D_H

#include "libXBase/XPt3D.h"

class XMat3D 
{
public:
	enum TypeRotation { KapaPhiOmega, KapaOmegaPhi, KPO, OPK};
//protected:
public:
    XPt3D A;	// 1er ligne de la matrice
    XPt3D B;	// 2em ligne de la matrice
    XPt3D C;	// 3em ligne de la matrice

public:
	// Constructeurs
	XMat3D(void) : A(), B(), C() { }
    XMat3D(double ax, double ay, double az,
           double bx, double by, double bz,
           double cx, double cy, double cz) : A(ax,ay,az), B(bx,by,bz), C(cx,cy,cz) { }
	XMat3D(XPt3D a, XPt3D b, XPt3D c) : A(a), B(b), C(c) { }
	XMat3D(XPt3D);	// Initialisation d'un axiateur
	XMat3D(double omega, double phi, double kapa, TypeRotation type = KPO);	// Matrice rotation

	// Acces au donnees membres
	XPt3D lig(int n) const;	// Renvoie la ligne n
	XPt3D col(int n) const;	// Renvoie la colonne n

	// Operations
	XMat3D& operator+=(XMat3D);
	XMat3D& operator-=(XMat3D);
	XMat3D& operator*=(XMat3D);
	XMat3D& operator*=(double);
	XMat3D& operator/=(double);

	double Det() const;	// Determinant
	XMat3D  Trn() const;	// Transposee

	// Calcul des angles d'une matrice rotation
	double Omega(TypeRotation = KPO) const;
	double Phi(TypeRotation = KPO) const;
	double Kapa(TypeRotation = KPO) const;

	//Calcul Axiateur et angle a partir de R
	void R2Axe_Angle(XPt3D* p,double* angle)const;
	//Calcul R a partir de Axiateur et angle  
	static XMat3D Axe_Angle2R(XPt3D p,double angle);
	//Calcul R=this + R(da,db,dc)
  XMat3D R_plus_dR(double da,double db,double dc);
	
	virtual bool XmlWrite(std::ostream* out);
	virtual bool WriteTxt(std::ostream* out);
	virtual bool ReadTxt(std::ifstream* in);

	void BinaryRead(std::ifstream& in);
	void BinaryWrite(std::ofstream& out) const;

	static XMat3D MPI(){return XMat3D(XPt3D(-1,0,0),XPt3D(0,-1,0),XPt3D(0,0,0));}
	static XMat3D Identite(){return XMat3D(XPt3D(1,0,0),XPt3D(0,1,0),XPt3D(0,0,1));}
	static XMat3D Null(){return XMat3D(XPt3D(0,0,0),XPt3D(0,0,0),XPt3D(0,0,0));}

	  void Normalise();	
};

XMat3D operator+(XMat3D, XMat3D);
XMat3D operator-(XMat3D, XMat3D);

XMat3D operator*(XMat3D, XMat3D);
XMat3D operator*(XMat3D, double);
XMat3D operator*(double, XMat3D);
XPt3D operator*(XMat3D, XPt3D);
XPt3D operator*(XPt3D, XMat3D);

XMat3D prod_Pa_tPb(XPt3D Pa, XPt3D Pb);

bool operator==(XMat3D, XMat3D);
bool operator!=(XMat3D, XMat3D);

std::istream& operator>>(std::istream&, XMat3D&);
std::ostream& operator<<(std::ostream&, XMat3D);


#endif //_XMAT3D_H
