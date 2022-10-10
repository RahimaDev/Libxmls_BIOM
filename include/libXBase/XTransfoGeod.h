//Francois Becispahic/IGN/2000-2003/Projet Camera Numerique
//JP Papelard/IGN/2000-2015/Projet camera numerique - Pôle maintenance - SR - MATIS
#ifndef _XTRANSFOGEOD_H
#define _XTRANSFOGEOD_H

#include "libXBase/XBase.h"
#include "libXBase/XAlgoGeod.h"
#include "libXBase/XPt3D.h"
#include "libXBase/XMat3D.h"
#include "libXBase/XGrilleGeoide.h"

class XTransfoGeod
{
protected:

    XMat3D m_ENH2GeoCart;

    XGrilleGeoide* m_grid;
	virtual XGrilleGeoide* CreateGrille() = 0;

	double TransfoGrid(double h, double lambdaRad,double phiRad);
	double TransfoGridInverse(double z, double lambdaRad,double phiRad);

public :
	XTransfoGeod() ;
	~XTransfoGeod();

	virtual std::string AltiGridName(){ if(m_grid==NULL) return std::string("None"); return m_grid->Name();}
	bool InitGrille(std::string path, XError* error);
	void UnloadGrille();
	XGrilleGeoide* GrilleGeoide(){return m_grid;}

	virtual int		NbTransfo() { return 1;}
	virtual const char*	Name(int n = -1) { return "Identite";}
	virtual void	InitParameters(int n) {;}
	virtual bool	InitParameters(std::string name) {return true;}
	virtual void	Transfo(double lambda,double phi,double h,UnitAng ua,double *x,double *y,double *z)
													{*x = lambda; *y = phi; *z = h;}
	virtual void	TransfoInverse(double x,double y,double z,UnitAng ua,double *lambda,double *phi,double *h)
													{*lambda = x; *phi = y; *h = z;}

	void Transfo(XPt3D LambdaPhiH, UnitAng ua, XPt3D* XYZ)
		{ Transfo(LambdaPhiH.X, LambdaPhiH.Y, LambdaPhiH.Z, ua, &(XYZ->X), &(XYZ->Y), &(XYZ->Z));}
	void TransfoInverse(XPt3D XYZ, UnitAng ua, XPt3D* LambdaPhiH)
		{ TransfoInverse(XYZ.X, XYZ.Y, XYZ.Z, ua, &(LambdaPhiH->X), &(LambdaPhiH->Y), &(LambdaPhiH->Z));}

	//Transfo combinee d'une matrice rotation
	double ConvergenceMeridien(XPt2D LambdaPhi, UnitAng ua);
	XMat3D TransfoMatrice(XMat3D in,double lambda,double phi);
    XMat3D  TransfoInverseMatrice(XMat3D in,double lambda,double phi);
	virtual void TransfoInverse(double x,double y,double z,XMat3D in,UnitAng ua,double *lambda,double *phi,double *h,XMat3D* out);
	virtual void Transfo(double lambda,double phi,double h,XMat3D in,UnitAng ua,double *x,double *y,double *z,XMat3D* out);

	virtual double ComputeHeading(XPt2D p1,XPt2D p2,UnitAng ua);
	virtual bool Compute_ModConv(double lambda, double phi, UnitAng ua, double *conv, double *modlin, double *K, double *alter){ return true;}

};

#endif //_XTRANSFOGEOD_H
