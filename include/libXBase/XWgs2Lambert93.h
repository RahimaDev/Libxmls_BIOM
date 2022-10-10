//Jean-Pierre Papelard/IGN/2000-2015/Projet camera numerique/Pôle maintenance/MATIS
#ifndef _XWGS2LAMBERT93_H
#define _XWGS2LAMBERT93_H

#include  "XTransfoGeod.h"
#include  "XGrilleGeoide.h"

class XWgs2Lambert93 : public XTransfoGeod
{
protected :

	double m_dElg1a;//"GRS 80"
	double m_dElg1e2;

	double m_dx0;
	double m_dy0;
	double m_dlambda0; //meridien de Paris
	double m_dphi0;    //valeur initiale en grades

	double m_dk0;// inutiles??

	double m_dphi1;
	double m_dphi2;

	double m_dn;//Parametres calcules de la projection
	double m_dc;
	double m_dlambdac;
	double m_dxs;
	double m_dys;

	virtual XGrilleGeoide* CreateGrille() { m_grid = new XRaf98(); return m_grid;}

public :

	XWgs2Lambert93();

	virtual int		NbTransfo() { return 1;}
	virtual const char*	Name(int n = -1) {return "Lambert93";}
	virtual void  Transfo(double lambda,double phi,double h,UnitAng ua,double *x,double *y,double *z);
	virtual void  TransfoInverse(double x,double y,double z,UnitAng ua,double *lambda,double *phi,double *h);
	virtual bool  Compute_ModConv(double lambda, double phi, UnitAng ua, double *conv, double *modlin, double *K, double *alter);
};


#endif //_XWGS2LAMBERT93_H
