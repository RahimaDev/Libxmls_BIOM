//Jean-Pierre Papelard/IGN/2000-2015/Projet camera numerique/Pôle maintenance/MATIS
#ifndef _XWGS2LAMBERT_H
#define _XWGS2LAMBERT_H

#include  "XTransfoGeod.h"

class XWgs2Lambert : public XTransfoGeod 
{
public :
	enum Projection { None = 0, Lambert1 = 1, Lambert2 = 2, Lambert3 = 3, Lambert4 = 4, Lambert2E = 5};

protected :
  Projection m_Proj;

  double m_dElg1a;//"GRS 80"
  double m_dElg1e2;

  double m_dElg2a;//"Clarke 1880 FR"
  double m_dElg2e2;

  double m_dTx;//"WGS84 vers NTF"
  double m_dTy;
  double m_dTz;
  double m_dd;
  double m_dRx;
  double m_dRy;
  double m_dRz;

  double m_dx0;
  double m_dy0;
  double m_dlambda0; //meridien de Paris
  double m_dphi0;    //valeur initiale en grades
  double m_dk0;

  double m_dn;//Parametres calcules de la projection
  double m_dc;
  double m_dlambdac;
  double m_dxs;
  double m_dys;

	virtual XGrilleGeoide* CreateGrille() { m_grid   = new XRaf98(); return m_grid;}

public :

	XWgs2Lambert(Projection proj = None);

	virtual int		NbTransfo() { return 6;}
	virtual const char*	Name(int n = -1);
	virtual void	InitParameters(int n);
	virtual bool	InitParameters(std::string name) ;
	virtual void  Transfo(double lambda,double phi,double h,UnitAng ua,double *x,double *y,double *z);
	virtual void	TransfoInverse(double x,double y,double z,UnitAng ua,double *lambda,double *phi,double *h);
	virtual bool Compute_ModConv(double lambda, double phi, UnitAng ua, double *conv, double *modlin, double *K, double *alter);

};


#endif //_XWGS2LAMBERT_H
