#pragma once

#define ECART_MINIMAL_INTERPOLATION   0.000000001

class XInterpolLin
{
protected:
	double m_a;
	double m_t0;
	double m_X0;
	double m_dT;
	double m_dT0;

public:
	XInterpolLin(void);
	~XInterpolLin(void);

	//utilisation avec une initialisation uniquement sur les valeurs encadrantes
	//optimisee pour interpolation de plusieurs attributs sur la meme position d'interpolation T
	bool Intitialize_T(double t,double t0, double t1);
	double Interpol_T(double X0, double X1);

	//utilisation et initialisation complete
	//optimisee pour interpolations multiples de X dans le meme interval avec variations de T
	bool Intitialize_XT(double X0, double X1, double t0, double t1);
	double Interpol_XT(double t);
};
