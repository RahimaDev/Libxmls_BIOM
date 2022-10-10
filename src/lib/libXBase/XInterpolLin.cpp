#include "libXBase/XInterpolLin.h"
#include <math.h>

XInterpolLin::XInterpolLin(void)
{
}
//-----------------------------------------------------------------------------
XInterpolLin::~XInterpolLin(void)
{
}
//-----------------------------------------------------------------------------
bool XInterpolLin::Intitialize_T(double t, double t0, double t1)
{
	m_dT = t1-t0;
	if(fabs(m_dT)<  ECART_MINIMAL_INTERPOLATION)
		return false;

	m_dT0 = t-t0;
	return true;
}
//-----------------------------------------------------------------------------
double XInterpolLin::Interpol_T(double X0, double X1)
{
	double a =(X1-X0)/m_dT;
	return a*m_dT0 + X0;
}

//-----------------------------------------------------------------------------
bool XInterpolLin::Intitialize_XT(double X0,double X1, double t0, double t1)
{
	m_dT = t1-t0;
	if(fabs(m_dT)<  ECART_MINIMAL_INTERPOLATION)
		return false;

	m_a =(X1-X0)/m_dT;
	m_t0 = t0;
	m_X0 = X0;
	return true;
}
//-----------------------------------------------------------------------------
double XInterpolLin::Interpol_XT(double t)
{
	return m_a*(t-m_t0)+ m_X0;
}
