//Francois Becispahic/IGN/2000-2003/Projet Camera Numerique
#ifndef _XFINTERPOL_H
#define _XFINTERPOL_H

#include <limits>

#include "libXBase/XBase.h"
#include "libXBase/XPt3D.h"

//-----------------------------------------------------------------------------
// XFInterpol : interpolation au plus proche voisin
//-----------------------------------------------------------------------------
template<class T> class XFInterpol {
protected:
	uint16	m_nWin;		// Taille de la fenetre
	T*			m_Y;			// Tableau des valeurs en Y
	uint32	m_nOffX;	// Offset en X

	void Init(uint16 win);
	virtual inline double ComputeDouble(T* value, double x)
								{ if (x <= 0.5) return value[0]; else return value[m_nOffX];}

public:
	XFInterpol<T>() { Init(1);}
	virtual ~XFInterpol() { delete[] m_Y;}

	inline uint16 Win() const { return m_nWin;}
	inline double Compute(T* value, double x) { return ComputeDouble(value, x);}
	double BiCompute(T* value, double x, double y, uint32 offset = 0);
	void OffsetX(uint32 off = 0) { m_nOffX = off + 1;}
};

//-----------------------------------------------------------------------------
// Constructeur
//-----------------------------------------------------------------------------
template<class T> void XFInterpol<T>::Init(uint16 win)
{
	m_nWin = win;
	m_Y = new T[2 * m_nWin];
	m_nOffX = 1;
}

//-----------------------------------------------------------------------------
// XInterlLin : interpolation lineaire
//-----------------------------------------------------------------------------
template<class T> class XFInterLin : public XFInterpol<T> {
public :
	XFInterLin<T>() { XFInterpol<T>::Init(1);}
protected:
	virtual inline double ComputeDouble(T* value, double x)
								{ return (double)(value[XFInterpol<T>::m_nOffX] - value[0]) * x + value[0];}
};

//-----------------------------------------------------------------------------
// XInterlCub : interpolation cubique
//-----------------------------------------------------------------------------
template<class T> class XFInterCub : public XFInterpol<T> {
public :
	XFInterCub<T>() { XFInterpol<T>::Init(2);}
protected:
	virtual double ComputeDouble(T* value, double x);
};

//-----------------------------------------------------------------------------
// Calcul sur un tableau de reels
//-----------------------------------------------------------------------------
/*template<class T> double XFInterCub<T>::ComputeDouble(T* value, double x)
{
	double a, b, c, d = value[m_nOffX];
	double dz1 = value[0] - d, dz2 = value[2*m_nOffX] - d, dz3 = value[3*m_nOffX] - d;

	a = (-1. * dz1 - 3. * dz2 + 1. * dz3) / 6.;
	b = ( 3. * dz1 + 3. * dz2) / 6.;
	c = (-2. * dz1 + 6. * dz2 - 1. * dz3) / 6.;
	return (a*x*x*x + b*x*x + c*x + d);
}*/

#ifdef max
#undef max
#endif //max

#ifdef min
#undef min
#endif //min

template<class T> double XFInterCub<T>::ComputeDouble(T* value, double x)
{
	float a, b, c, d = value[XFInterpol<T>::m_nOffX];
	float dz1 = value[0] - d, dz2 = value[2*XFInterpol<T>::m_nOffX] - d, dz3 = value[3*XFInterpol<T>::m_nOffX] - d;

	a = (-1. * dz1 - 3. * dz2 +      dz3) / 6.;
	b = (      dz1 +      dz2) / 2.;
	c = (-2. * dz1 + 6. * dz2 -      dz3) / 6.;
	float res = ((a*x + b)*x + c)*x + d;

	if (res > std::numeric_limits<T>::max())
		return std::numeric_limits<T>::max();
	if(res < std::numeric_limits<T>::min())
		return std::numeric_limits<T>::min();
	return res;
}

//-----------------------------------------------------------------------------
// Calcul sur un tableau de type quelconque
//-----------------------------------------------------------------------------
template<class T> double XFInterpol<T>::BiCompute(T* value, double x, double y, uint32 offset)
{
	for (uint16 i = 0; i < 2*m_nWin; i++)
		m_Y[i] = Compute(&value[i*(2*m_nWin + offset)], x);
	return ComputeDouble(m_Y, y);
}




#endif //_XFINTERPOL_H
