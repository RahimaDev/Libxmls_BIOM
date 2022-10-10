//Francois Becispahic/IGN/2000-2003/Projet Camera Numerique
#ifndef _XBASE_H
#define _XBASE_H

#include <cstdlib>
#include <cstring>
#include <cstdio>

// #ifndef _USE_MATH_DEFINES
// 	#define _USE_MATH_DEFINES
// #endif

#include <cmath>
#include <cassert>
#include <iostream>
#include <string>
#include <vector>
// #include "math.h"

// Definitions des types standards
#ifndef uint32
	typedef	unsigned int uint32;
#endif

#ifndef uint16
	typedef	unsigned short uint16;
#endif

#ifndef byte
	typedef unsigned char byte	;
#endif

#ifndef uint8
	typedef unsigned char uint8	;
#endif

// Definitions de XMin et XMax (en general min et max causent des problemes de
//								definitions multiples)
template<class T> const T& XMin(const T& x, const T& y)
								{ if (x < y) return x; else return y;}
template<class T> const T& XMax(const T& x, const T& y)
								{ if (x > y) return x; else return y;}

// Fonction XRint : arrondi a l'entier le plus proche
inline int XRint( double x ) { return (x - floor(x)<0.5)? (int)floor(x):(int)ceil(x); }
inline long Xround (double x){	double X = floor(x + 0.5);return (long) X;}	
template<class T> inline const T& XAbs(const T& x) { return (x >=0)? x : -x;}

// Macro pour les assertions
#ifndef _DEBUG
	#define XAssert(expr)	(void(0))
#else
	#define XAssert(expr)	assert(expr)
#endif
template<class T> T MetValAuPasInf(T val, T pas)
	{
		double sign = (val>0) ? 1. : -1.;
		double nb;
		double reste = modf(fabs(val)/pas,&nb);
		if(reste ==0)
			return val;

		if(sign>0)
			return nb*pas;
		else
			return sign*(nb+1)*pas;
	}
template<class T> T MetValAuPasSup(T val, T pas)
	{
		double sign = (val>0) ? 1. : -1.;
		double nb;
		double reste = modf(fabs(val)/pas,&nb);
		if(reste ==0)
			return val;

		if(sign>0)
			return (nb+1)*pas;
		else
			return sign*nb*pas;
	}

	enum TypeXInterpol { Interpol_PPVoisin, Interpol_Lineaire, Interpol_BiCubique, Interpol_SinXY};
	//enum class incompris par VS2010
	//enum class TypeAppui { PXYZ, PXY, PZ, PLiaison, PSift};//enum compatible avec une conversion directe vers stereoprep.xml (!= APP)
	enum TypeAppui { PXYZ, PXY, PZ, PLiaison, PSift};//enum compatible avec une conversion directe vers stereoprep.xml (!= APP)
	//TODO : voir consequence ajout PSift

	enum UnitAng    { Radian = 0,SecDec = 1, DegDec = 2, DegMinSecDec = 3, DegMinDec = 4};
	//-------------------------------------------------------------------------
	// Conversion degre.minute seconde decimale vers radian
	//-------------------------------------------------------------------------
	static double sexadms_to_rad(double val)
	{
		double sign = (val>0) ? 1. : -1.;

		double deg;
		double reste = modf(fabs(val),&deg);
		int    min   = (int)floor(reste*100);
		long double sec   = reste*10000 - min*100;

	return sign * M_PI*(deg + min/60. + sec/3600.)/180.;
	}

	//-------------------------------------------------------------------------
	// Conversion degre minute. decimale vers radian
	//-------------------------------------------------------------------------
	static double sexadm_to_rad(double val)
	{
		double sign = (val>0) ? 1. : -1.;

		int    deg = (int)floor(sign * fabs(val)/100);
		double min = val- deg*100.;

	return M_PI*(deg + min/60.)/180.;
	}

	//-------------------------------------------------------------------------
	// Conversion de radian vers degre.minute seconde decimale 
	//-------------------------------------------------------------------------
	static double rad_to_sexadms(double val)
	{
		double sign = (val>0) ? 1. : -1.;

		double deg = 180.*val/M_PI;
		double valdeg, valmin;

		double reste = modf(fabs(deg),&valdeg);
		reste = modf(reste*60,&valmin);

		if(fabs(reste - 1.0) < 0.00001){
			valmin +=1;
			reste = 0.;
		}

		
	return sign * (valdeg + valmin/100. + reste*60/10000.);
	}

	//-------------------------------------------------------------------------
	// Conversion de radian vers degre minute. decimale
	//-------------------------------------------------------------------------
	static double rad_to_sexadm(double val)
	{
		double sign = (val>0) ? 1. : -1.;

		double deg = 180.*val/M_PI;
		double valdeg;

		double reste = modf(fabs(deg),&valdeg);
		
	return sign * (valdeg + reste*60/100.);
	}
		
	//-------------------------------------------------------------------------
	// Conversion en radians de la valeur reçue  
	//-------------------------------------------------------------------------
	static double ConvertUnit2Rad(double val,UnitAng ua)
	{
		switch (ua){	
		case Radian : //Radians  
			return val;
		case SecDec : //Seconde decimale   
			return M_PI_2 *val/5400.;
		case DegDec ://Degres decimaux
			return M_PI*(val/180.);
		case DegMinSecDec ://Degres.MinuteSecondeDecimale
			return sexadms_to_rad(val);
		case DegMinDec ://DegresMinutes.Decimale
			return sexadm_to_rad(val);
	default: //dans les autres cas on renvoie les coordonnees originelles 
			return val;
		}
	}
	//-------------------------------------------------------------------------
	// Conversion de radians vers les unites demandees
	//-------------------------------------------------------------------------
	static double Rad2ConvertUnit(double val,UnitAng ua)
	{
		switch (ua){	
		case Radian : //Radians  
			return val;
		case SecDec : //Seconde decimale   
			return 5400.*val/M_PI_2;
		case DegDec :	//Degres decimaux
			return 180.*val/M_PI;
		case DegMinSecDec :	//Degres.MinuteSecondeDecimale
			return rad_to_sexadms(val);
		case DegMinDec :		//DegresMinutes.Decimale
			return rad_to_sexadm(val);
	default: //dans les autres cas on renvoie les radians 
			return val;
		}
	}	

	//-----------------------------------------------------------------------------
    template<class T> T ConvertUnitAng(T val, UnitAng start, UnitAng end)
	{ 
		return Rad2ConvertUnit(ConvertUnit2Rad(val, start), end);
	}

    static uint16 YY(uint32 date_YYMMDD)
    {
        double x = (double)date_YYMMDD * 0.0001;
        return (uint16)floor(x);
    }
    //-----------------------------------------------------------------------------
    static uint16 MM(uint32 date_YYMMDD)
    {
        double x = (double)date_YYMMDD * 0.01;
        x = floor(x) * 0.01;
        return (uint16)XRint((x - floor(x)) * 100.);
    }
    //-----------------------------------------------------------------------------
    static uint16 DD(uint32 date_YYMMDD)
    {
        double x = (double)date_YYMMDD * 0.01;
        return (uint16)XRint((x - floor(x)) * 100.);
    }
    struct XDate
    {
        uint16 year, month, day;
        XDate(uint32 date_YYMMDD):year(YY(date_YYMMDD)), month(MM(date_YYMMDD)), day(DD(date_YYMMDD)){}
    };

#endif //_XBASE_H
