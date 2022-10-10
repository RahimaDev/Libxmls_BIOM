//Jean-Pierre Papelard/IGN/2006-2015/SR/MATIS

#ifndef _XRAF98
#define _XRAF98

#include "libXImage2/XImage.h"
#include "libXBase/XFInterpol.h"
//-----------------------------------------------------------------------------
class XGrilleGeoide{
protected :
	XImage<double> * m_Image;

protected:
	double GetFastFloatPix(double x, double y, XFInterpol<double>* inter);
public :

	XGrilleGeoide()  {m_Image = NULL;}
	virtual~XGrilleGeoide();

	bool IsLoaded(){return (m_Image != NULL);}
	virtual std::string Name() = 0;

	virtual bool Load(std::ifstream& in) = 0;
	virtual bool Interpol(double longitude, double latitude, double* Hwgs84) = 0;
};

//-----------------------------------------------------------------------------
class XRaf98 : public XGrilleGeoide {
public :
	XRaf98() : XGrilleGeoide()  {;}

	std::string Name(){return  std::string("RAF98");}

	bool Load(std::ifstream& in);
	bool Interpol(double longitude, double latitude, double* Hwgs84);
};

//-----------------------------------------------------------------------------
class XEgm96 : public XGrilleGeoide {
public :

	XEgm96():XGrilleGeoide()  {}

	std::string Name(){return  std::string("Egm96");}

	bool Load(std::ifstream& in);
	bool Interpol(double longitude, double latitude, double* Hwgs84);
};

#endif//_XRAF98
