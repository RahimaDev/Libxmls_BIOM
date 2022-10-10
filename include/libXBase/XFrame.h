//Francois Becispahic/IGN/2000-2003/Projet Camera Numerique
#ifndef _XFRAME_H
#define _XFRAME_H

#include "libXBase/XPt2D.h"
#include "libXBase/XPt3D.h"

class XRect;
class XFrame 
{
public:
	double	Xmin;
	double	Ymin;
	double	Xmax;
	double	Ymax;

	XFrame(double xmin=0, double ymin=0, double xmax=0, double ymax=0) :
					Xmin(xmin), Ymin(ymin), Xmax(xmax), Ymax(ymax) {;}

	XFrame(XRect rect);

	bool IsValid();

	// Operations
	XFrame& operator+=(XFrame r);
	XFrame& operator+=(XPt2D P);
	XFrame& operator+=(double marge);
	XFrame& operator+=(XPt3D P) { return *this += XPt2D(P.X, P.Y);}
	XFrame& operator*=(double k);
				
    bool Intersect(const XFrame& r) const;
    bool IsIn(const XPt2D& P) const;
	bool Include(XFrame& r);
	inline bool IsEmpty() { return ((Xmin == 0)&&(Ymin == 0)&&(Xmax == 0)&&(Ymax == 0));}
	inline double Width() { return Xmax - Xmin;}
	inline double Height(){ return Ymax - Ymin;}

	XPt2D Center() const{ return XPt2D((Xmax + Xmin)*0.5, (Ymax + Ymin)*0.5);}
	XPt2D NW() {  return XPt2D(Xmin,Ymax);}
	XPt2D NE() {  return XPt2D(Xmax,Ymax);}
	XPt2D SW() {  return XPt2D(Xmin,Ymin);}
	XPt2D SE() {  return XPt2D(Xmax,Ymin);}

	bool XmlWrite(std::ostream* out);
};
	// Operateurs logiques
	bool operator==(XFrame, XFrame);
	bool operator!=(XFrame, XFrame);


#endif //_XFRAME_H

