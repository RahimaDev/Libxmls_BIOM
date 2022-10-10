//Francois Becispahic/IGN/2000-2003/Projet Camera Numerique
#ifndef _XRECT_H
#define _XRECT_H

#include "libXBase/XBase.h"
#include "libXBase/XFrame.h"
class XPoint 
{
public:
	uint32	X;
	uint32	Y;
public:
	XPoint(uint32 x=0, uint32 y=0): X(x), Y(y) {;}
};


class XRect 
{
public:
	uint32	X;
	uint32	Y;
	uint32	W;
	uint32	H;

	XRect(uint32 x=0, uint32 y=0, uint32 w=0, uint32 h=0) : X(x), Y(y), W(w), H(h) {;}
	XRect(XFrame f) {
		X = (uint32)floor(f.Xmin) ;
		Y = (uint32)floor(f.Ymin) ;
		W = (uint32)ceil(f.Xmax) - X;
		H = (uint32)ceil(f.Ymax) - Y;
	}

	// Operations
	XRect& operator+=(XRect r);
	virtual bool XmlWrite(std::ostream* out);

	uint32 centerX(){return X + W/2;}
	uint32 centerY(){return Y + H/2;}
};

#endif //_XRECT_H
