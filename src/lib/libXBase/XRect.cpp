// XRect.cpp Version initiale : F.Becirspahic - Projet Camera Numerique //2000
#include "libXBase/XRect.h"

//-----------------------------------------------------------------------------
// Operations
//-----------------------------------------------------------------------------
XRect& XRect::operator+=(XRect r)
{
	uint32 x = XMin(X, r.X);
	uint32 y = XMin(Y, r.Y);

	W = XMax( (X + W) - x, (r.X + r.W) - x);
	H = XMax( (Y + H) - y, (r.Y + r.H) - y);
	X = x;
	Y = y;
	return *this;
}
//-----------------------------------------------------------------------------
// Ecriture dans un fichier XML
//-----------------------------------------------------------------------------
bool XRect::XmlWrite(std::ostream* out)
{
	*out << "<rect> " << std::endl;
	*out << "<x> " << X << " </x>" << std::endl;
	*out << "<y> " << Y << " </y>" << std::endl;
	*out << "<w> " << W << " </w>" << std::endl;
	*out << "<h> " << H << " </h>" << std::endl;
	*out << "</rect>" << std::endl;
	return out->good();
}

