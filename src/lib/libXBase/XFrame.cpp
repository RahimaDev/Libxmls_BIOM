#include "libXBase/XFrame.h"
#include "libXBase/XRect.h"

//-----------------------------------------------------------------------------
//Constructeurs
//-----------------------------------------------------------------------------
XFrame::XFrame(XRect rect)
{ 
	Xmin = rect.X ; 
	Xmax = rect.X + rect.W; 
	Ymin = rect.Y; 
	Ymax = rect.Y+rect.H;
}
//-----------------------------------------------------------------------------
// Operateurs logiques
//-----------------------------------------------------------------------------
bool operator==(XFrame A, XFrame B)
{
	return A.Xmin == B.Xmin && A.Xmax==B.Xmax && A.Ymin ==B.Ymin && A.Ymax ==B.Ymax; 
}

bool operator!=(XFrame A, XFrame B)
{
	return !(A==B);
}
//-----------------------------------------------------------------------------
// Operations
//-----------------------------------------------------------------------------
XFrame& XFrame::operator+=(XFrame r)
{
	if (r.IsEmpty())
		return *this;
	if (IsEmpty()) {
		*this = r;
		return *this;
	}

	Xmin = XMin(Xmin, r.Xmin);
	Ymin = XMin(Ymin, r.Ymin);
	Xmax = XMax(Xmax, r.Xmax);
	Ymax = XMax(Ymax, r.Ymax);
	return *this;
}

//-----------------------------------------------------------------------------
XFrame& XFrame::operator+=(XPt2D P)
{
	if (IsEmpty()) {
		Xmin = Xmax = P.X;
		Ymin = Ymax = P.Y;
		return *this;
	}

	Xmin = XMin(Xmin, P.X);
	Ymin = XMin(Ymin, P.Y);
	Xmax = XMax(Xmax, P.X);
	Ymax = XMax(Ymax, P.Y);
	return *this;
}
//-----------------------------------------------------------------------------
XFrame& XFrame::operator+=(double marge)
{
	Xmin -= marge;
	Ymin -= marge;
	Xmax += marge;
	Ymax += marge;
	return *this;
}

//-----------------------------------------------------------------------------
XFrame& XFrame::operator*=(double k)
{
	Xmin -= k*Width();
	Ymin -= k*Height();
	Xmax += k*Width();
	Ymax += k*Height();
	return *this;
}

//-----------------------------------------------------------------------------
// Intersection
//-----------------------------------------------------------------------------
bool XFrame::Intersect(const XFrame& r) const
{
	if ((Xmin > r.Xmax)||(Xmax < r.Xmin))
		return false;
	if ((Ymin > r.Ymax)||(Ymax < r.Ymin))
		return false;
	return true;
}

//-----------------------------------------------------------------------------
// Teste si un point est dans le cadre
//-----------------------------------------------------------------------------
bool XFrame::IsIn(const XPt2D& P) const
{
	if ((P.X < Xmin)||(P.Y < Ymin))
		return false;
	if ((P.X > Xmax)||(P.Y > Ymax))
		return false;
	return true;
}
//-----------------------------------------------------------------------------
bool XFrame::IsValid()
{	if (Xmin >= Xmax)
	   return false;
	if (Ymin >= Ymax)
	   return false;
	return true;
}
//-----------------------------------------------------------------------------
// Teste si un cadre est dans le cadre
//-----------------------------------------------------------------------------
bool XFrame::Include(XFrame& r)
{
	if (!IsIn(XPt2D(r.Xmin,r.Ymax)))
		return false;
	if (!IsIn(XPt2D(r.Xmax,r.Ymin)))
		return false;
	return true;
}
//-----------------------------------------------------------------------------
// Ecriture dans un fichier XML
//-----------------------------------------------------------------------------
bool XFrame::XmlWrite(std::ostream* out)
{
	*out << "<frame> " << std::endl;
	*out << "<xmin> " << Xmin << " </xmin>" << std::endl;
	*out << "<ymin> " << Ymin << " </ymin>" << std::endl;
	*out << "<xmax> " << Xmax << " </xmax>" << std::endl;
	*out << "<ymax> " << Ymax << " </ymax>" << std::endl;
	*out << "</frame>" << std::endl;
	return out->good();
}
