#include "libXBase/XBBox.h"
#include <limits>

XBBox::XBBox(void)
{
	min.X = std::numeric_limits<double>::max();
	max.X = std::numeric_limits<double>::min();
	min.Y = std::numeric_limits<double>::max();
	max.Y = std::numeric_limits<double>::min();
	min.Z = std::numeric_limits<double>::max();
	max.Z = std::numeric_limits<double>::min();
}

XBBox::~XBBox(void)
{
}
