#include <ctime>
#include <time.h>
#include <fstream>
#include "libXBase/XBase.h"
#include "libXBase/XArchiXMLBaseTools.h"
#include "libXBase/XGpsTools.h"
#include "libXBase/XAbsoluteTime.h"
#include <boost/filesystem.hpp>

using namespace std;
using namespace boost::posix_time;

void ToYYYYMMDD(unsigned int & date_YYMMDD)
{
    if(date_YYMMDD < 1000000) date_YYMMDD += 20000000; // interpret 2 first digits as offset to 2000
}
ptime ToBoostPtime(unsigned int & date_YYMMDD)
{
    ToYYYYMMDD(date_YYMMDD);
    return ptime(boost::gregorian::date(YY(date_YYMMDD), MM(date_YYMMDD), DD(date_YYMMDD)));
}

/// Generic time correction for pivots
void XAbsoluteTime::Correct(XGpsTools * p_xgt, Zero zero, TimeRef time_ref)
{
    if(zero == WeekStart) m_ptime = m_ptime - boost::gregorian::days(p_xgt->WeekDay());
    //cout << "UTC " << m_ptime << endl;
    if(time_ref == GPS)
    {
        //cout << "p_xgt->CorrectionUtcToGps() " << p_xgt->CorrectionUtcToGps() << endl;
        m_ptime = m_ptime - seconds(p_xgt->CorrectionUtcToGps());
        //cout << "GPS " << m_ptime << endl;
    }
}

XAbsoluteTime::XAbsoluteTime(int year, int month, int day,
                             int hour, int min, float sec):
    m_ptime(boost::gregorian::date(year, month, day), time_duration(hour, min, sec)){}

XAbsoluteTime::XAbsoluteTime():m_ptime(boost::gregorian::date(1970, 1, 1), time_duration(0, 0, 0)){}

XAbsoluteTime::XAbsoluteTime(int year, int month, int day, Zero zero, TimeRef time_ref)
{
    m_ptime = ptime(boost::gregorian::date(year, month, day));
    XGpsTools xgt(10000*year + 100*month + day);
    Correct(&xgt, zero, time_ref);
}

XAbsoluteTime::XAbsoluteTime(unsigned int date_YYMMDD, Zero zero, TimeRef time_ref)
{
    m_ptime = ToBoostPtime(date_YYMMDD);
    XGpsTools xgt(date_YYMMDD);
    Correct(&xgt, zero, time_ref);
}

XAbsoluteTime::XAbsoluteTime(string filename)
{
    boost::filesystem::path path(filename);
    if(path.extension().string() == ".xml") XmlLoad(filename);
    else XmlLoad(PivotName(filename));
}

//-----------------------------------------------------------------------------
void XAbsoluteTime::SetTime(unsigned int date_YYMMDD, unsigned int corJour, unsigned int corSec)
{
    m_ptime = ToBoostPtime(date_YYMMDD);
    m_ptime = m_ptime - boost::gregorian::days(corJour) - seconds(corSec);
}

void XAbsoluteTime::SetTimeTrajGps(unsigned int date_YYMMDD)
{
    XGpsTools xgt(date_YYMMDD);
    SetTime(date_YYMMDD, xgt.WeekDay(), xgt.CorrectionUtcToGps());
}

void XAbsoluteTime::SetTimeTrajUtc(unsigned int date_YYMMDD)
{
    XGpsTools xgt(date_YYMMDD);
    SetTime(date_YYMMDD, xgt.WeekDay(), xgt.CorrectionUtcToGps());
}

void XAbsoluteTime::SetTimeRxpUtc(unsigned int date_YYMMDD)
{
    return SetTime(date_YYMMDD, 0, 0);
}

unsigned long XAbsoluteTime::Pivot()
{
    time_duration diff = m_ptime - ptime(boost::gregorian::date(1970,1,1));
    return diff.total_seconds();
}

bool XAbsoluteTime::XmlLoad(string xml_filename)
{
    TiXmlDocument doc( xml_filename.c_str() );
    if ( ! doc.LoadFile() )
    {
        cout << "Cannot open xml file " << xml_filename << endl;
        return false;
    }
    //TiXmlHandle hDoc( &doc );
    TiXmlElement* root = doc.RootElement();
    XArchiXML::AssertRoot(root,"pivot_ref_utc");
    m_ptime = time_from_string(XArchiXML::ReadAssertNodeAsString(root,"posix_time"));
    uint32 sec_1970 = XArchiXML::ReadAssertNodeAsUint32(root,"sec_1970");
    if(sec_1970 != Pivot())
        cout << "WARNING: inconsistent sec_1970: read " << sec_1970 << " expected " << Pivot() << " diff " << (int)sec_1970 - (int)Pivot() << endl;
    //cout << xml_filename << " -> " << ToString() << endl;
    return true;
}

//-----------------------------------------------------------------------------
bool XAbsoluteTime::XmlWrite(std::ostream* out)
{
    *out << "<pivot_ref_utc>\n";
    *out << "<posix_time>" << to_simple_string(m_ptime) << "</posix_time>\n";
    *out << "<sec_1970>" << Pivot() << "</sec_1970>\n";
    *out << "</pivot_ref_utc>\n";
    return out->good();
}

bool XAbsoluteTime::XmlWrite(string data_filename)
{
    std::ofstream ofs(PivotName(data_filename).c_str());
    cout << "Writing pivot: " << PivotName(data_filename) << " for file " << data_filename << endl;
    return XmlWrite(&ofs);
}

/// pivots should be in same folder as data, with same filename but with an added .pvt.xml extension
std::string XAbsoluteTime::PivotName(std::string data_filename)
{
    boost::filesystem::path path = data_filename + ".pvt.xml";
    return path.string();
}

std::string XAbsoluteTime::ToString()
{
    return to_simple_string(m_ptime);
}

double XAbsoluteTime::Diff(const XAbsoluteTime & t) const
{
    time_duration diff = m_ptime - t.m_ptime;
    return diff.total_seconds();
}
