#pragma once

#include <boost/date_time/posix_time/posix_time.hpp>
#include <iostream>

/// Define an absolute time (YY-MM-DD HH:MM:SS) in UTC
/// Used to define the zero for sensor times expressed in seconds
/// from week or day start and in UTC or GPS
/// Riegl: DayStart, UTC
/// Applanix: WeekStart, GPS
/// Landins: WeekStart, UTC
///

class XGpsTools;

class XAbsoluteTime
{
public:
    enum Zero {WeekStart, DayStart};
    enum TimeRef {UTC, GPS};
    boost::posix_time::ptime m_ptime;

    /// Constructor from boost posix_time
    XAbsoluteTime(boost::posix_time::ptime ptime):m_ptime(ptime){}

    /// Constructor from date_time parts
    XAbsoluteTime(int year, int month, int day,
                  int hour, int min, float sec);

    /// Constructor from nothing
    XAbsoluteTime();

    /// Constructor for pivots handling week_start or day start and UTC/GPS.
    XAbsoluteTime(int year, int month, int day, Zero zero, TimeRef time_ref);

    /// Constructor for pivots handling week_start or day start and UTC/GPS, date in YYMMDD format
    XAbsoluteTime(unsigned int date_YYMMDD, Zero zero, TimeRef time_ref);

    /// Constructor by reading an xml. If filename is not an xml,
    /// it is interpreted as a data_filename so PivotName(filename) is read instead
    XAbsoluteTime(std::string filename);

    ~XAbsoluteTime(void){}

    /// functions from original class from JPP preserved for backwards compatibility
    void SetTime(unsigned int date_YYMMDD, unsigned int corJour = 0, unsigned int corSec = 0);
    /// Use for Applanix
    void SetTimeTrajGps(unsigned int date_YYMMDD);
    /// Use for Landins
    void SetTimeTrajUtc(unsigned int date_YYMMDD);
    /// Use for Riegl (.rxp)
    void SetTimeRxpUtc(unsigned int date_YYMMDD);

    /// Load from an xml file
    bool XmlLoad(std::string xml_filename);
    /// Save to a stream
	bool XmlWrite(std::ostream* out);
    /// Save to an xml file
    bool XmlWrite(std::string data_filename);

    /// pivots should be in same folder as data, with same filename but .pvt.xml extension
    std::string PivotName(std::string data_filename);

    /// Add seconds to the time (use with care)
    XAbsoluteTime PlusSecs(int n_sec){return XAbsoluteTime(m_ptime + boost::posix_time::seconds(n_sec));}

    /// Add minutes to the time (use with care)
    XAbsoluteTime PlusMinutes(int n_min){return XAbsoluteTime(m_ptime + boost::posix_time::minutes(n_min));}

    /// Add hours to the time (use with care)
    XAbsoluteTime PlusHours(int n_hour){return XAbsoluteTime(m_ptime + boost::posix_time::hours(n_hour));}

    /// seconds from 01/01/1970
    unsigned long Pivot();

    /// Produce a human readable string expressing the date_time
    std::string ToString();

    /// what you have to add to a time expressed in this time frame to be expressed in the time frame of t
    double Diff(const XAbsoluteTime & t) const;

private:
    /// corrections for day/week start and GPS/UTC
    void Correct(XGpsTools * p_xgt, Zero zero=DayStart, TimeRef time_ref=UTC);
};
