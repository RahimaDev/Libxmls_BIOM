
#define WHAT "GenerateTimePivot: Generate a time pivot file for a given datafile given its date and time convention"

#include "libXBase/XAbsoluteTime.h"

using namespace std;

//-----------------------------------------------------------------------------
int main(int argc, char **argv)
{
    cout << WHAT << endl;
    if(argc < 5)
    {
        cout << "Usage: " << argv[0] << "  datafile1 [datafile2 ...] date Zero TimeRef" << endl;
        cout << "datafileX: name of the data files for which a pivot needs to be generated" << endl;
        cout << "date: date of acquisition in YYMMDD or YYYYMMDD format. WARNING: not sure what happens for acquisitions passing midnight" << endl;
        cout << "Zero: is the zero of that file at week_start (0) or day_start (1)" << endl;
        cout << "TimeRef: is time in UTC (0) or GPS (1)" << endl;
        cout << "Riegl: 1 0" << endl;
        cout << "Applanix: 0 1" << endl;
        cout << "Landins: 0 0" << endl;
        return 0;
    }

    int i_arg=argc-3;
    // required
    unsigned int date_YYMMDD = atoi(argv[i_arg++]);
    unsigned int Zero = atoi(argv[i_arg++]);
    unsigned int TimeRef = atoi(argv[i_arg++]);
    //cout << "date_YYMMDD=" << date_YYMMDD << " Zero=" << Zero << " TimeRef=" << TimeRef << endl;

    for(i_arg=1; i_arg<argc-3; i_arg++)
    {
        XAbsoluteTime abstime(date_YYMMDD,
                              (Zero?XAbsoluteTime::DayStart:XAbsoluteTime::WeekStart),
                              (TimeRef?XAbsoluteTime::GPS:XAbsoluteTime::UTC));
        abstime.XmlWrite(std::string(argv[i_arg]));
    }
    cout << argc-4 << " pivot(s) written" << endl;
    return 0;
}

