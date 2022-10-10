#ifndef _XARCHITIME_H
#define _XARCHITIME_H

#include "XBase.h"

#include <time.h>

typedef struct
{
	struct tm time_tm;
	double seconde; //pour un timming de precision
} UTC_time_t;

class XArchiTime
{
protected:
        // OTO:Bon, ca, c'est chiant
        // Sous windows ca passe nickel, pas sous un linux 32bits ...
        // Par ailleurs, la msdn dit (http://msdn.microsoft.com/en-us/library/1f4c8f33(VS.80).aspx):
        /*
        In Visual C++ 2005, time is a wrapper for _time64 and time_t is, by default, equivalent to __time64_t. If you need to force
        the compiler to interpret time_t as the old 32-bit time_t, you can define _USE_32BIT_TIME_T. This is not recommended because
        your application may fail after January 18, 2038; the use of this macro is not allowed on 64-bit platforms.
        */
        // Disons que ca laisse de la marge, donc on met un time_t ...
        #ifdef WIN32
            __time64_t ltime;
        #else
            time_t ltime;
        #endif // WIN32
public:
	XArchiTime(){
        #ifdef WIN32
            _time64( &ltime );
        #else
            time( &ltime );
        #endif // WIN32
	}
	uint32 Date()
	{
        #ifdef WIN32
            tm* today = _localtime64( &ltime );
        #else
            tm* today = localtime( &ltime );
        #endif // WIN32
		uint16 year =  today->tm_year -100;//normal : 0 en 1900
		uint16 month = today->tm_mon +1 ; ////????????? pourquoi +1 (0 pour janvier ???)
		uint16 day = today->tm_mday;
		uint32 ArchiDate = year*10000 +  month*100 + day;
		return ArchiDate;
	}
	void Date(uint16 &year,uint16 &month,uint16 &day)
	{
        #ifdef WIN32
            tm* today = _localtime64( &ltime );
        #else
            tm* today = localtime( &ltime );
        #endif // WIN32
		year =  today->tm_year -100;//normal : 0 en 1900
		month = today->tm_mon +1 ; ////????????? pourquoi +1 (0 pour janvier ???)
		day = today->tm_mday;
	}
	uint16 GmtMinuteInDay()
	{
        #ifdef WIN32
            tm* gmt =_gmtime64( &ltime );
        #else
            tm* gmt =gmtime( &ltime );
        #endif // WIN32
		return gmt->tm_hour * 60 + gmt->tm_min;
	}
	uint32 GmtSecondInDay()
	{
        #ifdef WIN32
            tm* gmt =_gmtime64( &ltime );
        #else
            tm* gmt =gmtime( &ltime );
        #endif // WIN32
		return gmt->tm_hour * 3600 + gmt->tm_min *60 +  gmt->tm_sec ;
	}

	static void ConvertStringToHMS(char * data,int &hour, int &min, int &second,double &second_dec,unsigned int separator_size=0);
    #ifdef WIN32
        __time64_t TimeInSecondes(){return ltime;}
        __time64_t DiffTimeInSecondes(XArchiTime* time2){return time2->TimeInSecondes() - ltime; }
    #else
        time_t TimeInSecondes(){return ltime;}
        time_t DiffTimeInSecondes(XArchiTime* time2){return time2->TimeInSecondes() - ltime; }
    #endif // WIN32

};



#endif //_XARCHITIME_H
