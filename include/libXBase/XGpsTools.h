#ifndef __XGPSTOOLS_H
#define __XGPSTOOLS_H

#include "libXBase/XBase.h"
#include "time.h"

#define NB_GPS_TO_UTC 18
static uint16  YearGpsToUtc[NB_GPS_TO_UTC]={1980,1981,1982,1983,1985,1988,1990,1991,1992,1993,1994,1996,1997,1999,2006,2009,2012,2015};
static uint16 MonthGpsToUtc[NB_GPS_TO_UTC]={1,7,7,7,7,1,1,1,7,7,7,1,7,1,1,1,7,7};//janvier ou juillet


//-----------------------------------------------------------------------------
class XGpsTools
{
protected:
	uint16 m_correctionUtcToGps;

	uint16 CorrectionUtcToGps(uint16 year, uint16 month)
	{
		uint16 correction = 0;
		for(uint16 i=0; i<NB_GPS_TO_UTC; i++)
		{
			if(year <YearGpsToUtc[i])
				return correction;

			if((year == YearGpsToUtc[i])&&(month < MonthGpsToUtc[i]))
					return correction;

			correction = i;
		}
		return correction;
	}

    uint16 m_week_day;

    uint16 CorrectionDayToWeek(uint16 year, uint16 month, uint16 day)
    {
       struct tm * timeinfo;
       //recupere une structure tm valid en prenant le temps courant
       time_t rawtime;
       time ( &rawtime );
       timeinfo = localtime ( &rawtime );
       //change les info de timeinfo avec year month day
       timeinfo->tm_year = year - 1900; //tm commence les annes en 1900
       timeinfo->tm_mon = month - 1; // janvier est le mois 0 en tm
       timeinfo->tm_mday = day;
       timeinfo->tm_hour = 12 ;

       /*appel mktime:  pour mettre a jour timeinfo->tm_wday*/
       mktime (timeinfo);
       return timeinfo->tm_wday;
    }


public:
	XGpsTools(uint16 year, uint16 month)
	{
		m_correctionUtcToGps = CorrectionUtcToGps(year, month);
        m_week_day=0;
	}

    XGpsTools(uint16 year, uint16 month, uint16 day)
    {
        m_correctionUtcToGps = CorrectionUtcToGps(year, month);
        m_week_day =CorrectionDayToWeek(year, month, day);
    }
    XGpsTools(uint32 date_YYMMDD)
    {
        XDate date(date_YYMMDD);
        m_correctionUtcToGps = CorrectionUtcToGps(date.year, date.month);
        m_week_day =CorrectionDayToWeek(date.year, date.month, date.day);
    }

	uint16 CorrectionUtcToGps() {return m_correctionUtcToGps;}
    double CorrectionDayToWeek() {return m_week_day*24*3600.0;}
    uint16 WeekDay() {return m_week_day;}

    //conversion d'un temps en seconde GPS depuis debut semaine UTC en temps GPS jour heure minute seconde
	void HMSDecUTC(double GpsTimeS, uint16& H, uint16&M , double& sec, uint16& JourSemaine)
	{
		double secUTC = GpsTimeS - m_correctionUtcToGps;
		SToHMSDec(secUTC, H, M, sec, JourSemaine);
	}

    //conversion temps en seconde depuis le debut de semaine en jour heure min sec
    //independant du fait que l'on soit en temps GPS ou UTC
	void SToHMSDec(double SDec, uint16& H, uint16&M , double& sec, uint16& JourSemaine)
	{
		double heure = SDec/3600.;
		JourSemaine = (uint16)(floor(heure /24.));
		heure = heure - (JourSemaine * 24.);
		JourSemaine = JourSemaine +1; //pour la journee en cours
		double valheure, valmin;

		double reste = modf(heure,&valheure);
		reste = modf(reste*60.,&valmin);

		H = (uint16) valheure;
		M = (uint16) valmin;
		sec = reste*60.;
	}

    //conversion de seconde GPS depuis debut semaine GPS ->  UTC HHMMSS,SSS + jour
	void HMSDecUTC(double GpsTimeS, double& HMSDecUTC, uint16& JourSemaine)
	{
		double secUTC = GpsTimeS - m_correctionUtcToGps;
		uint16 H;
		uint16 M;
		double sec;
		SToHMSDec(secUTC,H,M,sec,JourSemaine);
		HMSDecUTC = H*10000. + M*100. + sec;
	}
};

#endif //__XGPSTOOLS_H
