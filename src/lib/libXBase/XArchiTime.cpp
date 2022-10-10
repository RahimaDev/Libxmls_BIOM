#include "libXBase/XArchiTime.h"

void XArchiTime::ConvertStringToHMS(char * data,int &hour, int &min, int &second, double &second_dec,unsigned int separator_size)
{
	char *cursor=data;
	char temp[3];

	//recuperation de l'heure
	memcpy(temp,cursor,2);
	temp[2]='\0';
	hour=atoi(temp);

	cursor+=2+separator_size;

	//recuperation des minutes
	memcpy(temp,cursor,2);
	temp[2]='\0';
	min=atoi(temp);

	cursor+=2+separator_size;

	//recuperation des secondes
	memcpy(temp,cursor,2);
	temp[2]='\0';
	second=atoi(temp);

	int offset=4+2*separator_size;

	char *toto=new char [strlen(data)-offset+1];
	strncpy(toto,data+offset,strlen(data)-offset);
	toto[strlen(data)-offset]='\0';
	second_dec=atof(toto);
	delete [] toto;

}
