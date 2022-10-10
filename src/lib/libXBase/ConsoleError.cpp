#include "libXBase/ConsoleError.h"

#ifdef WIN32
	#include <windows.h>
# endif
#include <iostream>

namespace ConsoleError
{

	void Error(const char* origine, const char* mes,  const char* data,bool ShowOrigine)
	{
		#ifdef WIN32	
			SetConsoleTextAttribute( GetStdHandle(STD_OUTPUT_HANDLE),FOREGROUND_RED);
			if(ShowOrigine)
				std::cout << "Erreur " << origine <<std::endl;
			char strmes[1024];
			CharToOem (mes,strmes);
			std::cout <<  strmes << data <<std::endl;
			SetConsoleTextAttribute( GetStdHandle(STD_OUTPUT_HANDLE),FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE);
		#else
			std::cout << RED << "Erreur " << origine << mes << data << DEFAULT_COLOR <<std::endl;
		#endif
	}

	void Alert(const char* origine, const char* mes,  const char* data,bool ShowOrigine)	
	{
		#ifdef WIN32	
			SetConsoleTextAttribute( GetStdHandle(STD_OUTPUT_HANDLE),FOREGROUND_RED | FOREGROUND_GREEN );
			if(ShowOrigine)
				std::cout << "Alerte " << origine <<std::endl;
			char strmes[1024];
			CharToOem (mes,strmes);
			std::cout << strmes << data <<std::endl;
			SetConsoleTextAttribute( GetStdHandle(STD_OUTPUT_HANDLE),FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE);
		#else
			std::cout << YELLOW <<  "Alerte " << origine << mes << data << DEFAULT_COLOR <<std::endl;
		#endif
	}

	void Info(const char* origine, const char* mes,  const char* data,bool ShowOrigine)
	{
		#ifdef WIN32	
			SetConsoleTextAttribute( GetStdHandle(STD_OUTPUT_HANDLE), FOREGROUND_GREEN );
			if(ShowOrigine)
				std::cout << "Info " << origine <<std::endl;
			char strmes[1024];
			CharToOem (mes,strmes);
			std::cout << strmes << data <<std::endl;
			SetConsoleTextAttribute( GetStdHandle(STD_OUTPUT_HANDLE),FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE);
		#else
			std::cout << GREEN << "Info " << origine << mes << data << DEFAULT_COLOR <<std::endl;
		#endif
	}
	void Commentaire(const char* origine, const char* mes, bool ShowOrigine)
	{
			if(ShowOrigine)
				std::cout << origine <<'\n';
			std::cout << mes <<std::endl;
	}
}
