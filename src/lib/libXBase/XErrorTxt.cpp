#include "libXBase/XErrorTxt.h"

XErrorTxt::XErrorTxt(std::ostream* log)
{
	m_log = log;
}

XErrorTxt::~XErrorTxt(void)
{
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void XErrorTxt::Error(const char* origine, const char* mes, Type t)
{
	m_nError++;
	if(m_log != NULL)
		*m_log << OrigineMessageData((std::string("ERREUR ")+ TypeString(t)).c_str(),origine,mes,"");
}
//-----------------------------------------------------------------------------
void XErrorTxt::Error(const char* origine, const char* mes, const char* data)
{
	m_nError++;
	if(m_log != NULL)
		*m_log << OrigineMessageData("ERREUR",origine,mes,data);
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void XErrorTxt::Alert(const char* origine, const char* mes, Type t)
{
	m_nAlert++;
	if(m_log != NULL)
		*m_log << OrigineMessageData((std::string("ALERTE")+ TypeString(t)).c_str(), origine,mes,"");
}
//-----------------------------------------------------------------------------
void XErrorTxt::Alert(const char* origine, const char* mes, const char* data)
{
	m_nAlert++;
	if(m_log != NULL)
		*m_log << OrigineMessageData("ALERTE",origine,mes,data);
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void XErrorTxt::Info(const char* origine, const char* mes, Type t)	
{
	if(m_log != NULL)
		*m_log << OrigineMessageData((std::string("INFO ")+ TypeString(t)).c_str(),origine,mes,"");
}
//-----------------------------------------------------------------------------
void XErrorTxt::Info(const char* origine, const char* mes, const char* data)
{
	if(m_log != NULL)
		*m_log << OrigineMessageData( "INFO",origine,mes,data);
}
//-----------------------------------------------------------------------------
 void XErrorTxt::Commentaire(const char* origine, const char* mes)
 {
	if(m_log != NULL)
		*m_log << origine << ' ' << mes << '\n';
}
