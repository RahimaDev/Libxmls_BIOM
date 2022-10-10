//-----------------------------------------------------------------------------
//								XErrorXml.cpp
//								=============
//
// Auteur : F.Becirspahic - Projet Camera Numerique
//
// 22/08/00
//-----------------------------------------------------------------------------

#include "libXBase/XErrorXml.h"
#include <sstream>

/*
NOTE:
L'ecriture de donnees dans le fichier de log se fait en effectuant un append dans le fichier
juste avant le tag de fin...

Cette operation qui peut sembler extremement contre performante ne l'est pas
car le fichier est alors entierement en cache du disque

Cependant, si le fichier de log prend une taille considerable ça peut etre lent

Par contre, le disque du coup se retrouve tres sollicite a la fois en lecture et
ecriture...

Il ne faut donc pas trop le faire bosser en meme temps
*/



//-----------------------------------------------------------------------------
// Message d'erreur
//-----------------------------------------------------------------------------
void XErrorXml::Error(const char* origine, const char* mes, Type t)
{
	m_nError++;
	std::string s(origine);
	s += TypeString(t) + "\n" + mes;

	if(m_bFluxOnOutput)
	{
		*m_log << "<Error>" << s << "</Error>" << std::endl;
		return;
	}
	m_log->seekp((uint32)m_log->tellp() - (uint32)m_tagEnd.length());
	*m_log << "<Error>" << s << "</Error>" << std::endl << m_tagEnd;
}
//-----------------------------------------------------------------------------
void XErrorXml::Error(const char* origine, const char* mes, const char* data)
{
	m_nError++;
	std::ostringstream oss;
	oss << "<Error>" << std::endl;
	oss << "<origine>" << origine << "</origine>" << std::endl;
	oss << "<donnee>" << data << "</donnee>" << std::endl;
	oss << "<message>" << mes << "</message>" << std::endl;
	oss << "</Error>" << std::endl;

	if(m_bFluxOnOutput)
	{
		*m_log << oss.str();
		return;
	}
	m_log->seekp((uint32)m_log->tellp() - (uint32)m_tagEnd.length());
	*m_log << oss.str() << m_tagEnd;
}
//-----------------------------------------------------------------------------
void XErrorXml::Alert(const char* origine, const char* mes, const char* data)
{
	m_nAlert++;
	std::ostringstream oss;
	oss << "<Alert>" << std::endl;
	oss << "<origine>" << origine << "</origine>" << std::endl;
	oss << "<donnee>" << data << "</donnee>" << std::endl;
	oss << "<message>" << mes << "</message>" << std::endl;
	oss << "</Alert>" << std::endl;

	if(m_bFluxOnOutput)
	{
		*m_log << oss.str();
		return;
	}
	m_log->seekp((uint32)m_log->tellp() - (uint32)m_tagEnd.length());
	*m_log << oss.str() << m_tagEnd;
}
//-----------------------------------------------------------------------------
void XErrorXml::Info(const char* origine, const char* mes, const char* data)
{
	std::ostringstream oss;
	oss << "<Info>" << std::endl;
	oss << "<origine>" << origine << "</origine>" << std::endl;
	oss << "<donnee>" << data << "</donnee>" << std::endl;
	oss << "<message>" << mes << "</message>" << std::endl;
	oss << "</Info>" << std::endl;

	if(m_bFluxOnOutput)
	{
		*m_log << oss.str();
		return;
	}
	m_log->seekp((uint32)m_log->tellp() - (uint32)m_tagEnd.length());
	*m_log << oss.str() << m_tagEnd;
}

//-----------------------------------------------------------------------------
// Message d'information
//-----------------------------------------------------------------------------
void XErrorXml::Info(const char* origine, const char* mes, Type t)	
{
	std::string s(origine);
	s += TypeString(t) + "\n" + mes;

	if(m_bFluxOnOutput)
	{
		*m_log << "<Info>" << s << "</Info>" << std::endl;
		return;
	}

	m_log->seekp((uint32)m_log->tellp() - (uint32)m_tagEnd.length());
	*m_log << "<Info>" << s << "</Info>" << std::endl << m_tagEnd;
}
//-----------------------------------------------------------------------------
// Message d'alerte
//-----------------------------------------------------------------------------
void XErrorXml::Alert(const char* origine, const char* mes, Type t)
{
	m_nAlert++;
	std::string s(origine);
	s += TypeString(t) + "\n" + mes;

	if(m_bFluxOnOutput)
	{
		*m_log << "<Alert>" << s << "</Alert>" << std::endl;
		return;
	}

	m_log->seekp((uint32)m_log->tellp() - (uint32)m_tagEnd.length());
	*m_log << "<Alert>" << s << "</Alert>" << std::endl << m_tagEnd;
}

//-----------------------------------------------------------------------------
//Ecriture d'un tag complet
//-----------------------------------------------------------------------------
void XErrorXml::Tag(const char* mes, const char* tag)	
{
	if(!m_bFluxOnOutput)
		m_log->seekp((uint32)m_log->tellp() - (uint32)m_tagEnd.length());
	
	*m_log << "<";
	*m_log << tag;
	*m_log << ">";
	*m_log << mes;
	*m_log << "</";
	*m_log << tag;
	*m_log << ">";
	*m_log << std::endl ;

	if(!m_bFluxOnOutput)
		*m_log << m_tagEnd ;
}

//-----------------------------------------------------------------------------
//Commentaire XML
//-----------------------------------------------------------------------------
void XErrorXml::Commentaire(const char* origine, const char* mes, Type t)	
{
	std::string s(origine);
	s += TypeString(t) + "\n" + mes;

	if(m_bFluxOnOutput)
	{
		*m_log << "<!--" << s << "-->" << std::endl;
		return;
	}

	m_log->seekp((uint32)m_log->tellp() - (uint32)m_tagEnd.length());
	*m_log << "<!--" << s << "-->" << std::endl << m_tagEnd;
}

//-----------------------------------------------------------------------------
// Sortie d'un flux
//-----------------------------------------------------------------------------
void XErrorXml::WriteToOutput(std::ostream* out)
{
	m_bFluxOnOutput = true;
	m_log->seekp((uint32)m_log->tellp() - (uint32)m_tagEnd.length());
	*m_log << out->rdbuf() << m_tagEnd;
	m_bFluxOnOutput = false;
}

void XErrorXml::BeginOutput()
{
	m_bFluxOnOutput = true;
	m_log->seekp((uint32)m_log->tellp() - (uint32)m_tagEnd.length());
}

void XErrorXml::EndOutput()
{
	m_bFluxOnOutput = false;
	*m_log << m_tagEnd;
}

//-----------------------------------------------------------------------------
// Debut et fin de tag XML
//-----------------------------------------------------------------------------
void XErrorXml::StartTag(std::string& s)
{
	if(m_bFluxOnOutput)
		return;
	m_log->seekp((uint32)m_log->tellp() - (uint32)m_tagEnd.length());
	m_tagEnd = "</" + s + ">" + m_tagEnd;
	*m_log << "<" << s << ">" << m_tagEnd;
}

void XErrorXml::EndTag()
{
	if(m_bFluxOnOutput)
		return;
	int pos = (int) m_tagEnd.find('>');
	m_tagEnd = m_tagEnd.substr(pos+1);
}

