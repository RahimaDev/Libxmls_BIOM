//Francois Becispahic/IGN/2000-2003/Projet Camera Numerique
#ifndef _XERROR_H
#define _XERROR_H

#include "libXBase/XBase.h"
#include <sstream>

//-----------------------------------------------------------------------------
// XError : Gestion du routage des messages d'erreur
//-----------------------------------------------------------------------------
class XError 
{
public:
	uint32 m_nError;	// Compteur d'erreurs
	uint32 m_nAlert;	// Compteur d'alertes
public:
	enum Type { eNull, eAllocation, eRange, eIOOpen, eIOSeek, eIORead, eIOWrite,
							eBadFormat, eUnsupported, eIllegal, eBadData};

	XError() {m_nError = 0; m_nAlert = 0;}
	virtual ~XError() {;}

	virtual void Error(const char* origine, const char* mes, Type t = eNull) {m_nError++;}	// Message d'erreur
	virtual void Alert(const char* origine, const char* mes, Type t = eNull) {m_nAlert++;}	// Message d'alerte
	virtual void Info(const char* origine, const char* mes, Type t = eNull) {;}	// Message d'information

	//--------------------------------------------------------------------------------------------------------
	//nouvelles methodes sans le type enum mais avec le nom de la donn�es concern�e par l'erreur
	//pour un message d'erreur on aurait les 3 infos suivantes :
	//	origine =  nom de classe+methode ayant detect� l'errur
	//	mes = message d'erreur
	//	data = nom de la variable concern�e par l'erreur
	virtual void Error(const char* origine, const char* mes,  const char* data) =0;//{m_nError++;}	// Message d'erreur
	virtual void Alert(const char* origine, const char* mes,  const char* data) {m_nAlert++;}	// Message d'alerte
	virtual void Info(const char* origine, const char* mes, const char* data) {;}	// Message d'information
/*	std::string OrigineMessageData(const char* prefixe,const char* origine, const char* mes, const char* data)
	{ 
		std::stringstream OMD;
		OMD <<  prefixe;
		if(strlen(origine)>0)
			OMD <<  " ORIGINE: " << origine ;
		if(strlen(mes)>0)
			OMD <<  " MESSAGE: " << mes ;
		if(strlen(data)>0)
			OMD <<  " DATA: " << data ;
		OMD << std::endl;
		return OMD.str();
	}
*/
	std::string OrigineMessageData(const char* prefixe,const char* origine, const char* mes, const char* data)
	{ 
		std::stringstream OMD;
		OMD <<  prefixe << '\t' << origine << '\t' << mes << '\t' << data << '\n';
		return OMD.str();
	}

	// Renvoi d'une chaine de caractere explicitant le type d'erreur
	std::string TypeString(Type t)
	{
		switch(t) {
		case eNull:
			return "";
		case eAllocation:
			return "Erreur d'allocation";
		case eRange:
			return "Erreur d'intervale";
		case eIOOpen:
			return "Erreur d'ouverture de fichier";
		case eIOSeek:
			return "Erreur de positionnement dans un fichier";
		case eIORead:
			return "Erreur de lecture dans un fichier";
		case eIOWrite:
			return "Erreur d'�criture dans un fichier";
		case eBadFormat:
			return "Erreur de format de fichier";
		case eUnsupported:
			return "Format de fichier non support�";
		case eIllegal:
			return "Action ill�gale";
		case eBadData:
			return "Erreur de coh�rence des donn�es";
		default:
			return "Erreur inconnue";
		}
	}
	//--------------------------------------------------------------------------------------------------------

	virtual void Commentaire(const char* origine, const char* mes) {;}	
	virtual void Tag(const char* mes, const char* tag) {;}	// ecriture d'un tag xml

	virtual void Reset() { m_nError = 0; m_nAlert = 0;}			// Remise a 0 des compteurs
	inline uint32 NbError() const { return m_nError;}
	inline uint32 NbAlert() const { return m_nAlert;}

	virtual void Output(std::ostream* out) {;}				// Sortie d'un flux
	virtual std::ostream* Output() { return NULL;}
	virtual void BeginOutput() {;}
	virtual void EndOutput() {;}

	friend bool XErrorError(XError* error,const char* origine,  const char* mes, const char* data)
								{ if (error != NULL) error->Error( origine,mes, data); return false; }
	friend bool XErrorInfo(XError* error,const char* origine,  const char* mes, const char* data)
								{ if (error != NULL) error->Info( origine,mes, data); return true; }
	friend bool XErrorAlert(XError* error,const char* origine,  const char* mes,const char* data)
								{ if (error != NULL) error->Alert( origine,mes, data); return true; }

	friend bool XErrorError(XError* error,const char* origine,  const char* mes, XError::Type t = XError::eNull)
								{ if (error != NULL) error->Error( origine,mes, t); return false; }
	friend bool XErrorInfo(XError* error,const char* origine,  const char* mes, XError::Type t = XError::eNull)
								{ if (error != NULL) error->Info( origine,mes, t); return true; }
	friend bool XErrorAlert(XError* error,const char* origine,  const char* mes, XError::Type t = XError::eNull)
								{ if (error != NULL) error->Alert( origine,mes, t); return true; }
	friend bool XErrorCommentaire(XError* error,const char* origine,  const char* mes)
								{ if (error != NULL) error->Commentaire( origine,mes); return true; }
	friend bool XErrorTag(XError* error, const char* tag, const char* mes)
								{ if (error != NULL) error->Tag(mes, tag); return true; }

	friend void XErrorReset(XError* error) { if (error != NULL) error->Reset();}
	friend uint32 XErrorNbError(XError* error) { if (error != NULL) return error->NbError(); return 0;}
	friend uint32 XErrorNbAlert(XError* error) { if (error != NULL) return error->NbAlert(); return 0;}

	friend void XErrorOutput(XError* error, std::ostream* out)
								{ if (error != NULL) error->Output(out);}
	friend std::ostream* XErrorOutput(XError* error)
								{ if (error == NULL) return NULL; else return error->Output();}
	friend void XErrorBeginOutput(XError* error)
								{ if (error != NULL) error->BeginOutput();}
	friend void XErrorEndOutput(XError* error)
								{ if (error != NULL) error->EndOutput();}
};

#endif //_XERROR_H
