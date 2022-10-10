//Jean-Pierre Papelard/IGN/2006-2015/SR/MATIS
#ifndef _X_ARCHI_XML_EXCEPTION_H_
#define _X_ARCHI_XML_EXCEPTION_H_

#include <exception>
#include <string>

namespace XArchiXML
{

	class XmlException  : public std::exception
	{
		public:
			enum TypeExceptionXml { eNull, eBadFormat, eTagObligatoire, eTagVide, eRacine, eBadData, eInaccessible, eAllocation};

		protected:
			TypeExceptionXml m_type;
			std::string m_erreur;

		public:
			XmlException(const char* erreur, TypeExceptionXml t = eNull) {m_erreur = erreur; m_type = t;}
			virtual ~XmlException() throw() {;}

			std::string Erreur(){return TypeString() + m_erreur;}
			std::string TypeString()
			{
				switch(m_type) {
				case eNull:
					return "";
				case eBadFormat:
					return "Erreur de format de fichier ";
				case eTagObligatoire:
					return "Tag obligatoire manquant :";
				case eTagVide:
					return "Champ obligatoire non rempli :";
				case eRacine:
					return "Racine de l'arborescence incompatible ";
				case eBadData:
					return "Erreur de coherence des donnees ";
				case eInaccessible:
					return "Donnee inaccessible  ";
				case eAllocation:
					return "Echec allocation  ";
				default:
					return "Erreur inconnue";
				}
			}

	};
}
#endif