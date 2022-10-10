//Francois Becispahic/IGN/2000-2003/Projet Camera Numerique
//Jean-Pierre Papelard/IGN/2006-2015/SR/MATIS

#ifndef _XPATH_H
#define _XPATH_H

#include <string>
	
class XPath  
{
public:
	static char	PathSep;

	XPath() {;}

	static std::string Convert(const std::string in,const char sep_avt = '\\',const char sep_ap = '/');
	std::string ConvertPathSep(const std::string path);
	std::string Path(const char* filename);
	std::string Folder(const char* filename);
	std::string LastSubFolder(const char* file_or_folder_name);
	std::string Name(const char* filename);
	std::string NameNoExt(const char* filename);
	std::string Extension(const char* filename);
	std::string Extension2(const char* filename, std::string& ext1);
	std::string RemoveExtension(const char* filename);
	std::string FullExtension(const char* filename);
	std::string ChangeExtension(const char* filename, const char* newExt);
	std::string InsertBeforeExt(const char* filename, const char* insertion);

	std::string Relative(const char* root, const char* path);
	std::string FolderRelative(const char* root, const char* path);//idem relative mais pour un repertoire
	std::string Absolute(const char* root, const char* path);
	std::string FolderAbsolute(const char* root, const char* path);

	static std::string ConvertWindows(const char* filename);
	std::string DoubleBackSlash(const char* filename);
	void AddPathSep(std::string & path);
	void SubPathSep(std::string & path);

};

#endif //_XPATH_H
