
//-----------------------------------------------------------------------------
//								SystemInfo.cpp
//								==============
//
//-----------------------------------------------------------------------------

#include "libXBase/XSystemInfo.h"
#include "libXBase/XStringTools.h"
#include "libXBase/XPath.h"
#include "libXBase/XArchiTime.h"
#include <fstream>
#include <iostream>
#include <sstream>
//#define BOOST_FILESYSTEM_VERSION 2
#define BOOST_FILESYSTEM_NO_DEPRECATED
#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/operations.hpp>

#define NOT_IMPLEMENTED std::cerr << __FILE__ << ":" << __LINE__ << std::endl << __FUNCTION__ << " not implemented in linux yet" << std::endl; // throw

// using namespace boost::filesystem;
namespace fs = boost::filesystem;

bool XSystemInfo::GetFileListInFolder(std::string folder, std::vector<std::string> &liste, std::string filtre, bool sort ){
    fs::directory_iterator dir_iter(folder), dir_end;
    if(filtre.at(0)!='*'){
        std::cout<<"Pas gere pour le moment";
        return false;
    }
    std::string filtre_bis=filtre.substr(1);
    for(;dir_iter != dir_end; ++dir_iter)
    {
        //std::cout<<dir_iter->leaf()<<std::endl;;
        if(boost::filesystem::is_directory(*dir_iter))
            continue;

        if(filtre=="*.*")
            liste.push_back(dir_iter->path().string());
        else if(dir_iter->path().extension() == filtre_bis)
            liste.push_back(dir_iter->path().string());
    }
    return true;
}

int XSystemInfo::CountFileInFolder(std::string folder, std::string filtre)
{
    if(filtre.at(0)!='*'){
        std::cout<<"Pas gere pour le moment";
        return 0;
    }
    std::string filtre_bis=filtre.substr(1);
    int n=0;
    fs::directory_iterator dir_iter(folder), dir_end;
    for(;dir_iter != dir_end; ++dir_iter)
        if(!boost::filesystem::is_directory(*dir_iter) &&
                (filtre=="*.*" || dir_iter->path().extension() == filtre_bis))
            n++;
    return n;
}

void XSystemInfo::GetFileListInAllSubFolder(std::string folder, std::vector<std::string> &liste, std::string filtre){
  fs::directory_iterator dir_iter(folder), dir_end;
  if(filtre.at(0)!='*'){
    std::cout<<"Pas gere pour le moment";
    return;
  }
  std::string filtre_bis=filtre.substr(1);
  for(;dir_iter != dir_end; ++dir_iter)
  {
    //std::cout<<dir_iter->leaf()<<std::endl;;
    if(boost::filesystem::is_directory(*dir_iter))
      GetFileListInAllSubFolder(dir_iter->path().string(), liste, filtre);
    
    if(filtre=="*.*")
      liste.push_back(dir_iter->path().string());
    else if(dir_iter->path().extension() == filtre_bis)
      liste.push_back(dir_iter->path().string());
  }
  return;
}

bool XSystemInfo::NumberofElementsInFolderIsGreaterThan(std::string folder, int nbElementsMax)
{
    int n=0;
    fs::directory_iterator dir_iter(folder), dir_end;
    for(;dir_iter != dir_end; ++dir_iter) if(n++ > nbElementsMax) return true;
    return false;
}


std::vector<std::string> XSystemInfo::GetSubFolderListInFolder(std::string folder, bool multiLevel)
{
    std::vector<std::string>  liste;
    fs::directory_iterator dir_iter(folder), dir_end;
    for(;dir_iter != dir_end; ++dir_iter)
    {
        if(fs::is_directory(*dir_iter))
        {
            if((dir_iter->path().filename() != ".")&&(dir_iter->path().filename() != ".."))
                liste.push_back(dir_iter->path().string());
        }
    }
    if(!multiLevel) return liste;
    /*
	unsigned int size_avt,size_apres;
	size_avt=liste.size();
	size_apres=liste.size()+1;
	while(size_avt!=size_apres){

		for(int i=0;i<	liste.size();i=i+1)

		size_avt=size_apres;
		size_apres=liste.size()
	}
	//for(int i=0;i<
*/
    for(unsigned int i=0; i < liste.size(); i++)
    {
        std::vector<std::string> tmp = GetSubFolderListInFolder(liste[i], true);
        //			for(int j=0;j<tmp.size();j=j+1){
        //				std::cout<<tmp.at(j)<<std::endl;
        //		}
        //std::cout<<"-----------------------------------------------"<<std::endl;
        if(!tmp.empty())
            liste.insert(liste.end(), tmp.begin(), tmp.end());
    }
    return liste;
}

_int64 XSystemInfo::GetFileSize1(const char* src)
{
    return fs::file_size(src);
}

double XSystemInfo::GetFileSize2(const char* src)
{
    return fs::file_size(src);
}

//-----------------------------------------------------------------------------
std::string XSystemInfo::FileSizeToString(double fileSizeDbl)
{
	std::string units("bytes");

	if (fileSizeDbl > 1024.0)
	{
		fileSizeDbl /= 1024.0;
		units = "Kb";
	}
	if (fileSizeDbl > 1024.0)
	{
		fileSizeDbl /= 1024.0;
		units = "Mb";
	}
	if (fileSizeDbl > 1024.0)
	{
		fileSizeDbl /= 1024.0;
		units = "Gb";
	}

	char fileSizeStr[1024];
	sprintf(fileSizeStr,"%.2f %s", fileSizeDbl, units.c_str());

	return std::string(fileSizeStr);
}

bool XSystemInfo::FindFile(const char* filename)
{
    fs::path p(filename);
    return fs::exists(p) && fs::is_regular_file(p);
}
bool XSystemInfo::FindFolder(std::string folder)
{
    fs::path p(folder);
    return fs::exists(p) && fs::is_directory(p);
}

bool XSystemInfo::Copy_File(const char* src, const char* dst, bool FailIfExist)
{
    fs::path from(src), to(dst);
    if(!fs::exists(from) || !fs::is_regular_file(from)) return false;
    if(FailIfExist && fs::exists(to)) return false;
    fs::copy_file(from, to);
    return true;
}

bool XSystemInfo::CopyDirContent(const char* src, const char* dst, bool FailIfExist)
{
    fs::path from(src), to(dst);
    if(fs::exists(to) && (FailIfExist || !fs::is_directory(to))) return false;
    // copy_directory(from, to); not implemented in boost 1.40
    NOT_IMPLEMENTED;
    return true;
}
bool XSystemInfo::CopyFileInFolder(const char* fullFilename, const char* folderDest)
{
    fs::path from(fullFilename), to(folderDest);
    if(!fs::exists(from) || !fs::is_regular_file(from) || (fs::exists(to) && !fs::is_directory(to))) return false;
    if(!fs::exists(to)) fs::create_directory(to);
    fs::copy_file(from, to);
    return true;
}
bool XSystemInfo::CreateMultiDirectory(const char* FullDirectoryPath )
{
    fs::path p(FullDirectoryPath);
    fs::create_directory(p);
    return true;
}
bool XSystemInfo::CreateDirectoryIfNotExist(const char* FullDirectoryPath )
{
    fs::path p(FullDirectoryPath);
    if(fs::exists(p)) return true;
    return fs::create_directory(p);
}

#include <boost/date_time.hpp>
std::string XSystemInfo::DateSysteme(std::string sep)
{
    std::ostringstream oss; oss << boost::gregorian::day_clock::local_day();
    return std::string(oss.str());
}
std::string XSystemInfo::HeureSysteme(std::string sep){
    std::ostringstream oss; oss << boost::posix_time::second_clock::local_time();
    return std::string(oss.str());
}
void XSystemInfo::System_Time(std::string &strDate, std::string & strHeure){
    std::ostringstream oss_date; oss_date << boost::gregorian::day_clock::local_day();
    strDate = oss_date.str();
    std::ostringstream oss_heure; oss_heure << boost::posix_time::second_clock::local_time();
    strHeure = oss_heure.str();
}

int XSystemInfo::NbProcessors(){
  return sysconf( _SC_NPROCESSORS_ONLN );
}

std::string XSystemInfo::ComputerName()
{
//   DWORD buf_size = 1024;
//   char buffer[1024];
//   if (::GetComputerName(buffer,&buf_size)== FALSE)
//     return std::string();
  char hostname[1024];
  hostname[1023] = '\0';
  ::gethostname(hostname, 1023);
  return std::string(hostname);
}

std::string XSystemInfo::GetExePath()
{
    // NOTE: not very robust, see remark on http://stackoverflow.com/a/5694199/4186941
    return fs::current_path().string();
}

// std::string get_selfpath()
std::string  XSystemInfo::GetCurrentExe()
{
  char buff[1024];
  ssize_t len = ::readlink("/proc/self/exe", buff, sizeof(buff)-1);
  if (len != -1) {
    buff[len] = '\0';
    return std::string(buff);
  } else {
    return "";
  }
}

std::string XSystemInfo::CreateLogFile(std::ofstream* log, const char * extension)
{
  XPath p;
  
  std::string AppPath = GetCurrentExe();
  std::string AppName = p.NameNoExt(AppPath.c_str());
  
  std::string folder = p.Path(AppPath.c_str()) + XPath::PathSep + "Journal_" + AppName;
  if (!CreateDirectoryIfNotExist(folder.c_str()))
    return std::string();
  folder += XPath::PathSep;
  
  //date
  XArchiTime AT;
  uint16 Y,M,D;
  AT.Date(Y,M,D);
  char buf[1024];
  sprintf(buf,"%04d_%02d_%02d",Y+2000,M,D);
  std::string strTime = buf;
  std::string filename;
  int nblogToday = 0;
  do {
    sprintf(buf,"%s%s_%s_%s_%d.%s",folder.c_str(),AppName.c_str(),ComputerName().c_str(),strTime.c_str(),nblogToday,extension);
    nblogToday++;
  } while( fs::exists(fs::path(buf)) );//enddo
  
  log->open(buf, std::ios_base::out);
  if(!log->good())
    return std::string();
  log->setf(std::ios_base::unitbuf);
  
  return std::string(buf);
}

void XSystemInfo::StripDirectoryPath(std::string &input)
{
    fs::path path(input);
    input = path.filename().string();
}


