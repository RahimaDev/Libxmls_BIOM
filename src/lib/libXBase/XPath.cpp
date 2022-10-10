#include "libXBase/XPath.h"
#include "libXBase/XBase.h"
#include <vector>

#ifdef WIN32
char	XPath::PathSep='\\';
#else
char	XPath::PathSep='/';
#endif
//-----------------------------------------------------------------------------
// Renvoie le nom du repertoire parent
//-----------------------------------------------------------------------------
std::string XPath::Path(const char* filename)
{
	std::string P = filename;
	std::string::size_type pos = P.rfind(PathSep);
	if(pos ==std::string::npos)
		pos = P.rfind('\\');
	if(pos ==std::string::npos)
		pos = P.rfind('/');
	return P.substr(0, pos);

}
//-----------------------------------------------------------------------------
// Renvoie le nom du repertoire parent (vide si pas de repertoire)
//-----------------------------------------------------------------------------
std::string XPath::Folder(const char* filename)
{
	std::string P = filename;
	std::string::size_type pos = P.rfind(PathSep);
	if(pos ==std::string::npos)
		pos = P.rfind('\\');
	if(pos ==std::string::npos)
		pos = P.rfind('/');
	if(pos ==std::string::npos)
		return std::string();
	return P.substr(0, pos);
}
//-----------------------------------------------------------------------------
// Renvoie le nom du dernier sous-repertoire
//-----------------------------------------------------------------------------
std::string XPath::LastSubFolder(const char* file_or_folder_name)
{//*****non teste !!!
	std::string P = file_or_folder_name;
	std::string::size_type pos = P.rfind('.'); 
	if(pos != std::string::npos)//on enleve le nom de fichier
		P= P.substr(0, P.rfind(PathSep));
	
	if (P[P.length()-1] == PathSep)
		P= P.substr(0, P.rfind(PathSep));
	return P.substr(P.rfind(PathSep)+1);
}
//-----------------------------------------------------------------------------
// Renvoie le nom du fichier ou du repertoire (sans le chemin)
//-----------------------------------------------------------------------------
std::string XPath::Name(const char* filename)
{
	std::string P = filename;
	std::string::size_type pos = P.rfind(PathSep);
	if(pos ==std::string::npos)
		pos = P.rfind('\\');
	if(pos ==std::string::npos)
		pos = P.rfind('/');
	P = P.substr(pos+1);
	return P;
}
//-----------------------------------------------------------------------------
// Renvoie le nom du fichier (ou du repertoire)  sans extension
//-----------------------------------------------------------------------------
std::string XPath::NameNoExt(const char* filename)
{
	std::string P = filename;
	std::string::size_type pos = P.rfind(PathSep);
	if(pos ==std::string::npos)
		pos = P.rfind('\\');
	if(pos ==std::string::npos)
		pos = P.rfind('/');
	P = P.substr(pos+1);
	P = P.substr(0,P.rfind('.')); 
	return P.substr(0,P.rfind('.'));//enleve les doubles extensions (.ori.xml par exemple)
}
//-----------------------------------------------------------------------------
std::string XPath::RemoveExtension(const char* filename)
{
	std::string P = filename;
	P = P.substr(0,P.rfind('.')); 
	return P.substr(0,P.rfind('.'));//enleve les doubles extensions (.ori.xml par exemple)
}
//-----------------------------------------------------------------------------
// Renvoie la premiere extension d'un fichier
//-----------------------------------------------------------------------------
std::string XPath::Extension(const char* filename)
{
	std::string P = filename;
	std::string ext = "";
	std::string::size_type pos = P.rfind('.');
	if(pos == std::string::npos)
		return ext;
	ext = P.substr(pos+1);
	return ext;
}
//-----------------------------------------------------------------------------
std::string XPath::Extension2(const char* filename, std::string& ext1)
{
	ext1 = "";
	std::string ext2 = "";

	std::string P = filename;
	std::string::size_type pos = P.rfind('.');
	if(pos == std::string::npos)
		return ext2;
	ext1 = P.substr(pos+1);
	P = P.substr(0,pos); 
	pos = (int)P.rfind('.');
	if(pos == std::string::npos)
		return ext2;
	return  P.substr(pos+1);;
}
//-----------------------------------------------------------------------------
// Renvoie l'extension complete d'un fichier (ori.xml)
std::string XPath::FullExtension(const char* filename)
{
	std::string P = filename;
	std::string ext = "";
	std::string::size_type pos = P.find('.');
	if(pos ==std::string::npos)
		return ext;
	ext = P.substr(pos+1);
	return ext;
}
//-----------------------------------------------------------------------------
// change l'extension ou ajoute  une extension a un fichier
//-----------------------------------------------------------------------------
std::string XPath::ChangeExtension(const char* filename, const char* newExt)
{
	std::string F = filename;
	std::string Fnew = F.substr(0,F.rfind('.'));
	Fnew += "." ;
	Fnew += newExt;
	return Fnew;
}
//-----------------------------------------------------------------------------
// insere une chaine dans le nom du fichier avant l'extension 
// toto.txt avec "_new" devient toto_new.txt
//-----------------------------------------------------------------------------
std::string XPath::InsertBeforeExt(const char* filename, const char* insertion)
{
	std::string F = filename;
	std::string ext = Extension(filename);	
	std::string Fnew = F.substr(0,F.rfind('.'));
	Fnew += insertion ;
	Fnew += "." ;
	Fnew += ext;
	return Fnew;
}
//-----------------------------------------------------------------------------
// Renvoie un chemin relatif
//-----------------------------------------------------------------------------
std::string XPath::Relative(const char* root, const char* path)
{
	std::string R = root;
	SubPathSep(R);//root ne doit pas avoir le \ a la fin

	std::string P = path;
	std::string Rel;

	if (P == R)
		return (std::string)"." + PathSep;

	// Cas simple
	// root = C:\toto
	// path = C:\toto\titi
	// => .\titi
	std::string::size_type pos = P.find(R);
	if (pos == 0) {
		if (P[R.length()] == PathSep) {
			Rel = "." + P.substr(R.length());
			return Rel;
		}
	}

	// Cas difficile
	// root = C:\toto\titi
	// path = C:\toto\tata\tutu
	// => ..\tata\tutu
	std::string root1 = R.substr(0, R.find(PathSep));
	std::string root2 = P.substr(0, P.find(PathSep));
	if (root1.compare(root2) != 0)	// Pas le meme disque
		return P;

	std::vector<std::string> T, U;
	std::string tmp = R;
	while(true) {
		std::string s = tmp.substr(0, tmp.find(PathSep));
		if (s.size() < 1)
			break;
		T.push_back(s);
		if (tmp.find(PathSep) == std::string::npos)
			break;
		tmp = tmp.substr(tmp.find(PathSep)+1);
		if (tmp.size() < 1)
			break;
	}
	tmp = P;
	while(true) {
		std::string s = tmp.substr(0, tmp.find(PathSep));
		if (s.size() < 1)
			break;
		U.push_back(s);
		if (tmp.find(PathSep) == std::string::npos)
			break;
		tmp = tmp.substr(tmp.find(PathSep)+1);
		if (tmp.size() < 1)
			break;
	}

    uint32 i, j;
	for (i = 0; i < T.size(); i++) {
		if (i < U.size())
			if (T[i].compare(U[i]) == 0)
				continue;
		for (j = i; j < T.size(); j++)
			Rel = Rel + ".." + PathSep;
		for (j = i; j < U.size() - 1; j++)
			Rel = Rel + U[j] + PathSep;
		Rel += U[U.size() - 1];
		break;
	}

	return Rel;
}
//-----------------------------------------------------------------------------
// Renvoie un chemin relatif dans le cas ou path est un repertoire
//a utiliser pour etre sur d'avoir le \ a la fin -> ex : ".\orthos_imagettes\"
//-----------------------------------------------------------------------------
std::string XPath::FolderRelative(const char* root, const char* path)
{
	std::string R = root;
	SubPathSep(R);//root ne doit pas avoir le \ a la fin
	std::string P = path;
	std::string Rel;

	int pos = (int)P.find(PathSep);
	if (pos == 0){//chemin reseau
		if (P[P.length()-1] == PathSep)
			return P;
		return P + PathSep;
	}

	if (P == R)
		return (std::string)"." + PathSep;

	// Cas simple
	// root = C:\toto
	// path = C:\toto\titi
        // => .\titi
	//ça c'est ok lorsque le repertoire est un sous repertoire, on a bien le \ a la fin
	pos = (int)P.find(R);
	if (pos == 0) {
		if (P[R.length()] == PathSep) {//pour verifier le cas e:\toto  e:\tototiti\tata
			Rel = "." + P.substr(R.length());
			return Rel;
		}
	}

	// Cas difficile
	// root = C:\toto\titi
	// path = C:\toto\tata\tutu
	// => ..\tata\tutu

	//ça c'etait pas bon dans la fonction "relative", on avait pas le \ a la fin
	std::string root1 = R.substr(0, R.find(PathSep));
	std::string root2 = P.substr(0, P.find(PathSep));
	if (root1.compare(root2) != 0)	// Pas le meme disque
		return P;

	std::vector<std::string> T, U;
	std::string tmp = R;
	while(true) {
		std::string s = tmp.substr(0, tmp.find(PathSep));
		if (s.size() < 1)
			break;
		T.push_back(s);
		if (tmp.find(PathSep) == std::string::npos)
			break;
		tmp = tmp.substr(tmp.find(PathSep)+1);
		if (tmp.size() < 1)
			break;
	}
	tmp = P;
	while(true) {
		std::string s = tmp.substr(0, tmp.find(PathSep));
		if (s.size() < 1)
			break;
		U.push_back(s);
		if (tmp.find(PathSep) == std::string::npos)
			break;
		tmp = tmp.substr(tmp.find(PathSep)+1);
		if (tmp.size() < 1)
			break;
	}

    uint32 i, j;
	for (i = 0; i < T.size(); i++) {
		if (i < U.size())
			if (T[i].compare(U[i]) == 0)
				continue;
		for (j = i; j < T.size(); j++)
			Rel = Rel + ".." + PathSep;
		for (j = i; j < U.size() - 1; j++)
			Rel = Rel + U[j] + PathSep;
		Rel += U[U.size() - 1];
		break;
	}
	//si parfois on a le \ et il ne faut pas le rajouter !!!!!
	if (Rel[Rel.length()-1] == PathSep)
		return Rel;

	return Rel + PathSep;
}

//-----------------------------------------------------------------------------
// Renvoie un chemin absolu
//-----------------------------------------------------------------------------
std::string XPath::Absolute(const char* root, const char* path)
{
	std::string R = root;
	SubPathSep(R);//root ne doit pas avoir le \ a la fin
	
	std::string P = path;
	std::string Abs;

	std::string RacineRelative = "." + PathSep;
	if (P.compare(RacineRelative) == 0)
		return R;

	if (P.length() < 3)
		return P;

	//ajout jpp pour compatibilite avec ancienne version qui rajoutait un \ en trop
	//avant le nom du fichier mnt ex ..\..\Mnt\\Fr200.xml
	for (unsigned int i = 1; i < P.size() - 1; i++)
	{
		if ((P[i] == PathSep) && (P[i + 1] == PathSep))
			P.erase(i, 1);
	}

	// Cas ./toto
	if ((P[0] == '.')&&(P[1] == PathSep)) {
		Abs = R + P.substr(1);
		return Abs;
	}

	// Cas ../toto
	if ((P[0] == '.')&&(P[1] == '.')&&(P[2] == PathSep)) {
		std::string subpath = P.substr(P.find(PathSep) + 1);
		std::string subroot = R.substr(0, R.rfind(PathSep));
		if ((subpath[0] == '.')&&(subpath[1] == '.')&&(subpath[2] == PathSep))
			return Absolute(subroot.c_str(), subpath.c_str());
		Abs = subroot + PathSep + subpath;
		return Abs;
	}

	return P;
}
//-----------------------------------------------------------------------------
// Renvoie un chemin absolu pour un repertoire
//a utiliser quand on est pas sûr d'avoir le \ a la fin -> ex : ".\orthos_imagettes\"
//-----------------------------------------------------------------------------
std::string XPath::FolderAbsolute(const char* root, const char* path)
{
	std::string R = root;
	std::string P = path;
	std::string Abs;
	if(P.empty())
		return R;
	char last = P[P.length()-1] ;

	std::string RacineRelative = "." + PathSep;
	if (P.compare(RacineRelative) == 0)
	{
		last = R[R.length()-1] ;
		if(last != PathSep)
			R += PathSep;
		return R;
	}

	if (P.length() < 2)
		return P;

	if(last != PathSep)
		P += PathSep;

	// Cas ./toto
	if ((P[0] == '.')&&(P[1] == PathSep)) {
		Abs = root + P.substr(1);
		return Abs;
	}

	// Cas ../toto
	if ((P[0] == '.')&&(P[1] == '.')&&(P[2] == PathSep)) {
		std::string subpath = P.substr(P.find(PathSep) + 1);
		std::string subroot = R.substr(0, R.rfind(PathSep));
		if ((subpath[0] == '.')&&(subpath[1] == '.')&&(subpath[2] == PathSep))
			return Absolute(subroot.c_str(), subpath.c_str());
		Abs = subroot + PathSep + subpath;
		return Abs;
	}

	return P;
}

//-----------------------------------------------------------------------------
// Traduit un nom de fichier a la norme Windows
// Les caracteres : \ / : * ? " < > | sont interdits
//-----------------------------------------------------------------------------
std::string XPath::ConvertWindows(const char* filename)
{
	std::string P = filename;
	for (uint32 i = 0; i < P.length(); i++)
		if ((P[i] == '\\') || (P[i] == '/') || (P[i] == ':') || (P[i] == '*') ||
				(P[i] == '?') || (P[i] == '"') || (P[i] == '<') || (P[i] == '>') ||
				(P[i] == '|'))
				P[i] = '_';

	return P;
}
//-----------------------------------------------------------------------------
// Transforme c:\\titi\toto\tata.txt  en c:\\titi\\toto\\tata.txt
//-----------------------------------------------------------------------------
std::string XPath::DoubleBackSlash(const char* filename)
{
	std::string Pini = filename;
	std::string Pnew;
	for (uint32 i = 0; i < Pini.length(); i++)
	{
		Pnew += Pini[i];
		if((Pini[i] == '\\')&& (i < Pini.size()-1)){
			Pnew += '\\';
			if(Pini[i+1] == '\\')
				i = i+1;
		}
	}

	return Pnew;
}


//-----------------------------------------------------------------------------
// Transforme c:\titi\toto\tata.txt  c:/titi/toto/tata.txt
//-----------------------------------------------------------------------------
std::string XPath::Convert(const std::string in,const char sep_avt,const char sep_ap)
{
  std::string out="";
  for(unsigned int i=0;i<in.size();i=i+1)
  {
    char en_cours=in.at(i);
    if(en_cours==sep_avt)
      en_cours=sep_ap;
    out=out+en_cours;
  }
  return out;
 }
//-----------------------------------------------------------------------------
std::string XPath::ConvertPathSep(const std::string path)
{
	char badSep = '/';
	if(PathSep == '/')
		badSep = '\\';

  std::string out="";
  if(path.empty())
	  return out;

  for(unsigned int i=0; i<path.size();i++)
  {
    char car = path.at(i);
	if(car == badSep)
		car = PathSep;
	out += car;
  }
   //on en profite pour enlever les separateurs en trop
  for(unsigned int i =1; i< out.size()-1; i++)
	if((out[i] == PathSep)&&(out[i+1] == PathSep))
		out.erase(i,1);

 return out;
 }
//-----------------------------------------------------------------------------
void XPath::AddPathSep(std::string & path)
{
	if((path.empty())||((*path.rbegin())!= XPath::PathSep))
		path.push_back(XPath::PathSep);
}
//-----------------------------------------------------------------------------
void XPath::SubPathSep(std::string & path)
{
	char last = path[path.length()-1] ;
	if(last == PathSep)
		path = path.substr(0,path.length()-1);
}
