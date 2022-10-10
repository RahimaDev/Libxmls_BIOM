
#include "libXMls/XAttribute.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

using namespace std;

// split: receives a char delimiter; returns a vector of strings
// By default ignores repeated delimiters, unless argument rep == 1.
vector<string> split(string s, char delim, int rep=0)
{
    vector<string> flds;
    if (!flds.empty()) flds.clear();  // empty vector if necessary
    string work = s.data();
    string buf = "";
    unsigned int i = 0;
    while (i < work.length()) {
        if (work[i] != delim)
            buf += work[i];
        else if (rep == 1) {
            flds.push_back(buf);
            buf = "";
        } else if (buf.length() > 0) {
            flds.push_back(buf);
            buf = "";
        }
        i++;
    }
    if (!buf.empty())
        flds.push_back(buf);
    return flds;
}

XAbstractAttribBlock::XAbstractAttribBlock(XAbstractAttrib * p_attrib, XSecond sec):
    m_size(0), mp_attrib(p_attrib), m_sec(sec) {}

bool XAbstractAttribBlock::SetSizeFromFileIfExists()
{
    if(mp_attrib && FileExists()) {m_size = NObjectsFromFile(); return true;}
    return false;
}

string XAbstractAttribBlock::Filename()
{
    ostringstream oss;
    oss << m_sec << '.' << Ext();
    return oss.str();
}

boost::filesystem::path XAbstractAttribBlock::Path()
{
    boost::filesystem::path full_path = (mp_attrib == NULL?".":mp_attrib->AttributePath());
    full_path /= Filename();
    return full_path;
}

void XAbstractAttribBlock::CreateFolderIfNeeded()
{
    boost::filesystem::path attrib_folder = mp_attrib->AttributePath();
    if(!boost::filesystem::exists(attrib_folder.parent_path()))
    {
        cout << "Creating " << attrib_folder.parent_path() << endl;
        boost::filesystem::create_directory(attrib_folder.parent_path());
    }
    if(!boost::filesystem::exists(attrib_folder))
    {
        cout << "Creating " << attrib_folder << endl;
        boost::filesystem::create_directory(attrib_folder);
    }
}

void XAbstractAttribBlock::SetSize(unsigned int size)
{
    if(m_size && m_size!=size) cout << "ERROR: called XAbstractAttribBlock::SetSize(" << size << ") on an attribute of size " << m_size << endl;
    else m_size = size;
}

void XLinearAttribBlock::Load()
{
    if(is_loaded) return; // already loaded
    boost::filesystem::path block_path = Path();
    ifstream ifs(block_path.c_str(), ios_base::in);
    if(!ifs.good())
    {
        cout << "WARN: XLinearAttribBlock::Load() cannot open " << block_path << endl;
        return;
    }
    string word1="", plus="", word2="";
    int n_objects;
    ifs >> n_objects >> word1 >> m_b >> plus >> m_a >> word2;
    //cout << "time attrib " << m_b << "+idx." << m_a << endl;
    if(word1!="objects") {cout << "ERR: First word is not objects in " << block_path << endl; return;}
    if(plus!="+") {cout << "ERR: No + (or bad position) in " << block_path << endl; return;}
    if(word2!="pulse_index") {cout << "ERR: last word is not pulse_index in " << block_path << endl; return;}
    is_loaded = true;
}

void XLinearAttribBlock::Save()
{
    if(!is_loaded) return;
    if(boost::filesystem::exists(Path()))
        std::cout << "WARN: save overwrites " << Path() << std::endl;
    CreateFolderIfNeeded();
    ofstream ofs(Path().c_str(), ios_base::out);
    if(!ofs.good())
    {
        cout << "ERR: cannot write to " << Path() << endl;
        return;
    }
    ofs.precision(15);
    ofs << Size() << " objects " << m_b << " + " << m_a << " pulse_index" << endl;
}

void XLinearAttribBlock::Free(){is_loaded = false;}

unsigned int XLinearAttribBlock::NLoadedObjects(){if(is_loaded) return Size(); return 0;}

unsigned int XLinearAttribBlock::NObjectsFromFile()
{
    boost::filesystem::path block_path = Path();
    if(!FileExists())
    {
        //cout << "WARN: NObjectsFromFile() " << block_path << " does not exist" << endl;
        return 0;
    } //else cout << "NObjectsFromFile() " << block_path << " exist" << endl;
    ifstream ifs(block_path.c_str(), ios_base::in);
    unsigned int n_objects;
    ifs >> n_objects;
    ifs.close();
    return n_objects;
}

//////////////////////////////////////////
string XAbstractAttrib::AttributeFolder()
{
    return m_object+"-"+Typename()+"-"+m_attribname;
}
string XAbstractAttrib::Filename(XSecond sec)
{
    ostringstream oss;
    oss << sec << '.' << Ext();
    return oss.str();
}
string XAbstractAttrib::Filename(XBlockIndex block_idx)
{
    ostringstream oss;
    oss << Second(block_idx) << '.' << Ext();
    return oss.str();
}
boost::filesystem::path XAbstractAttrib::AttributePath()
{
    boost::filesystem::path full_path = m_ept_path;
    full_path /= AttributeFolder();
    return full_path;
}
boost::filesystem::path XAbstractAttrib::BlockPath(XSecond sec)
{
    boost::filesystem::path full_path = AttributePath();
    full_path /= Filename(sec);
    //cout << "BlockPath " << full_path << endl;
    return full_path;
}

XAbstractAttrib * CreateAttribFromType(string type_name)
{
    XAbstractAttrib * p_attrib=NULL;
    if(type_name == "float32") p_attrib = new XFloatAttrib();
    else if(type_name == "float64") p_attrib = new XDoubleAttrib();
    else if(type_name == "int32") p_attrib = new XIntAttrib();
    else if(type_name == "uint32") p_attrib = new XUIntAttrib();
    else if(type_name == "int8") p_attrib = new XCharAttrib();
    else if(type_name == "uint8") p_attrib = new XUCharAttrib();
    else if(type_name == "bool") p_attrib = new XBoolAttrib();
    else if(type_name == "linear") p_attrib = new XLinearAttrib();
    return p_attrib;
}

XAbstractAttrib * CreateAttribFromFolder(string folder)
{
    vector<string> radical = split(folder, '-');
    if(radical.size() < 3) {cout << "ERROR: less than 3 radicals in " << folder << endl; return NULL;}
    if(radical.size() > 3) {cout << "WARNING: more than 3 radicals in " << folder << endl;}
    // object, T, attribname
    XAbstractAttrib * p_attrib = CreateAttribFromType(radical[1].c_str());
    p_attrib->m_object = radical[0];
    p_attrib->m_attribname = radical[2];
    return p_attrib;
}
