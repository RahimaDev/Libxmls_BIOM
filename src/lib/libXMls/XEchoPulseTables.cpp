
#include "libXMls/XEchoPulseTables.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

using namespace std;

XAttribTable::~XAttribTable()
{
    // todo: use shared_ptrs
    for(auto& it:m_map) delete it.second;
}

void XAttribTable::SetReferenceAttrib(std::string attrib_name)
{
    Map_t::iterator it = m_map.find(attrib_name);
    if(it == m_map.end())
        cout << "ERROR: setting as reference unexisting attribute " << attrib_name << endl;
    else mp_ref_attrib = it->second;
}

unsigned int XAttribTable::NSelectedBlocks()
{
    for(auto& it:m_map) if(it.second->m_is_accessible) return it.second->NBlock();
    return 0;
}

void XAttribTable::SetNBlocks(unsigned int n_blocks)
{
    for(auto& it:m_map) if(it.second->m_is_accessible) it.second->Resize(n_blocks);
}

unsigned int XAttribTable::BlockSize(XBlockIndex block_idx)
{
    for(auto& it:m_map)
    {
        unsigned int ret = it.second->NObjects(block_idx);
        if(ret > 0) return ret;
    }
    for(auto& it:m_map)
    {
        unsigned int ret = it.second->NObjectsFromFile(block_idx);
        if(ret > 0) return ret;
    }
    return 0;
}

void XAttribTable::MakeAccessible(XAbstractAttrib * p_attrib)
{
    //cout << "MakeAccessible " << p_attrib->m_attribname << endl;
    if(p_attrib->m_is_accessible) return;
    if(p_attrib->NBlock()>0)
        cout << "WARNING: called MakeAccessible() on non empty attribute " << p_attrib->m_attribname << endl;
    p_attrib->m_is_accessible = true; // always succeed
    if(mp_ref_attrib == NULL) return; // p_attrib is the first attribute so is the reference
    p_attrib->Resize(0);
    // add all selected blocks in the new attrib and load what is loaded
    for(XBlockIndex block_idx=0; block_idx<mp_ref_attrib->NBlock(); block_idx++)
    {
        XAbstractAttribBlock * ref_block = mp_ref_attrib->AbstractBlockPtr(block_idx);
        XBlockIndex new_block_idx = p_attrib->SelectBlock(ref_block->m_sec);
        if(new_block_idx != block_idx) cout << "ERROR: new_block_idx != block_idx" << endl;
        XAbstractAttribBlock * new_block = p_attrib->AbstractBlockPtr(block_idx);
        if(new_block->Size() == 0) new_block->SetSize(ref_block->Size());
        else if(new_block->Size() != ref_block->Size())
        {
            cout << "ERROR: incompatible block sizes " << new_block->Size() << "!=" <<
                    ref_block->Size() << ", something wrong with your block files" << endl;
        }
        if(ref_block->NLoadedObjects()>0) new_block->Load(); // loaded or created if no file
    }
}

void XAttribTable::MakeAllAccessible()
{
    for(auto& it:m_map)
        if(!it.second->m_is_accessible) MakeAccessible(it.second);
}

bool XAttribTable::AddAttrib(XAbstractAttrib * p_attrib, bool make_accessible)
{
    p_attrib->m_object = m_object_name;
    if(m_map.find(p_attrib->m_attribname) != m_map.end())
    {
        cout << "Warning: attempt to add already present attribute " << p_attrib->m_attribname << endl;
        cout << "Attribute will be loaded from block if present" << endl;
    }
    if(make_accessible) MakeAccessible(p_attrib);
    // default ept_path value is "."
    if(p_attrib->EptPath() == boost::filesystem::path(".") || p_attrib->EptPath().empty())
    {
        cout << "Taking ept path from ref attrib " << mp_ref_attrib->m_attribname
             << ":" << mp_ref_attrib->EptPath() << " for " << p_attrib->m_attribname << endl;
        p_attrib->SetEptPath(mp_ref_attrib->EptPath());
    }
    m_map[p_attrib->m_attribname] = p_attrib;
    return true;
}

// add a new table, size inferred from the others (loaded or not)
XAbstractAttrib * XAttribTable::AddAttrib(string attribname, string type_name, bool make_accessible)
{
    XAbstractAttrib * p_attrib = CreateAttribFromType(type_name);
    p_attrib->m_attribname = attribname;
    AddAttrib(p_attrib, make_accessible);
    return p_attrib;
}

void XAttribTable::AddAttribsFromFolder(boost::filesystem::path ept_path, bool make_accessible)
{
    using namespace boost::filesystem;
    if ( !exists( ept_path ) ) return;
    directory_iterator end_itr; // default construction yields past-the-end
    for ( directory_iterator itr( ept_path ); itr != end_itr; ++itr )
    {
        if ( is_directory(itr->status()) )
        {
            XAbstractAttrib * p_attrib = CreateAttribFromFolder(itr->path().leaf().string());
            if(p_attrib == NULL) continue; // other extension
            if(p_attrib->m_object != m_object_name) continue; // other object
            p_attrib->SetEptPath(ept_path);
            AddAttrib(p_attrib, make_accessible); // checks on block sizes done here
        }
    }
}

void XAttribTable::ListAttribs(std::vector<XAbstractAttrib *> & vp_attrib)
{
    for(auto& it: m_map) vp_attrib.push_back(it.second);
}

/// Select a time block. Note: you can select but not unselect (complicated, and do we really need it ?)
/// Number of objects in the block is estimated at selection (but block is not loaded/created)
void XAttribTable::Select(XSecond sec)
{
    XBlockIndex block_idx = 0;
    // select the block in all accessible attributes
    for(auto& it:m_map)
    {
        if(it.second->m_is_accessible)
        {
            block_idx = it.second->SelectBlock(sec);
        }
    }
    // TODO: check all are the same (they should be)
    size_t n_objects = ReferenceAttribute()->NObjectsFromFile(block_idx);
    for(auto& it:m_map) if(it.second->m_is_accessible)
    {
        XAbstractAttribBlock * p_block = it.second->AbstractBlockPtr(block_idx);
        size_t n_cur_objects = p_block->NObjectsFromFile();
        if(n_cur_objects>0 && n_cur_objects != n_objects) // we found a file (>0) and got a wrong number of objects from it
        {
            cout << "Inconsistent block size: " << n_cur_objects << "!=" << n_objects << endl;
        }
        p_block->SetSize(n_objects);
    }
}

/// Select time blocks in interval [sec_start, sec_end[
void XAttribTable::Select(XSecond sec_start, XSecond sec_end)
{
    for(XSecond sec = sec_start; sec < sec_end; sec++) Select(sec);
}

/// Load the i^th selected block for all accessible attibutes in the table
size_t XAttribTable::Load(XBlockIndex block_idx, boost::filesystem::path ept_path)
{
    size_t n_objects = 0;
    // load all you can and get number of objects
    for(auto& it:m_map)
    {
        if(it.second->m_is_accessible)
        {
            n_objects = it.second->Load(block_idx);
            //cout << "n_objects=" << n_objects << endl;
        }
    }
    return n_objects;
}
/// load all selected blocks
size_t XAttribTable::LoadSelected(boost::filesystem::path ept_path)
{
    size_t n_objects = 0;
    for(XBlockIndex block_idx=0; block_idx<NBlock(); block_idx++) n_objects += Load(block_idx, ept_path);
    return n_objects;
}
/// Selects and loads block sec
size_t XAttribTable::SelectAndLoad(XSecond sec, boost::filesystem::path ept_path)
{
    XBlockIndex block_idx=0;
    for(auto& it:m_map)
    {
        XBlockIndex cur_block_idx = it.second->SelectBlock(sec);
        if(block_idx == 0) block_idx = cur_block_idx;
        else if(block_idx != cur_block_idx)
            cout << "ERR: Inconsistent block indices " << block_idx << "!=" << cur_block_idx << endl;
    }
    return Load(block_idx, ept_path);
}
/// Selects and loads blocks in the time interval
size_t XAttribTable::SelectAndLoad(XSecond start_second, XSecond end_second, boost::filesystem::path ept_path)
{
    size_t n_objects = 0;
    for(XSecond sec = start_second; sec < end_second; sec++) n_objects += SelectAndLoad(sec, ept_path);
    return n_objects;
}

/// Free the i^th selected block
void XAttribTable::Free(XBlockIndex block_idx)
{
    for(auto& it:m_map) if(it.second->m_is_accessible) it.second->Free(block_idx);
}

/// Free all selected blocks
void XAttribTable::FreeSelected()
{
    for(auto& it:m_map) if(it.second->m_is_accessible) it.second->FreeSelected();
}

void XAttribTable::Save(boost::filesystem::path ept_path)
{
    for(auto& it:m_map) if(it.second->m_is_accessible) it.second->SaveSelected();
}

/// Select a time block (create it without loading)
/// Note: you can select but not unselect (complicated, and do we really need it ?)
void Select(XSecond sec);
/// Loads current selection
void LoadSelected();
/// Load the i^th selected block
void Load(XBlockIndex block_idx, boost::filesystem::path ept_path);
/// Load block sec
unsigned int Load(XSecond sec, boost::filesystem::path ept_path);
/// save all selected blocks
void Save(boost::filesystem::path ept_path);
/// Free the i^th selected block
void Free(XBlockIndex block_idx);
/// Free all selected blocks
void FreeSelected();


////////////////////////////////////////////////////    
int XPulseAttribTable::PulsePerLine()
{
    if(m_pulse_per_line) return m_pulse_per_line;
    Load(0); // load first block to estimate PulsePerLine
    // iteration on pulses
    float theta_0 = Theta(0,0), theta_jm1=theta_0, theta_j=theta_0;
    bool cycled = false;
    for(XPulseIndex j_pulse = 1; j_pulse < BlockSize(0); j_pulse++)
    {
        theta_j = Theta(0, j_pulse);
        if(theta_j < theta_jm1)
        {
            if(cycled) break; // cycled twice, happens if theta_0 is last before cycle
            cycled=true;
        }
        if(cycled && theta_j > theta_0)
        {
            m_pulse_per_line = j_pulse-1;
            return m_pulse_per_line;
        }
        theta_jm1 = theta_j;
    }
    // we should never reach this point
    cout << "ERROR: Not a full rotation in XPulseAttribTable" << endl;
    return m_pulse_per_line;
}

bool XPulseAttribTable::AddAttrib(XAbstractAttrib* p_attrib, bool make_accessible)
{
    bool force_accessible = true;
    if(p_attrib->m_attribname == "theta") mp_theta = dynamic_cast<XFloatAttrib *>(p_attrib);
    else if(p_attrib->m_attribname == "phi") mp_phi = dynamic_cast<XFloatAttrib *>(p_attrib);
    else if(p_attrib->m_attribname == "time") mp_time = dynamic_cast<XTAttrib<XLinearAttribBlock> *>(p_attrib);
    else if(p_attrib->m_attribname == "n_echo") mp_n_echo = dynamic_cast<XUCharAttrib *>(p_attrib);
    else if(p_attrib->m_attribname == "first_echo_index") mp_first_echo_index = dynamic_cast<XUIntAttrib *>(p_attrib);
    else force_accessible = false;
    return XAttribTable::AddAttrib(p_attrib, force_accessible || make_accessible);
}

// don't know why I can't use directly the inherited method
XAbstractAttrib * XPulseAttribTable::AddAttrib(string attribname, string type_name, bool make_accessible)
{
    return XAttribTable::AddAttrib(attribname, type_name, make_accessible);
}

////////////////////////////////////////////////////
bool XEchoAttribTable::AddAttrib(XAbstractAttrib* p_attrib, bool make_accessible)
{
    bool force_accessible = true;
    if(p_attrib->m_attribname == "range") mp_range = dynamic_cast<XFloatAttrib * >(p_attrib);
    else if(p_attrib->m_attribname == "pulse_index") mp_pulse_index = dynamic_cast<XUIntAttrib * >(p_attrib);
    else force_accessible = false; // physicals are not accessible by default, you must ask for them
    return XAttribTable::AddAttrib(p_attrib, force_accessible || make_accessible);
}
// don't know why I can't use directly the inherited method
XAbstractAttrib * XEchoAttribTable::AddAttrib(string attribname, string type_name, bool make_accessible)
{
    return XAttribTable::AddAttrib(attribname, type_name, make_accessible);
}

////////////////////////////////////////////////////
XEchoPulseTable::XEchoPulseTable(boost::filesystem::path ept_path, bool all_accessible)
{
    boost::filesystem::path info_path(ept_path);
    info_path /= "info.txt";
    ifstream ifs(info_path.c_str());
    if(!ifs.good())
    {
        cout << "ERR: XEchoPulseTable::XEchoPulseTable() cannot open " << info_path << endl;
        return;
    }
    string word1="", word2="", word3="";
    int delta;
    ifs >> word1 >> m_sec_min >> word2 >> m_sec_max >> word3 >> delta;
    cout << "Blocks are in range [" << m_sec_min << ", " << m_sec_max << "]" << endl;
    m_sec_max++; // consistency with load(min,max) loading in [min,max[
    if(word1!="Blocks") {cout << "ERR: First word is not Blocks in " << info_path << endl; return;}
    if(word2!="to") {cout << "ERR: second word is not to in " << info_path << endl; return;}
    if(word3!="=") {cout << "ERR: third word is not = in " << info_path << endl; return;}
    mts_pulse.AddAttribsFromFolder(ept_path);
    mts_echo.AddAttribsFromFolder(ept_path);
    // set reference attributes
    // !!! reference attributes are used to define block sizes !!!
    // other blocks should have the size of the reference or be absent (thus created at the correct size)
    ifs >> word1 >> word2 >> word3;
    if(word1!="pulse" || word2!="reference") {cout << "ERR: Second line should start with pulse reference in " << info_path << endl; return;}
    mts_pulse.SetReferenceAttrib(word3);
    ifs >> word1 >> word2 >> word3;
    if(word1!="echo" || word2!="reference") {cout << "ERR: Third line should start with echo reference in " << info_path << endl; return;}
    mts_echo.SetReferenceAttrib(word3);
    // add indexation attributes
    mts_pulse.AddAttrib("first_echo_index", "uint32");
    mts_echo.AddAttrib("pulse_index", "uint32");
    if(all_accessible)
    {
        mts_pulse.MakeAllAccessible();
        mts_echo.MakeAllAccessible();
    }
}

vector<XAbstractAttrib*> XEchoPulseTable::AttribList()
{
    vector<XAbstractAttrib*> vp_attrib;
    mts_pulse.ListAttribs(vp_attrib);
    mts_echo.ListAttribs(vp_attrib);
    return vp_attrib;
}

////////////////////////////////////////////////////
/// \brief XEchoPulseTable::AddAttribsFromFolder
/// \param ept_folder
void XEchoPulseTable::AddAttribsFromFolder(boost::filesystem::path ept_path)
{
    mts_pulse.AddAttribsFromFolder(ept_path);
    mts_echo.AddAttribsFromFolder(ept_path);
}

/// Select a time block. Note: you can select but not unselect (complicated, and do we really need it ?)
void XEchoPulseTable::Select(XSecond sec)
{
    if(sec == -1) sec = m_sec_min;
    mts_pulse.Select(sec);
    mts_echo.Select(sec);
}

/// Select time blocks in interval [sec_start, sec_end[
void XEchoPulseTable::Select(XSecond & sec_start, XSecond & sec_end)
{
    if(sec_start == -1) sec_start = m_sec_min;
    if(sec_end == -1) sec_end = m_sec_max;
    if(sec_start >= sec_end) return;
    if(sec_start < m_sec_min) sec_start = m_sec_min;
    if(sec_end > m_sec_max) sec_end = m_sec_max;
    mts_pulse.Select(sec_start, sec_end);
    mts_echo.Select(sec_start, sec_end);
}
/// Select all time blocks
void XEchoPulseTable::SelectAll()
{
    /// m_sec_max should be selected
    mts_pulse.Select(m_sec_min, m_sec_max);
    mts_echo.Select(m_sec_min, m_sec_max);
}
/// Load the i^th selected block
void XEchoPulseTable::Load(XBlockIndex block_idx)
{
    mts_pulse.Load(block_idx);
    mts_echo.Load(block_idx);
    EchoPulseIndexation(block_idx); // pulse<->echo links
}
/// Loads current selection
void XEchoPulseTable::LoadSelected()
{
    mts_pulse.LoadSelected();
    mts_echo.LoadSelected();
    // TODO: indexation
}
/// Selects and loads blocks in the time interval, start/end_second snapped to actual time interval
void XEchoPulseTable::SelectAndLoad(XSecond & sec_start, XSecond & sec_end)
{
    if(sec_start == -1) sec_start = m_sec_min;
    if(sec_end == -1) sec_end = m_sec_max;
    if(sec_start >= sec_end) return;
    if(sec_start < m_sec_min) sec_start = m_sec_min;
    if(sec_end < m_sec_max) sec_end = m_sec_max;
    mts_pulse.SelectAndLoad(sec_start, sec_end);
    mts_echo.SelectAndLoad(sec_start, sec_end);
}
/// Selects and loads everything (dangerous)
void XEchoPulseTable::SelectAndLoad()
{
    XSecond start_second=-1, end_second=-1;
    SelectAndLoad(start_second, end_second);
}
void XEchoPulseTable::Save(boost::filesystem::path ept_path)
{
    mts_pulse.Save(ept_path);
    mts_echo.Save(ept_path);
}
/// Free the i^th selected block
void XEchoPulseTable::Free(XBlockIndex block_idx)
{
    mts_pulse.Free(block_idx);
    mts_echo.Free(block_idx);
}
/// Free all selected blocks
void XEchoPulseTable::FreeSelected()
{
    mts_pulse.FreeSelected();
    mts_echo.FreeSelected();
}

void XEchoPulseTable::EchoPulseIndexation(XBlockIndex block_idx)
{
    XTAttribBlock<unsigned char> * p_n_echo_block = mts_pulse.mp_n_echo->BlockPtr(block_idx);
    XTAttribBlock<unsigned int> * p_first_echo_index_block = mts_pulse.mp_first_echo_index->BlockPtr(block_idx);
    XTAttribBlock<unsigned int> * p_pulse_index_block = mts_echo.mp_pulse_index->BlockPtr(block_idx);
    if(p_first_echo_index_block->FileExists() && p_pulse_index_block->FileExists()) return; // both attributes were read from disk
    //cout << p_n_echo_block->NLoadedObjects() << " pulses to link" << endl;
    //cout << "p_first_echo_index_block " << p_first_echo_index_block->NLoadedObjects() << endl;
    //cout << "p_pulse_index_block " << p_pulse_index_block->NLoadedObjects() << endl;
    XEchoIndex i_echo=0;
    for(XPulseIndex pulse_idx=0; pulse_idx < p_first_echo_index_block->NLoadedObjects(); pulse_idx++)
    {
        p_first_echo_index_block->at(pulse_idx) = i_echo;
        for(unsigned int i_echo_in_pulse=0; i_echo_in_pulse < p_n_echo_block->at(pulse_idx); i_echo_in_pulse++)
            p_pulse_index_block->at(i_echo++) = pulse_idx;
    }
}

void XEchoPulseTable::PrecomputeXYZ()
{
    // TODO
}

// access
// WARNING: no checks for performance but these will raise null pointer if corresponding attrib has not been loaded
XPt3D XEchoPulseTable::Ray(XBlockIndex block_idx, XPulseIndex pulse_idx)
{
    float theta = Theta(block_idx, pulse_idx);
    float phi = Phi(block_idx, pulse_idx);
    float cos_phi = cos(phi);
    return XPt3D(cos_phi*cos(theta), cos_phi*sin(theta), sin(phi));
}
float XEchoPulseTable::X(XBlockIndex block_idx, XEchoIndex echo_idx)
{
    return Range(block_idx, echo_idx)*cos(Theta(block_idx, echo_idx))*cos(Phi(block_idx, echo_idx));
}
float XEchoPulseTable::Y(XBlockIndex block_idx, XEchoIndex echo_idx)
{
    return Range(block_idx, echo_idx)*sin(Theta(block_idx, echo_idx))*cos(Phi(block_idx, echo_idx));
}
float XEchoPulseTable::Z(XBlockIndex block_idx, XEchoIndex echo_idx)
{
    return Range(block_idx,echo_idx)*sin(Phi(block_idx,echo_idx));
}
XPt3D XEchoPulseTable::P(XBlockIndex block_idx, XEchoIndex echo_idx)
{
    float theta = Theta(block_idx,echo_idx);
    float phi = Phi(block_idx,echo_idx);
    float r = Range(block_idx,echo_idx), rcp = r*cos(phi);
    return XPt3D(rcp*cos(theta), rcp*sin(theta), r*sin(phi));
}
