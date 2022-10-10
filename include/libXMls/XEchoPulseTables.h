#pragma once

/// Echo/pulse tables store attibutes per object (echo or pulse) with a block structure to handle large point clouds
/// each attribute is stored in an attribute_folder which is a sub folder of an ept_folder
/// An attribute_folder contains one .bin or .txt file for each block of one second, which name is the (rounded) second itself
/// Blocks, Pulses and Echos are indexed by XBlockIndex, XPulseIndex and XEchoIndex which are typed integers.
/// Supported attribute types are:
/// -vectors of common types (stored in .bin)
/// -linear attributes (value is computed linearily from pulse_index, used mainly for time)
/// Philosophy:
///  * load only what you need (attributes AND time blocks).
///  * Except for "structural" attributes that are required to be always accessible:
///   + XLinearAttrib time;
///   + XFloatAttrib theta;
///   + XFloatAttrib phi;
///   + XUCharAttrib n_echo;
///   + XFloatAttrib range;
///   + XUIntAttrib first_echo_index; (computed from n_echo)
///   + XUIntAttrib pulse_index; (computed from n_echo)
/// Standard usage:
/// -Select the time blocks that you want to work with (does not load).
///  Each selected time block can be loaded (from a file), created (from scratch), freed (free corrresponding memory) and saved (to a file)
/// -Make the attributes that you want to work with accessible. Only accessible attributes are loaded/created/freed/saved.
///  * Structural attributes are always accessible (time, theta, phi, n_echo, range) and have direct accessors
///    ex: ep_table.Range(block_idx, echo_idx) (see end of this file)
///  * If a (non structural) attribute is present in the ept_folder, you can access it with
///    XFloatAttrib * p_attrib = ep_table.GetEchoAttrib<XFloatAttrib>("attrib_name");
///  * To create and access a new attribute, use:
///    XUCharAttrib * p_attrib = ep_table.AddPulseAttrib<XUCharAttrib>("attrib_name");
///  The two syntaxes above make the attributes accessible, meaning they have the same selected/loaded blocks pattern
///  so attribute values for loaded blocks can be accessed with p_attrib->at(block_idx)[echo_idx]
/// -XEchoPulseTable support attributes in multiple folders (for instance, exported attributes on a read only server and computed attributes locally)
///  * List the present attributes with XEchoPulseTable::AddAttribsFromFolder(boost::filesystem::path ept_path);
///  * Save an attribute in any folder with
///   + p_attrib->Save(XBlockIndex block_idx, boost::filesystem::path ept_path="")
///   + p_attrib->SaveSelected(boost::filesystem::path ept_path="")

#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <iostream>
#include <typeinfo>
#include <boost/filesystem.hpp>
#include "libXBase/XPt3D.h"
#include "XAttribute.h"

/** \brief Attribute Table class
*
* Table of attributes refering to the same object (pulse or echo)
* If you work with single echo, use directly this class, with object=echo
*
**/
class XAttribTable
{
public:    
    /// name of the object to which the attribute refers (pulse, echo)
    std::string m_object_name;

    /// type of attribute map. todo: use shared_ptrs
    typedef std::map<std::string, XAbstractAttrib*> Map_t;
private:
    /// map storing all attributes relative to this object
    Map_t m_map;
    XAbstractAttrib * mp_ref_attrib;

public:
    XAttribTable(std::string object_name):m_object_name(object_name),mp_ref_attrib(NULL){}
    ~XAttribTable();
    /// Accessible attributes have NBlock blocks that can be loaded/freed/saved/created(=allocated)
    /// Non accessible attributes have no block (just meta-info)
    inline XAbstractAttrib * ReferenceAttribute() {return mp_ref_attrib;}
    void SetReferenceAttrib(std::string attrib_name);
    virtual unsigned int NSelectedBlocks();
    /// Resize all accessible attributes
    virtual void SetNBlocks(unsigned int size);
    /// size of the block from memory (if loaded) or file
    unsigned int BlockSize(XBlockIndex block_idx);
    /// give the same pattern (number of blocks and of objects per block) to p_attrib as the FirstAccessibleAttribute()
    void MakeAccessible(XAbstractAttrib * p_attrib);
    /// make all added attributes accessible
    void MakeAllAccessible();
    /// add a (non accessible, else it is made non accessible and cleared) existing attrib
    virtual bool AddAttrib(XAbstractAttrib * p_attrib, bool make_accessible=true);
    /// add a new (non accessible) abstract attrib
    virtual XAbstractAttrib * AddAttrib(std::string attribname, std::string type_name, bool make_accessible=true);
    /// add a new (non accessible) typed attrib
    template <class T> T * AddAttrib(std::string attribname, bool make_accessible=true)
    {
        T * p_attrib = new T(m_object_name, attribname);
        AddAttrib(p_attrib, make_accessible);
        return p_attrib;
    }
    /// add (non accessible=meta-info only) attributes from folder
    virtual void AddAttribsFromFolder(boost::filesystem::path ept_path, bool make_accessible=false);

    inline Map_t AttribMap(){return m_map;}
    void ListAttribs(std::vector<XAbstractAttrib *> & vp_attrib);
    /// Get the attribute with the specified name and type and make it accessible except if not required
    template <class T> T * GetAttrib(const std::string & attribname, bool make_accessible=true)
    {
        Map_t::iterator it = m_map.find(attribname);
        if(it == m_map.end())
        {
            std::cout << "Warning: no attribute found with name " << attribname << std::endl;
            return NULL;
        }
        T * p_attrib = dynamic_cast<T *>(it->second);
        if(p_attrib == NULL)
        {
            std::cout << "Warning: asked for wrong type for attribute " << attribname << std::endl;
        }
        if(make_accessible) MakeAccessible(p_attrib); // if user wants it its to use it
        return p_attrib;
    }
    inline unsigned int NAttrib(){return m_map.size();}
    inline unsigned int NBlock(){return mp_ref_attrib->NBlock();}
    inline unsigned int NLoadedObjects(){return mp_ref_attrib->NLoadedObjects();}
    inline unsigned int NObjects(){return mp_ref_attrib->NObjects();}
    inline unsigned int NObjects(XBlockIndex block_idx){return mp_ref_attrib->NObjects(block_idx);}
    inline XSecond Second(XBlockIndex block_idx){return mp_ref_attrib->Second(block_idx);}

    /// Select a time block (create it without loading)
    /// Note: you can select but not unselect (complicated, and do we really need it ?)
    void Select(XSecond sec);
    /// Select time blocks in interval [sec_start, sec_end[
    void Select(XSecond sec_start, XSecond sec_end);

    /// Load the i^th selected block for all accessible attibutes in the table
    size_t Load(XBlockIndex block_idx, boost::filesystem::path ept_path="");
    /// load all selected blocks
    size_t LoadSelected(boost::filesystem::path ept_path="");
    /// Selects and loads block sec
    size_t SelectAndLoad(XSecond sec, boost::filesystem::path ept_path="");
    /// Selects and loads blocks in the time interval
    size_t SelectAndLoad(XSecond start_second, XSecond end_second, boost::filesystem::path ept_path="");
    /// save all selected blocks
    void Save(boost::filesystem::path ept_path="");
    /// Free the i^th selected block
    void Free(XBlockIndex block_idx);
    /// Free all selected blocks
    void FreeSelected();
};

/** \brief Pulse Attribute Table
*
* Table of pulse attributes
* Mainly shortcuts to structural pulse attributes
**/
class XPulseAttribTable:public XAttribTable
{
public:
    /// shortcuts to structural attribs
    XLinearAttrib * mp_time;
    XFloatAttrib * mp_theta;
    XFloatAttrib * mp_phi;
    XUCharAttrib * mp_n_echo;
    XUIntAttrib * mp_first_echo_index;
    int m_pulse_per_line;

    XPulseAttribTable():XAttribTable("pulse"),m_pulse_per_line(0){}

    /// Number of pulse per line, rounded towards 0.
    /// The 6 neighbors of i in sensor topology are
    /// i-PulsePerLine-1, i-PulsePerLine, i-1, i+1, i+PulsePerLine, i+PulsePerLine+1
    int PulsePerLine();

    /// adding an attribute to the table
    virtual bool AddAttrib(XAbstractAttrib* p_attrib, bool make_accessible=true);
    /// adding an attribute to the table from its name and type
    virtual XAbstractAttrib * AddAttrib(std::string attribname, std::string type_name, bool make_accessible=true);
    /// add a new (non accessible) typed attrib
    template <class T> T * AddAttrib(std::string attribname, bool make_accessible=true)
    {
        T * p_attrib = new T(m_object_name, attribname);
        AddAttrib(p_attrib, make_accessible);
        return p_attrib;
    }

    /// shortcuts
    inline float Theta(XBlockIndex block_idx, XPulseIndex pulse_idx){return mp_theta->at(block_idx).at(pulse_idx);}
    inline float Phi(XBlockIndex block_idx, XPulseIndex pulse_idx){return mp_phi->at(block_idx).at(pulse_idx);}
    inline float NbOfEcho(XBlockIndex block_idx, XPulseIndex pulse_idx){return mp_n_echo->at(block_idx).at(pulse_idx);}
    inline int FirstEchoIdx(XBlockIndex block_idx, XPulseIndex pulse_idx){return mp_first_echo_index->at(block_idx).at(pulse_idx);}
};

/** \brief Echo Attribute Table
*
* Table of echo attributes
* Mainly shortcuts to structural echo attributes
**/
class XEchoAttribTable:public XAttribTable
{
public:
    /// shortcuts to structural attribs
    XFloatAttrib * mp_range;
    XUIntAttrib * mp_pulse_index;

    XEchoAttribTable():XAttribTable("echo"){}

    virtual bool AddAttrib(XAbstractAttrib* p_attrib, bool make_accessible=true);
    virtual XAbstractAttrib * AddAttrib(std::string attribname, std::string type_name, bool make_accessible=true);
    /// add a new (non accessible) typed attrib
    template <class T> T * AddAttrib(std::string attribname, bool make_accessible=true)
    {
        T * p_attrib = new T(m_object_name, attribname);
        AddAttrib(p_attrib, make_accessible);
        return p_attrib;
    }

    // shortcuts
    inline double Range(XBlockIndex block_idx, XEchoIndex echo_idx){return mp_range->at(block_idx).at(echo_idx);}
    inline unsigned int PulseIdx(XBlockIndex block_idx, XEchoIndex echo_idx){return mp_pulse_index->at(block_idx).at(echo_idx);}
};


/** \brief Echo and Pulse Tables
*
* Aggregate the echo and pulse tables for a point cloud
**/
class XEchoPulseTable
{
public:
    XPulseAttribTable mts_pulse;
    XEchoAttribTable mts_echo;

    /// metainfos
    XSecond m_sec_min, m_sec_max; // first and last exported block

    //XEchoPulseTable(){}
    /// creation: need an ept_folder with info.txt
    XEchoPulseTable(boost::filesystem::path ept_path=".", bool all_accessible=false);

    /// get all attributes in a vector
    std::vector<XAbstractAttrib*> AttribList();

    /// add attributes from another folder (user attributes)
    /// if current attributes have loaded blocks, loads the same
    void AddAttribsFromFolder(boost::filesystem::path ept_path);
    /// Select time block of index sec
    void Select(XSecond sec);
    /// Select time blocks in interval [sec_start, sec_end[ (-1 = bound)
    void Select(XSecond & sec_start, XSecond & sec_end);
    /// Select all time blocks
    void SelectAll();
    /// Loads current selection
    void LoadSelected();
    /// Selects and loads everything (dangerous)
    void SelectAndLoad();
    /// Selects and loads blocks in the time interval, start/end_second snapped to actual time interval
    /// Pass -1 to get bounds
    void SelectAndLoad(XSecond & sec_start, XSecond & sec_end);
    /// Load the i^th selected block
    void Load(XBlockIndex block_idx);
    /// Free the i^th selected block
    void Free(XBlockIndex block_idx);
    /// Free all selected blocks
    void FreeSelected();

    /// Save all attributes in the specified folder
    /// You should never need this function but only save new attributes instead
    /// Only interest is if you create an XEchoPulseTable from scratch (good luck)
    void Save(boost::filesystem::path ept_path);
    /// iterate on blocks with for(XBlockIndex block_idx=0; block_idx<NBlock(); block_idx++)
    inline XBlockIndex NBlock(){return mts_pulse.NBlock();}
    inline unsigned int PulsePerLine(){return mts_pulse.PulsePerLine();}

    void AddEchoAttrib(XAbstractAttrib* p_attrib) {mts_echo.AddAttrib(p_attrib);}
    void AddPulseAttrib(XAbstractAttrib* p_attrib) {mts_pulse.AddAttrib(p_attrib);}
    inline XAbstractAttrib * AddEchoAttrib(std::string attribname, std::string type_name, bool make_accessible=true)
    {
        return mts_echo.AddAttrib(attribname, type_name, make_accessible);
    }
    inline XAbstractAttrib * AddPulseAttrib(std::string attribname, std::string type_name, bool make_accessible=true)
    {
        return mts_pulse.AddAttrib(attribname, type_name, make_accessible);
    }
    template <class T> T * AddEchoAttrib(std::string attribname, bool make_accessible=true)
    {
        return mts_echo.AddAttrib<T>(attribname, make_accessible);
    }
    template <class T> T * AddPulseAttrib(std::string attribname, bool make_accessible=true)
    {
        return mts_pulse.AddAttrib<T>(attribname, make_accessible);
    }
    template <class T> T * GetPulseAttrib(const std::string & attribname, bool load=true)
    {
        return mts_pulse.GetAttrib<T>(attribname, load);
    }
    template <class T> T * GetEchoAttrib(const std::string & attribname, bool load=true)
    {
        return mts_echo.GetAttrib<T>(attribname, load);
    }

    inline int NEchoAttrib(){return mts_echo.NAttrib();}
    inline int NPulseAttrib(){return mts_pulse.NAttrib();}
    inline XEchoIndex NEcho(XBlockIndex block_idx){return mts_echo.NObjects(block_idx);}
    inline XPulseIndex NPulse(XBlockIndex block_idx){return mts_pulse.NObjects(block_idx);}
    inline unsigned int NTotalEcho(){return mts_echo.NObjects();}
    inline unsigned int NTotalPulse(){return mts_pulse.NObjects();}
    inline XSecond Second(XBlockIndex block_idx){if(mts_echo.NBlock()) return mts_echo.Second(block_idx); return 0;}
    inline XSecond FirstSecond(){return Second(0);}
    inline XSecond LastSecond(){return Second(NBlock()-1);}

    /// echo/pulse indexation (indices stored in echo/pulse attributes)
    void EchoPulseIndexation(XBlockIndex block_idx);
    /// XYZ coordinates precomputation
    void PrecomputeXYZ();

    // ACCESSORS. WARNING: no checks for performance but these will raise null pointer if corresponding attrib has not been loaded
    // echo<->pulse
    inline XPulseIndex IdxPulse(XBlockIndex block_idx, XEchoIndex echo_idx){return XPulseIndex(mts_echo.mp_pulse_index->at(block_idx).at(echo_idx));}
    inline int NumEcho(XBlockIndex block_idx, XEchoIndex echo_idx){return echo_idx - IdxFirstEcho(block_idx, IdxPulse(block_idx, echo_idx));}
    inline int NbOfEcho(XBlockIndex block_idx, XPulseIndex pulse_idx){return mts_pulse.mp_n_echo->at(block_idx).at(pulse_idx);}
    inline int NbOfEcho(XBlockIndex block_idx, XEchoIndex echo_idx){return NbOfEcho(block_idx, IdxPulse(block_idx, echo_idx));}
    inline XEchoIndex IdxFirstEcho(XBlockIndex block_idx, XPulseIndex pulse_idx){return XEchoIndex(mts_pulse.mp_first_echo_index->at(block_idx).at(pulse_idx));}
    inline XEchoIndex IdxLastEcho(XBlockIndex block_idx, XPulseIndex pulse_idx){return XEchoIndex(IdxFirstEcho(block_idx, pulse_idx) + NbOfEcho(block_idx, pulse_idx) - 1);}

    // pulse specific access
    inline double Time(XBlockIndex block_idx, XPulseIndex pulse_idx){return mts_pulse.mp_time->at(block_idx).at(pulse_idx);}
    inline double Time(XBlockIndex block_idx, XEchoIndex echo_idx){return Time(block_idx, IdxPulse(block_idx, echo_idx));}
    inline float Theta(XBlockIndex block_idx, XPulseIndex pulse_idx){return mts_pulse.mp_theta->at(block_idx).at(pulse_idx);}
    inline float Theta(XBlockIndex block_idx, XEchoIndex echo_idx){return Theta(block_idx, IdxPulse(block_idx, echo_idx));}
    inline float Phi(XBlockIndex block_idx, XPulseIndex pulse_idx){return mts_pulse.mp_phi->at(block_idx).at(pulse_idx);}
    inline float Phi(XBlockIndex block_idx, XEchoIndex echo_idx){return Phi(block_idx, IdxPulse(block_idx, echo_idx));}
    XPt3D Ray(XBlockIndex block_idx, XPulseIndex pulse_idx);

    // echo specific access
    inline float Range(XBlockIndex block_idx, XEchoIndex echo_idx){return mts_echo.mp_range->at(block_idx).at(echo_idx);}
    inline unsigned int PulseIdx(XBlockIndex block_idx, XEchoIndex echo_idx){return mts_echo.mp_pulse_index->at(block_idx).at(echo_idx);}

    // geometry access
    float X(XBlockIndex block_idx, XEchoIndex echo_idx);
    float Y(XBlockIndex block_idx, XEchoIndex echo_idx);
    float Z(XBlockIndex block_idx, XEchoIndex echo_idx);
    XPt3D P(XBlockIndex block_idx, XEchoIndex echo_idx);

};

