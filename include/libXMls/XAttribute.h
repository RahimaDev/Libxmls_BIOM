#pragma once

#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <typeinfo>
#include <boost/filesystem.hpp>

/// string typenames
template <typename T> std::string Type2Str()
{
    if(typeid(T) == typeid(float)) return "float32";
    if(typeid(T) == typeid(double)) return "float64";
    if(typeid(T) == typeid(int)) return "int32";
    if(typeid(T) == typeid(unsigned int)) return "uint32";
    if(typeid(T) == typeid(short)) return "int16";
    if(typeid(T) == typeid(unsigned short)) return "uint16";
    if(typeid(T) == typeid(char)) return "int8";
    if(typeid(T) == typeid(unsigned char)) return "uint8";
    if(typeid(T) == typeid(bool)) return "bool";
    return "unknown";
}

/// NamedType<T, Name> behaves like T, but has a name to ensure no confusion with another T
template<typename T, class Name> class NamedType
{
private:
    T value;
public:
    NamedType() { }
    NamedType (T v) : value(v) { }  // allow (*this) = v
    operator T & () { return value; }  // allow v = (*this)
    // type preserving addition
    template <typename U> NamedType<T, Name> Plus(U val) {return NamedType<T, Name>(value+val);}
};
// names for the types
struct EchoIndexName{};
struct PulseIndexName{};
struct BlockIndexName{};
struct SecondName{};

typedef NamedType<unsigned int, EchoIndexName> XEchoIndex;
typedef NamedType<unsigned int, PulseIndexName> XPulseIndex;
typedef NamedType<unsigned int, BlockIndexName> XBlockIndex;
typedef NamedType<int, SecondName> XSecond;

class XAbstractAttrib;

/** \brief Abstract attribute block class
*
* Abstract base class for all attributes blocks
* a member at(int i) should be defined in child classes to access i^th element
* but it cannot be in the mother class because return type depends on child
* An attribute block corresponds to a unique file on disk from which it is loaded an/or to which it is saved
* The name of this file is fixed by convention to $attrib_folder$/$m_sec$.$Ext()$
**/
class XAbstractAttribBlock
{
    /// number of objects in the block.
    /// If block is to be loaded from a file, the file size defines the number of objects.
    /// If block is created from scratch, this is set by the Resize() function
    size_t m_size;
public:
    /// pointer to the attribute, used to gather meta-info
    XAbstractAttrib * mp_attrib;
    /// second (identifier) of the block = integer part of the time value for all objects of this block
    XSecond m_sec;

    /// Ctor: needs a pointer to the attrib for meta-info access and the (second) identifier
    XAbstractAttribBlock(XAbstractAttrib * p_attrib, XSecond sec=0);

    /// name of the type of the attribute stored in the block
    virtual std::string Typename()=0;
    /// extention of files handling this type of attribute
    virtual std::string Ext()=0;
    /// name of the file where the attribute block is stored
    virtual std::string Filename();
    /// full path to the file corresponding to the block
    virtual boost::filesystem::path Path();
    /// create the folder containing the file (Path()) if it does not exist
    void CreateFolderIfNeeded();
    /// does the file corresponding to this block exist ?
    /// (true=attribute read from file, false=attribute created from scratch)
    inline bool FileExists(){return boost::filesystem::exists(Path());}
    /// load the attribute block from a file, if file does not exist, allocate memory for m_size objects
    /// Call NLoadedObjects() to know how many objects were loaded
    /// Call FileExists() to know if the block was loaded from a file or created from scratch
    virtual void Load()=0;
    /// save the attribute block to a file
    virtual void Save()=0;
    /// free the attribute (free the memory, you cannot access attribute values after a free)
    virtual void Free()=0;
    /// get the number of objects
    inline unsigned int Size(){return m_size;}
    /// set the number of objects (should only be called for attributes created from scratch)
    void SetSize(unsigned int size);
    /// use file to set the number of objects. Should be called by derived Ctors
    bool SetSizeFromFileIfExists();
    /// get the number of loaded objects (pulses or echos)
    /// use as a bool to know if the block is loaded or not
    virtual unsigned int NLoadedObjects()=0;
    /// get the number of objects from file size
    virtual unsigned int NObjectsFromFile()=0;
};

/** \brief typed attribute block class
*
* Implements an XAbstractAttribBlock as a vector of a common type.

**/
template <class T> class XTAttribBlock : public XAbstractAttribBlock, public std::vector<T>
{
public:
    typedef T value_type;
    XTAttribBlock(XAbstractAttrib * p_attrib, XSecond sec=0):XAbstractAttribBlock(p_attrib, sec) {SetSizeFromFileIfExists();}
    virtual std::string Typename(){return Type2Str<T>();}
    virtual std::string Ext(){return "bin";}
    virtual void Load()
    {
        if(NLoadedObjects()>0) return; // something is already loaded
        std::vector<T>::resize(Size());
        if(!FileExists())
        {
            //std::cout << "No file " << Path() << ", Load() simply allocates memory for n_objects=" << Size() << std::endl;
            return;
        }
        if(Size() != NObjectsFromFile())
        {
            std::cout << "ERR: Size()=" << Size() << " != NObjectsFromFile()= " << NObjectsFromFile() << std::endl;
            return;
        }
        std::ifstream ifs(Path().c_str(), std::ios_base::in | std::ios_base::binary);
        typename std::vector<T>::const_reference front = std::vector<T>::front();
        ifs.read((char*)&front, Size()*sizeof(T));
        ifs.close();
    }
    virtual void Save()
    {
        if(NLoadedObjects() == 0) return; // not loaded
        if(boost::filesystem::exists(Path()))
            std::cout << "WARN: save overwrites " << Path() << std::endl;
        CreateFolderIfNeeded();
        std::ofstream ofs(Path().c_str(), std::ios_base::out | std::ios_base::binary);
        if(!ofs.good())
        {
            std::cout << "ERR: cannot write to " << Path() << std::endl;
            return;
        }
        int filesize = std::vector<T>::size()*sizeof(T);
        typename std::vector<T>::const_reference front = std::vector<T>::front();
        ofs.write((char*)&front, filesize);
        ofs.close();
    }
    virtual void Free()
    {
        if(NLoadedObjects() == 0) return;
        std::vector<T>().swap(*this);
    }
    virtual unsigned int NLoadedObjects(){return std::vector<T>::size();}
    virtual unsigned int NObjectsFromFile()
    {
        if(!FileExists())
        {
            //std::cout << "NObjectsFromFile() " << Path() << " does not exist" << std::endl;
            return 0;
        }
        std::ifstream ifs(Path().c_str(), std::ios_base::in | std::ios_base::binary | std::ios::ate);
        unsigned int filesize = ifs.tellg();
        unsigned int n_objects = filesize/(sizeof(T));
        if(filesize != n_objects*sizeof(T))
        {
            std::cout << "ERR: " << Path() << " size " << filesize << " not multiple of type size " << sizeof(T) << std::endl;
            return 0;
        }
        //std::cout << Path() << " size " << filesize << " objects " << n_objects << std::endl;
        return n_objects;
    }
};

/** \brief lineat attribute block class
*
* Implements an XAbstractAttribBlock as a linear function of the index: a*index + b
**/
class XLinearAttribBlock : public XAbstractAttribBlock
{
    /// linear factors
    double m_a, m_b;
    /// because we do not really load for linear attribs, we need a bool to remember the state
    bool is_loaded;
public:
    XLinearAttribBlock(XAbstractAttrib * p_attrib, XSecond sec=0, double a=1., double b=0.):
        XAbstractAttribBlock(p_attrib, sec), m_a(a), m_b(b), is_loaded(false) {SetSizeFromFileIfExists();}
    virtual std::string Typename(){return "linear";}
    virtual std::string Ext(){return "txt";}
    virtual void Load();
    virtual void Save();
    virtual void Free();
    virtual unsigned int NLoadedObjects();
    virtual unsigned int NObjectsFromFile();
    template<class T> double at(T v){return m_a*v+m_b;}
};

/** \brief abstract attribute class
*
* An attribute corresponds to a subfolder of an ept_folder (which name follows a convention, cf AttributeFolder())
* and consists of a selection of attribute blocks (files in this folder)
* Selected blocks can be loaded/freed/saved
* If the file does not exist, the block is created from scratch and can be saved
*
**/
class XAbstractAttrib
{
public: // protected:
    /// default ept_path (from last load or save)
    boost::filesystem::path m_ept_path;
public:
    /// object to which the attribute is attached ("pulse" or "echo"), name of the attribute
    std::string m_object, m_attribname;
    /// only accessible attributes should be accessed.
    bool m_is_accessible;
    /// Ctor for attributes. default ept_path value is "." such that attributes created from scratch
    /// with no ept_path not provided will be saved in the CWD
    XAbstractAttrib(std::string object="", std::string attribname="", std::string ept_path="."):
        m_ept_path(ept_path), m_object(object), m_attribname(attribname), m_is_accessible(false) {}
    virtual ~XAbstractAttrib(){}

    // ---------- Pure virtuals --------- //
    /// name of the attribute type
    virtual std::string Typename()=0;
    /// name of the attribute extention
    virtual std::string Ext()=0;
    /// change the number of blocks
    virtual void Resize(unsigned int size)=0;
    /// number of blocks
    virtual unsigned int NBlock()=0;
    /// Pointer to the specified block
    virtual XAbstractAttribBlock * AbstractBlockPtr(XBlockIndex block_idx)=0;
    /// select a block (add it to the attribute but do not load the data)
    /// resizes to n_objects if given, returns the selected block index
    virtual XBlockIndex SelectBlock(XSecond sec, unsigned int n_objects=0)=0;
    virtual unsigned int BlockSizeFromFile(XBlockIndex block_idx)=0;

    // ---------- virtuals (filesystem) --------- //
    /// (relative) folder where the attribute is read/written
    virtual std::string AttributeFolder();
    /// name of the file where the attribute block of second sec is stored
    virtual std::string Filename(XSecond sec);
    virtual std::string Filename(XBlockIndex block_idx);
    /// set the ept path for this attribute. This is where the attribute folder is located.
    /// this folder will be used by load/save functions.
    /// If attribute is read from the disk this is the corresponding folder
    /// If attribute is added to attributes read from disk it is the same
    /// If attribute is created from scratch, this is empty and should be set before calling save
    inline void SetEptPath(boost::filesystem::path ept_path){m_ept_path = ept_path;}
    inline boost::filesystem::path EptPath(){return m_ept_path;}
    /// full path to this attribute's folder
    boost::filesystem::path AttributePath();
    /// full path to the file where the attribute block of second sec is stored
    virtual boost::filesystem::path BlockPath(XSecond sec);

    // ---------- inline Helpers --------- //
    /// Set the size of the specified block !!! Should only be done for attributes created from scratch !!!
    inline void SetBlockSize(XBlockIndex block_idx, unsigned int n_objects=0)
    {
        AbstractBlockPtr(block_idx)->SetSize(n_objects);
    }
    /// load the (selected) block
    inline unsigned int Load(XBlockIndex block_idx)
    {
        AbstractBlockPtr(block_idx)->Load();
        //std::cout << "XAbstractAttrib::Load(" << block_idx << ") " << m_attribname << " " << AbstractBlockPtr(block_idx)->NLoadedObjects() << std::endl;
        return AbstractBlockPtr(block_idx)->NLoadedObjects();
    }
    /// save the (selected) block
    inline void Save(XBlockIndex block_idx)
    {
        AbstractBlockPtr(block_idx)->Save();
    }
    /// save the loaded attribute blocks in the folder
    inline void SaveSelected()
    {
        for(XBlockIndex block_idx=0; block_idx<NBlock(); block_idx++) Save(block_idx);
    }
    /// free the (selected) block
    inline void Free(XBlockIndex block_idx)
    {
        if(m_is_accessible) AbstractBlockPtr(block_idx)->Free();
    }
    /// free selected time blocks
    inline void FreeSelected()
    {
        for(XBlockIndex block_idx=0; block_idx<NBlock(); block_idx++) Free(block_idx);
    }
    /// access to size (in memory) of specified block
    inline unsigned int NLoadedObjects(XBlockIndex block_idx){return AbstractBlockPtr(block_idx)->NLoadedObjects();}
    /// access to size of specified block from file size
    inline unsigned int NObjectsFromFile(XBlockIndex block_idx){return AbstractBlockPtr(block_idx)->NObjectsFromFile();}
    /// access to size of specified block from memory or file if not loaded
    inline unsigned int NObjects(XBlockIndex block_idx)
    {
        unsigned int n_loaded_objects = NLoadedObjects(block_idx);
        if(n_loaded_objects) return n_loaded_objects;
        return NObjectsFromFile(block_idx);
    }

    /// total number of loaded objects (sum over blocks)
    inline unsigned int NLoadedObjects()
    {
        unsigned int ret = 0;
        for(XBlockIndex block_idx = 0; block_idx < NBlock(); block_idx++) ret += NLoadedObjects(block_idx);
        return ret;
    }
    /// total number of objects (sum over blocks, loaded or not)
    inline unsigned int NObjects()
    {
        unsigned int ret = 0;
        for(XBlockIndex block_idx = 0; block_idx < NBlock(); block_idx++)
        {
            //std::cout << "Block " << block_idx << " size " << NObjects(block_idx) << std::endl;
            ret += NObjects(block_idx);
        }
        return ret;
    }
    /// get the second
    inline XSecond Second(XBlockIndex block_idx){return AbstractBlockPtr(block_idx)->m_sec;}
};

/** \brief template attribute class
*
* Implements XAbstractAttrib as a vector of attribute blocks (TAttribBlock should be a child of XAbstractAttribBlock)
*
**/
template <class TAttribBlock> class XTAttrib : public XAbstractAttrib, public std::vector<TAttribBlock>
{
public:
    XTAttrib(std::string object="", std::string attribname=""):
        XAbstractAttrib(object, attribname) {}

    /// name of the attribute type
    virtual std::string Typename(){TAttribBlock block(NULL); return block.Typename();}
    /// name of the attribute extention
    virtual std::string Ext(){TAttribBlock block(NULL); return block.Ext();}
    virtual void Resize(unsigned int size){std::vector<TAttribBlock>::resize(size, TAttribBlock(this));}
    virtual unsigned int NBlock(){return std::vector<TAttribBlock>::size();}
    /// get the block pointer with proper type
    TAttribBlock * BlockPtr(XBlockIndex block_idx)
    {
        return &(std::vector<TAttribBlock>::at(block_idx));
    }
    /// Pointer to the block_idx^th block
    virtual XAbstractAttribBlock * AbstractBlockPtr(XBlockIndex block_idx)
    {
        return (XAbstractAttribBlock *)BlockPtr(block_idx);
    }
    /// select time block at second sec. Optionally force the number of objects (dangerous)
    virtual XBlockIndex SelectBlock(XSecond sec, unsigned int n_objects=0)
    {
        XBlockIndex block_idx = std::vector<TAttribBlock>::size();
        //std::cout << "block_idx=" << block_idx << std::endl;
        std::vector<TAttribBlock>::push_back(TAttribBlock(this, sec));
        if(n_objects) std::vector<TAttribBlock>::at(block_idx).SetSize(n_objects);
        return block_idx;
    }
    virtual unsigned int BlockSizeFromFile(XBlockIndex block_idx)
    {
        return std::vector<TAttribBlock>::at(block_idx).NObjectsFromFile();
    }
};

/** \brief template attribute class
*
* Shifted attributes are float attribs with a pivot point to handle large values (use for geographic coords or time)
*
**/
/*class XPivotAttrib : public XTAttrib< XTAttribBlock<float> >
{
public:
    XTAttrib(std::string object="", std::string attribname=""):
        XAbstractAttrib(object, attribname) {}

    double at(XBlockIndex block_idx){return }
};*/

/// Factory for attributes from a type name
XAbstractAttrib * CreateAttribFromType(std::string type_name);
/// Factory for attributes from an attribute folder name
XAbstractAttrib * CreateAttribFromFolder(std::string folder);

/// attribute pairs
template <class TAttribBlock> class X2TAttrib
{
    XTAttrib<TAttribBlock> * m_a1, * m_a2;
    X2TAttrib(XTAttrib<TAttribBlock> * a1, XTAttrib<TAttribBlock> * a2):m_a1(a1), m_a2(a2){}
    X2TAttrib(std::string object="", std::string attrib1_name="", std::string attrib2_name="")
    {
        m_a1=new XTAttrib<TAttribBlock>(object, attrib1_name);
        m_a2=new XTAttrib<TAttribBlock>(object, attrib2_name);
    }
};

/// attribute triplets
template <class TAttribBlock> class X3TAttrib
{
    XTAttrib<TAttribBlock> *m_a1, *m_a2, *m_a3;
    X3TAttrib(XTAttrib<TAttribBlock> a1, XTAttrib<TAttribBlock> & a2, XTAttrib<TAttribBlock> & a3):m_a1(a1), m_a2(a2), m_a3(a3){}
    X3TAttrib(std::string object="", std::string attrib1_name="", std::string attrib2_name="", std::string attrib3_name="")
    {
        m_a1=new XTAttrib<TAttribBlock>(object, attrib1_name);
        m_a2=new XTAttrib<TAttribBlock>(object, attrib2_name);
        m_a3=new XTAttrib<TAttribBlock>(object, attrib3_name);
    }
};

/// typedef standard attrib types
typedef XTAttrib<XTAttribBlock<double> > XDoubleAttrib; // time
typedef XTAttrib<XLinearAttribBlock> XLinearAttrib; // time (linear approx)
typedef XTAttrib<XTAttribBlock<float> > XFloatAttrib; // range, amplitude, reflectance
typedef X2TAttrib<XTAttribBlock<float> > X2FloatAttrib; // (theta, phi)
typedef X3TAttrib<XTAttribBlock<float> > X3FloatAttrib; // XYZ(_sensor), normal, dimentionality
typedef XTAttrib<XTAttribBlock<unsigned int> > XUIntAttrib; // index, class, id
typedef XTAttrib<XTAttribBlock<int> > XIntAttrib;
typedef XTAttrib<XTAttribBlock<unsigned short> > XUShortAttrib;
typedef XTAttrib<XTAttribBlock<short> > XShortAttrib;
typedef XTAttrib<XTAttribBlock<unsigned char> > XUCharAttrib; // nb_of_echo, num_echo, deviation
typedef X3TAttrib<XTAttribBlock<unsigned char> > X3UCharAttrib; // RGB
typedef XTAttrib<XTAttribBlock<char> > XCharAttrib;
typedef XTAttrib<XTAttribBlock<bool> > XBoolAttrib; // pruning
