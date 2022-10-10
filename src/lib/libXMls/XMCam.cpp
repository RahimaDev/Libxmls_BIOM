
#include "libXMls/XMCam.h"
#include "libXBase/XArchiXMLBaseTools.h"
#include "libXBase/XArchiXMLTools.h"
#include "libXBase/XArchiXMLException.h"
#include "libXBase/XArchiGeorefXML.h"
#include "tinyxml.h"
#include <boost/filesystem.hpp>

using namespace std;

XMCam::XMCam(std::string instrinsic_file, std::string extrinsic_file,
             std::string trigger_file, XTrajecto * p_trajecto):
    XCam(instrinsic_file),
    mp_trajecto(p_trajecto)
{
    // TODO: read trigger_file
    cout << "Reading " << instrinsic_file << endl;

    // read extrinsic calib
    XArchiXML::XArchiGeoref_LoadFromNode(&m_camera_georef, extrinsic_file);

    // read time pivots
    m_trigger_time_pivot.XmlLoad(trigger_file+".pvt.xml");

    boost::filesystem::path sbet_pivot_name = p_trajecto->Filename() + ".pvt.xml";
    m_trajecto_time_pivot.XmlLoad(sbet_pivot_name.string());

    // compute shift
    m_time_shift = m_trigger_time_pivot.Diff(m_trajecto_time_pivot);
}
