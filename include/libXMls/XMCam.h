#pragma once

/// Mobile camera class consisting of
/// - A camera (XCam) known by its intrinsic parameters
/// - An extrinsic calibration of the camera (ins ref->camera ref = XArchiRIEGLSensorXML, member)
/// - A trajectory for the vehicle ref (sbet = N Hz sampling of ins ref = trajecto_reader, member)
/// - A list of trigger times when the camera was actionned
/// - a time shift between time expressed in the camera trigs and the trajectory

#include "libXMls/XCam.h"
#include "libXMls/XTrajecto.h"
#include "libXBase/XAbsoluteTime.h"

class XMCam:public XCam
{
public:
    double m_time_shift; // time in secs to add to time in the point cloud to make it coherent with time in trajecto
    XTrajecto * mp_trajecto;
    XArchiGeoref m_camera_georef;
    XAbsoluteTime m_trigger_time_pivot, m_trajecto_time_pivot;

    XMCam(std::string instrinsic_file, std::string extrinsic_file,
         std::string trigger_file, XTrajecto * p_trajecto);

};
