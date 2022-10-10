#pragma once

/// Mobile laser scan class consisting of
/// - A (multiecho) point cloud in sensor frame (XEchoPulseTable, inherited)
/// - An extrinsic calibration of the laser sensor (ins ref->sensor ref = XArchiRIEGLSensorXML, member)
/// - A trajectory for the vehicle ref (sbet = N Hz sampling of ins ref = trajecto_reader, member)
/// - a time shift between time expressed in the point cloud and the trajectory

#include "XEchoPulseTables.h"
#include "libXMls/XTrajecto.h"
//#include "libXBase/XArchiRIEGLSensorXML.h"
//#include "libXBase/XArchiRIEGLSensor.h"
#include "libXBase/XAbsoluteTime.h"

class XMls:public XEchoPulseTable
{
public:
    double m_time_shift; // time in secs to add to time in the point cloud to make it coherent with time in trajecto
    XTrajecto * mp_trajecto;
    XArchiGeoref m_sensor_georef;
    XAbsoluteTime m_ept_time_pivot, m_trajecto_time_pivot;
    struct Param{
        std::string ept_folder, laser_calib, sbet, acc;
        Param(){}
        Param(std::string ept_folder_, std::string laser_calib_,
              std::string sbet_, std::string acc_=""):
            ept_folder(ept_folder_), laser_calib(laser_calib_),
            sbet(sbet_), acc(acc_){}
    };
    XMls(std::string ept_folder, std::string laser_calib,
         XTrajecto * p_trajecto);

    /// select all the rxp (to avoid RAM issued, you should then load block by block)
    void Select();

    /// select a single 1 second block
    void Select(XSecond & start_second);

    /// select points and traj in time interval [start_second,end_second] (clamped to acquisition time interval)
    /// in seconds from rxp file (trajecto time is transformed)
    void Select(XSecond & start_second, XSecond & end_second);

    /// select points and traj in time interval [floor(start_time), ceil(end_time)]
    void Select(const XAbsoluteTime & start_time, const XAbsoluteTime & end_time);

    /// Echo coords in INS ref
    inline XPt3D Pins(XBlockIndex block_idx, XEchoIndex echo_idx){return m_sensor_georef.Applique_transfo(P(block_idx, echo_idx));}

    /// Beam origin coords in INS ref
    inline XPt3D Cins(){return m_sensor_georef.Translation();}

    /// Ray direction in INS ref
    inline XPt3D RayIns(XBlockIndex block_idx, XPulseIndex pulse_idx){return m_sensor_georef.Rotation()*Ray(block_idx, pulse_idx);}

    /// Ins ref in Lamb93
    XPt3D Oworld(double ins_time, bool precomputed=true);
    inline XPt3D Oworld(XBlockIndex block_idx, XPulseIndex pulse_idx, bool precomputed=true)
    {
        return Oworld(Time(block_idx, pulse_idx) + m_time_shift, precomputed);
    }
    inline XPt3D Oworld(XBlockIndex block_idx, XEchoIndex echo_idx, bool precomputed=true)
    {
        return Oworld(Time(block_idx, echo_idx) + m_time_shift, precomputed);
    }
    XArchiGeoref Ins(double ins_time, bool precomputed=true);
    inline XArchiGeoref Ins(XBlockIndex block_idx, XEchoIndex echo_idx, bool precomputed=true)
    {
        return Ins(Time(block_idx, echo_idx) + m_time_shift, precomputed);
    }
    inline XArchiGeoref Ins(XBlockIndex block_idx, XPulseIndex pulse_idx, bool precomputed=true)
    {
        return Ins(Time(block_idx, pulse_idx) + m_time_shift, precomputed);
    }

    /// Echo coords in Lamb93
    XPt3D Pworld(XBlockIndex block_idx, XEchoIndex echo_idx);
    XPt3D Pworld_interpol_frame(XBlockIndex block_idx, XEchoIndex echo_idx);
    XPt3D Pworld_interpol_angles(XBlockIndex block_idx, XEchoIndex echo_idx);
    XPt3D Pworld(XBlockIndex block_idx, XPulseIndex pulse_idx, double range);

    /// Beam origin coords in Lamb93
    XPt3D Cworld(XBlockIndex block_idx, XPulseIndex pulse_idx);
    XPt3D Cworld(XBlockIndex block_idx, XEchoIndex echo_idx);

    /// Ray direction in Lamb93
    XPt3D RayWorld(XBlockIndex block_idx, XPulseIndex pulse_idx);

    /// interpolated INS trajectory information at instant echo was acquired
    SbetEvent Sbet(XBlockIndex block_idx, XPulseIndex pulse_idx);
    SbetEvent Sbet(XBlockIndex block_idx, XEchoIndex echo_idx);

    /// interpolated INS accuracy information at instant echo was acquired
    AccuracyEvent Accuracy(XBlockIndex block_idx, XPulseIndex pulse_idx);
    AccuracyEvent Accuracy(XBlockIndex block_idx, XEchoIndex echo_idx);

};
