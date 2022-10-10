
#include "libXMls/XMls.h"
#include "libXBase/XArchiGeorefXML.h"
#include <boost/filesystem.hpp>

using namespace std;

XMls::XMls(std::string ept_folder, std::string laser_calib,
           XTrajecto * p_trajecto):
    XEchoPulseTable(ept_folder),
    mp_trajecto(p_trajecto)
{
    // read laser calib
    XArchiXML::XArchiGeoref_LoadFromNode(&m_sensor_georef, laser_calib);

    // read time pivots
    m_ept_time_pivot.XmlLoad(ept_folder+".pvt.xml");

    boost::filesystem::path sbet_pivot_name = p_trajecto->Filename() + ".pvt.xml";
    m_trajecto_time_pivot.XmlLoad(sbet_pivot_name.string());

    // compute shift
    m_time_shift = m_ept_time_pivot.Diff(m_trajecto_time_pivot);
}

// in seconds from rxp file, trajecto time is transformed
void XMls::Select()
{
    XSecond start_second=-1, end_second=-1;
    Select(start_second, end_second);
}

void XMls::Select(XSecond & start_second)
{
    XSecond end_second=start_second+1;
    Select(start_second, end_second);
}

// in seconds from rxp file, trajecto time is transformed
void XMls::Select(XSecond & start_second, XSecond & end_second)
{
    // load echos
    XEchoPulseTable::Select(start_second, end_second);
    cout << "Select (" << start_second << "," << end_second << ")+" << m_time_shift << endl;

    // load trajecto
    mp_trajecto->Load(start_second+m_time_shift, end_second+m_time_shift); // todo: look for all trajectory parts of the same trajectory
    mp_trajecto->PreComputeGeoref(); // todo: optionnally change coord system, precompute with laser calib
}

// in absolute UTC time
void XMls::Select(const XAbsoluteTime & start_time, const XAbsoluteTime & end_time)
{
    XSecond start_second=floor(start_time.Diff(m_ept_time_pivot));
    XSecond end_second=ceil(end_time.Diff(m_ept_time_pivot));
    cout << "start_second " << start_second << endl;
    cout << "end_second " << end_second << endl;
    Select(start_second, end_second);
}

XPt3D XMls::Oworld(double ins_time, bool precomputed)
{
    XPt3D O;
    bool OK = mp_trajecto->GetTranslation(ins_time, O, precomputed);
    if(!OK)
        cout << "Warning: Oworld(" << ins_time << ") outside bounds " << mp_trajecto->StartTime() << "," << mp_trajecto->EndTime() << endl;
    return O;
}

XArchiGeoref XMls::Ins(double ins_time, bool precomputed)
{
    XArchiGeoref sbet_georef;
    bool OK = false;
    if(precomputed)
        OK = mp_trajecto->GetGeoref_precomputed(ins_time, sbet_georef);
    else OK = mp_trajecto->GetGeoref(ins_time, sbet_georef);
    if(!OK)
        cout << "Warning: Ins(" << ins_time << ") outside bounds " << mp_trajecto->StartTime() << "," << mp_trajecto->EndTime() << endl;
    return sbet_georef;
}

XPt3D XMls::Pworld(XBlockIndex block_idx, XEchoIndex echo_idx)
{
    return mp_trajecto->ApplyGeoref_precomputed(Time(block_idx, echo_idx) + m_time_shift, Pins(block_idx, echo_idx));
}

XPt3D XMls::Pworld(XBlockIndex block_idx, XPulseIndex pulse_idx, double range)
{
    return Cworld(block_idx, pulse_idx) + range * RayWorld(block_idx, pulse_idx);
}

XPt3D XMls::Pworld_interpol_frame(XBlockIndex block_idx, XEchoIndex echo_idx)
{
    XArchiGeoref sbet_georef = Ins(block_idx, echo_idx, true);
    return sbet_georef.Applique_transfo(Pins(block_idx, echo_idx));
}

XPt3D XMls::Pworld_interpol_angles(XBlockIndex block_idx, XEchoIndex echo_idx)
{
    XArchiGeoref sbet_georef = Ins(block_idx, echo_idx, false);
    return sbet_georef.Applique_transfo(Pins(block_idx, echo_idx));
}

XPt3D XMls::Cworld(XBlockIndex block_idx, XEchoIndex echo_idx)
{
    XArchiGeoref sbet_georef = Ins(block_idx, echo_idx);
    return sbet_georef.Applique_transfo(Cins());
}

XPt3D XMls::Cworld(XBlockIndex block_idx, XPulseIndex pulse_idx)
{
    XArchiGeoref sbet_georef = Ins(block_idx, pulse_idx);
    return sbet_georef.Applique_transfo(Cins());
}

/// Ray direction in world coords
XPt3D XMls::RayWorld(XBlockIndex block_idx, XPulseIndex pulse_idx)
{
    XArchiGeoref sbet_georef = Ins(block_idx, pulse_idx);
    return sbet_georef.Rotation()*RayIns(block_idx, pulse_idx);
}

SbetEvent XMls::Sbet(XBlockIndex block_idx, XPulseIndex pulse_idx)
{
    SbetEvent sbet_event;
    double trajecto_time = Time(block_idx, pulse_idx) + m_time_shift;
    //TEventSeries<SbetEvent> * p_event_series = dynamic_cast<TEventSeries<SbetEvent> *>(mp_trajecto);
    if(!mp_trajecto->Interpol_event(sbet_event, trajecto_time))
        cout << "Warning: Interpol_event outside bounds for time " << trajecto_time << endl;
    return sbet_event;
}

SbetEvent XMls::Sbet(XBlockIndex block_idx, XEchoIndex echo_idx)
{
    SbetEvent sbet_event;
    double trajecto_time = Time(block_idx, echo_idx) + m_time_shift;
    if(!mp_trajecto->Interpol_event(sbet_event, trajecto_time))
        cout << "Warning: Interpol_event outside bounds for time " << trajecto_time << endl;
    return sbet_event;
}

AccuracyEvent XMls::Accuracy(XBlockIndex block_idx, XPulseIndex pulse_idx)
{
    AccuracyEvent acc_event;
    if(!mp_trajecto->HasAcc())
    {
        cout << "Warning: Asked for an accuracy but no accuracy file provided, "
             << "returning uninitialized accuracy event" << endl;
        return acc_event;
    }
    double trajecto_time = Time(block_idx, pulse_idx) + m_time_shift;
    if(!mp_trajecto->Accuracy().Interpol_event(acc_event, trajecto_time))
        cout << "Warning: Interpol_event outside bounds for time " << trajecto_time << endl;
    return acc_event;
}

AccuracyEvent XMls::Accuracy(XBlockIndex block_idx, XEchoIndex echo_idx)
{
    AccuracyEvent acc_event;
    if(!mp_trajecto->HasAcc())
    {
        cout << "Warning: Asked for an accuracy but no accuracy file provided, "
             << "returning uninitialized accuracy event" << endl;
        return acc_event;
    }
    double trajecto_time = Time(block_idx, echo_idx) + m_time_shift;
    if(!mp_trajecto->Accuracy().Interpol_event(acc_event, trajecto_time))
        cout << "Warning: Interpol_event outside bounds for time " << trajecto_time << endl;
    return acc_event;
}

