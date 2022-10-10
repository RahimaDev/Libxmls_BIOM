#include "libXMls/XTrajecto.h"

#ifndef _USE_MATH_DEFINES
#define  _USE_MATH_DEFINES
#endif 
#include <math.h>

using namespace std;

// should be continuous and in [0, 2pi] for all alpha in [0,1], assuming a0 and a1 are in [0, 2pi]
inline double InterpolAngle(double a0, double a1, double alpha)
{
    double da = a1-a0;
    while(da < -M_PI) da += 2*M_PI;
    while(da > M_PI) da -= 2*M_PI;
    double ret = a0+alpha*da;
    while(ret < 0) ret += 2*M_PI;
    while(ret > 2*M_PI) ret -= 2*M_PI;
    return ret;
}

// interpolation factor for a time between 2 sampling times
inline double Alpha(const double & t0, const double & t1, double time)
{
    double dt = t1 - t0;
    return (dt<1.e-10?0.5:(time-t0)/dt);
}

// linear interpolation in time
SbetEvent InterpolPosOnly(const SbetEvent & s0, SbetEvent & s1, double time)
{
    double alpha=Alpha(s0.m_time, s1.m_time, time), uma = 1.-alpha;
    SbetEvent s;
    s.m_time = time;
    s.m_latitude = uma*s0.m_latitude + alpha*s1.m_latitude;
    s.m_longitude = uma*s0.m_longitude + alpha*s1.m_longitude;
    s.m_altitude = uma*s0.m_altitude + alpha*s1.m_altitude;
    s.m_roll = InterpolAngle(s0.m_roll, s1.m_roll, alpha);
    s.m_pitch = InterpolAngle(s0.m_pitch, s1.m_pitch, alpha);
    s.m_plateformHeading = InterpolAngle(s0.m_plateformHeading, s1.m_plateformHeading, alpha);
    s.m_wanderAngle = InterpolAngle(s0.m_wanderAngle, s1.m_wanderAngle, alpha);
    return s;
}

// linear interpolation in time
SbetEvent Interpol(const SbetEvent & s0, SbetEvent & s1, double time)
{
    double alpha=Alpha(s0.m_time, s1.m_time, time), uma = 1.-alpha;
    SbetEvent s;
    s.m_time = time;
    s.m_xVelocity = uma*s0.m_xVelocity + alpha*s1.m_xVelocity;
    s.m_yVelocity = uma*s0.m_yVelocity + alpha*s1.m_yVelocity;
    s.m_zVelocity = uma*s0.m_zVelocity + alpha*s1.m_zVelocity;
    s.m_xAcceleration = uma*s0.m_xAcceleration + alpha*s1.m_xAcceleration;
    s.m_yAcceleration = uma*s0.m_yAcceleration + alpha*s1.m_yAcceleration;
    s.m_zAcceleration = uma*s0.m_zAcceleration + alpha*s1.m_zAcceleration;
    s.m_xBodyAngularRate = uma*s0.m_xBodyAngularRate + alpha*s1.m_xBodyAngularRate;
    s.m_yBodyAngularRate = uma*s0.m_yBodyAngularRate + alpha*s1.m_yBodyAngularRate;
    s.m_zBodyAngularRate = uma*s0.m_zBodyAngularRate + alpha*s1.m_zBodyAngularRate;
    s.m_latitude = uma*s0.m_latitude + alpha*s1.m_latitude;
    s.m_longitude = uma*s0.m_longitude + alpha*s1.m_longitude;
    s.m_altitude = uma*s0.m_altitude + alpha*s1.m_altitude;
    s.m_roll = InterpolAngle(s0.m_roll, s1.m_roll, alpha);
    s.m_pitch = InterpolAngle(s0.m_pitch, s1.m_pitch, alpha);
    s.m_plateformHeading = InterpolAngle(s0.m_plateformHeading, s1.m_plateformHeading, alpha);
    s.m_wanderAngle = InterpolAngle(s0.m_wanderAngle, s1.m_wanderAngle, alpha);
    return s;
}

ostream& operator<<(ostream& out, const SbetEvent& s)
{
    return out << "sol"                 <<"\t"
               << s.m_time              <<"\t"
               << s.m_latitude          <<"\t"
               << s.m_longitude         <<"\t"
               << s.m_altitude          <<"\t"
               << s.m_xVelocity         <<"\t"
               << s.m_yVelocity         <<"\t"
               << s.m_zVelocity         <<"\t"
               << s.m_roll              <<"\t"
               << s.m_pitch             <<"\t"
               << s.m_plateformHeading  <<"\t"
               << s.m_wanderAngle       <<"\t"
               << s.m_xAcceleration     <<"\t"
               << s.m_yAcceleration     <<"\t"
               << s.m_zAcceleration     <<"\t"
               << s.m_xBodyAngularRate  <<"\t"
               << s.m_yBodyAngularRate  <<"\t"
               << s.m_zBodyAngularRate;
}

// linear interpolation in time
AccuracyEvent Interpol(const AccuracyEvent & a0, AccuracyEvent & a1, double time)
{
    double alpha = Alpha(a0.m_time, a1.m_time, time), uma = 1.-alpha;
    AccuracyEvent a;
    a.m_time = time;
    a.m_northPositionRMSError = uma*a0.m_northPositionRMSError + alpha*a1.m_northPositionRMSError;
    a.m_eastPositionRMSError = uma*a0.m_eastPositionRMSError + alpha*a1.m_eastPositionRMSError;
    a.m_downPositionRMSError = uma*a0.m_downPositionRMSError + alpha*a1.m_downPositionRMSError;

    a.m_northVelocityRMSError = uma*a0.m_northVelocityRMSError + alpha*a1.m_northVelocityRMSError;
    a.m_eastPositionRMSError = uma*a0.m_eastPositionRMSError + alpha*a1.m_eastPositionRMSError;
    a.m_downPositionRMSError = uma*a0.m_downPositionRMSError + alpha*a1.m_downPositionRMSError;

    a.m_RollRMSError = uma*a0.m_RollRMSError + alpha*a1.m_RollRMSError;
    a.m_PitchRMSError = uma*a0.m_PitchRMSError + alpha*a1.m_PitchRMSError;
    a.m_headingRMSError = uma*a0.m_headingRMSError + alpha*a1.m_headingRMSError;
    return a;
}

ostream& operator<<(ostream& out, const AccuracyEvent& a)
{
    return out << "acc"                    <<"\t"
               << a.m_time                 <<"\t"
               << a.m_northPositionRMSError<<"\t"
               << a.m_eastPositionRMSError <<"\t"
               << a.m_downPositionRMSError <<"\t"
               << a.m_northVelocityRMSError<<"\t"
               << a.m_eastVelocityRMSError <<"\t"
               << a.m_downVelocityRMSError <<"\t"
               << a.m_RollRMSError         <<"\t"
               << a.m_PitchRMSError        <<"\t"
               << a.m_headingRMSError;
}

XGeorefEventSeries::XGeorefEventSeries(string filename, string strProj4_out, string strProj4_in):
    TEventSeries<SbetEvent>(filename)
{
    SetInputProjection(strProj4_in);
    SetOutputProjection(strProj4_out);
}

void XGeorefEventSeries::SetInputProjection(std::string strProj4_in)
{
    if(!(m_proj_in = pj_init_plus(strProj4_in.c_str())))
    {
        cout<<"ERROR: wrong proj_in initialization " << strProj4_in << endl;
        cout<<"maybe the PROJ_LIB variable is not set to the NAD dir (something like /home/bvallet/dev/extern/proj4/nad)" << endl;
    }
}

void XGeorefEventSeries::SetOutputProjection(std::string strProj4_out)
{
    if(!(m_proj_out = pj_init_plus(strProj4_out.c_str())))
    {
        cout<<"ERROR: wrong proj_out initialization " << strProj4_out << endl;
        cout<<"maybe the PROJ_LIB variable is not set to the NAD dir (something like /home/bvallet/dev/extern/proj4/nad)" << endl;
    }
}

bool XGeorefEventSeries::PreComputeGeoref()
{
    if(mv_event.empty())
    {
        cout << "ERROR: Called PreComputeGeoref() with empty trajectory" << endl;
        return false;
    }
    mv_georef.resize(m_n_event);
    if(true) // copy coords in a tab to accelerate proj transforms
    {
        std::vector<double> x(m_n_event), y(m_n_event), z(m_n_event);
        for(uint i=0; i<m_n_event; i++)
        {
            x[i]=mv_event[i].m_longitude;
            y[i]=mv_event[i].m_latitude;
            z[i]=mv_event[i].m_altitude;
        }
        int tmp = pj_transform(m_proj_in, m_proj_out, m_n_event, 1, &(x[0]), &(y[0]), &(z[0]));
        if(tmp) cout << "Warning: proj error " << tmp << endl;
        for(uint i=0; i<m_n_event; i++)
        {
            mv_georef[i].Translation(XPt3D(x[i],y[i],z[i]));
            mv_georef[i].Rotation(GetRotation(mv_event[i]));
        }
    }
    else if(false) // version with only 1 copy, does not work with Miloud
    {
        //cout << "Copying..." << endl;
        for(uint i=0; i<m_n_event; i++)
            mv_georef[i].Translation(XPt3D(mv_event[i].m_longitude,
                                           mv_event[i].m_latitude,
                                           mv_event[i].m_altitude));

        //cout << "Transform..." << endl;
        int tmp = pj_transform(m_proj_in, m_proj_out, m_n_event,
                               sizeof(XArchiGeoref)/sizeof(double),
                               &(mv_georef[0].Translation().X),
                &(mv_georef[0].Translation().Y),
                &(mv_georef[0].Translation().Z));
        if(tmp) cout << "Warning: proj error " << tmp << endl;

        //cout << "Rotation..." << endl;
        for(uint i=0; i<m_n_event; i++) mv_georef[i].Rotation(GetRotation(mv_event[i]));
    }
    else for(uint i=0; i<m_n_event; i++) GetGeoref(mv_event[i], mv_georef[i]); // secure but slow
    double conv = GetConvMeridien(mv_event[0]);
    cout<<"convergence of merdien: "<<conv<<" rad = "<<conv/M_PI*180<<" dd"<<endl;
    return true;
}

void XGeorefEventSeries::GetGeoref(const SbetEvent & event, XArchiGeoref & georef)
{
    georef.Translation( GetTranslation(event) );
    georef.Rotation( GetRotation(event) );
}

/// Get the translation only
XPt3D XGeorefEventSeries::GetTranslation(const SbetEvent & event)
{
    double x=event.m_longitude, y=event.m_latitude, z=event.m_altitude;
    int tmp = pj_transform(m_proj_in, m_proj_out, 1, 1, &x, &y, &z);
    if(tmp) cout << "Warning: proj error " << tmp << endl;
    return XPt3D(x, y, z);
}
bool XGeorefEventSeries::GetTranslation(double time, XPt3D & P, bool precomputed)
{
    if(precomputed)
    {
        if(!m_meta_info_present) return false;
        if(mv_georef.size() != m_n_event) PreComputeGeoref();
        int index = PrevIndex(time);
        if(index < 0) {P = mv_georef.front().Translation(); return false;}
        if(index >= (int)m_n_event-1) {P = mv_georef.back().Translation(); return false;}
        double t0 = mv_event[index].m_time, t1 = mv_event[index+1].m_time;
        P = XPt3D(mv_georef[index].Translation(), mv_georef[index+1].Translation(), (time-t0)/(t1-t0));
        return true;
    }
    SbetEvent event;
    bool OK = Interpol_event(event, time);
    P = GetTranslation(event);
    return OK;
}

double XGeorefEventSeries::GetConvMeridien(SbetEvent event)
{
    return -0.72537437089 * (event.m_longitude-0.0523598775598); // sin(0.811578102)=0.72537437089
}

XMat3D XGeorefEventSeries::GetRotation(const SbetEvent & event)
{
    double conv_meridien = -GetConvMeridien(event);
    double cp = cos(event.m_pitch), sp = sin(event.m_pitch);
    double cr = cos(event.m_roll), sr = sin(event.m_roll);
    if(false) // reference impl
    {
        double cc = cos(conv_meridien), sc = sin(conv_meridien);
        double heading = event.m_plateformHeading-event.m_wanderAngle;
        double ch = cos(heading), sh = sin(heading);
        XMat3D M_ROLL=XMat3D(     XPt3D( 1,  0,   0 ),
                                  XPt3D( 0, cr, -sr ),
                                  XPt3D( 0, sr,  cr ) );

        XMat3D M_PITCH=XMat3D(    XPt3D(  cp, 0, sp ),
                                  XPt3D(   0, 1,  0 ),
                                  XPt3D( -sp, 0, cp ) );

        XMat3D M_HEADING=XMat3D(  XPt3D( ch,-sh, 0),
                                  XPt3D( sh, ch, 0),
                                  XPt3D(  0,  0, 1) );

        XMat3D M_NEB_ENH= XMat3D( XPt3D( 0, 1, 0),
                                  XPt3D( 1, 0, 0),
                                  XPt3D( 0, 0,-1) );

        /*XMat3D M_NWU_ENU = XMat3D( XPt3D( 0,-1, 0),
                                   XPt3D( 1, 0, 0),
                                   XPt3D( 0, 0, 1) );*/

        XMat3D M_CONV= XMat3D(  XPt3D( cc,-sc, 0),
                                XPt3D( sc, cc, 0),
                                XPt3D(  0,  0, 1) );

        return M_CONV * M_NEB_ENH * M_HEADING * M_PITCH * M_ROLL;
    }
    else // optimized impl
    {
        double heading = event.m_plateformHeading-event.m_wanderAngle-conv_meridien;
        double ch = cos(heading), sh = sin(heading);
        return XMat3D(cp*sh, sr*sp*sh+cr*ch, cr*sp*sh-sr*ch,
                      cp*ch, sr*sp*ch-cr*sh, cr*sp*ch+sr*sh,
                      sp   ,-sr*cp         ,-cr*cp         );
    }
}

bool XGeorefEventSeries::GetGeoref(double time, XArchiGeoref & georef)
{
    SbetEvent event;
    bool OK = Interpol_event(event, time);
    GetGeoref(event, georef);
    return OK;
}

bool XGeorefEventSeries::GetGeoref_precomputed(double time, XArchiGeoref & georef)
{
    if(!m_meta_info_present) return false;
    if(mv_georef.size() != m_n_event) PreComputeGeoref();
    int index = PrevIndex(time);
    if(index < 0) {georef = mv_georef.front(); return false;}
    if(index >= (int)m_n_event-1) {georef = mv_georef.back(); return false;}
    double t0 = mv_event[index].m_time, t1 = mv_event[index+1].m_time;
    georef = XArchiGeoref(mv_georef[index], mv_georef[index+1], (time-t0)/(t1-t0));
    return true;
}

XPt3D XGeorefEventSeries::ApplyGeoref_precomputed(double time, const XPt3D & P)
{
    if(!m_meta_info_present) return false;
    if(mv_georef.size() != m_n_event) PreComputeGeoref();
    int index = PrevIndex(time);
    if(index < 0) index=0;
    if(index >= (int)m_n_event-1) index = m_n_event-2;
    double t0 = mv_event[index].m_time, t1 = mv_event[index+1].m_time;
    double alpha = (time-t0)/(t1-t0);
    return (1-alpha)*mv_georef[index].Applique_transfo(P)
            +alpha*mv_georef[index+1].Applique_transfo(P);
}

/// apply a translation to the georefevent
void XGeorefEventSeries::ApplyTranslation(unsigned int i, const XPt3D & P)
{
    // new_translation is in m_proj_out
    XPt3D new_translation = mv_georef.at(i).Translation()+P;
    mv_georef.at(i).Translation() = new_translation;
    int tmp = pj_transform(m_proj_out, m_proj_in, 1, 1,
                           &new_translation.X, &new_translation.Y, &new_translation.Z);
    // new_translation is now in m_proj_in
    if(tmp) cout << "Warning: proj error " << tmp << endl;
    mv_event.at(i).m_latitude = new_translation.X;
    mv_event.at(i).m_longitude = new_translation.Y;
    mv_event.at(i).m_altitude = new_translation.Z;
}

template <typename T> void Write(char * & it, T data)
{
    *reinterpret_cast<T*>(it) = data;
    it += sizeof(T);
}

void XGeorefEventSeries::ExportPly(std::string ply_filename)
{
    vector<bool> keep;
    ExportPly(ply_filename, keep);
}

void XGeorefEventSeries::ExportPly(std::string ply_filename, vector<bool> & keep)
{
    int n_event = 0;
    if(keep.size() == m_n_event) {for(uint i=0; i<m_n_event; i++) if(keep.at(i)) n_event++;}
    else n_event = m_n_event;
    ofstream fileOut(ply_filename.c_str());
    if(!fileOut.good())
    {
        cout << "Cannot open " + ply_filename + " for writing\n";
        return;
    }
    // write text header
    fileOut << "ply\nformat binary_little_endian 1.0" << endl;
    fileOut << "comment Generated from " << m_filename << endl;
    XPt3D pivot = mv_georef[0].Translation();
    int pivot_E = 100*(int)(pivot.X/100);
    int pivot_N = 100*(int)(pivot.Y/100);
    fileOut << "comment IGN offset Pos " << pivot_E << " " << pivot_N << " 0" << endl;
    fileOut << "element vertex " << n_event << endl;
    fileOut << "property float GPS_time" << endl;
    fileOut << "property float latitude" << endl;
    fileOut << "property float longitude" << endl;
    fileOut << "property float altitude" << endl;
    fileOut << "property float xVelocity" << endl;
    fileOut << "property float yVelocity" << endl;
    fileOut << "property float zVelocity" << endl;
    fileOut << "property float roll" << endl;
    fileOut << "property float pitch" << endl;
    fileOut << "property float plateformHeading" << endl;
    fileOut << "property float wanderAngle" << endl;
    fileOut << "property float xAcceleration" << endl;
    fileOut << "property float yAcceleration" << endl;
    fileOut << "property float zAcceleration" << endl;
    fileOut << "property float xBodyAngularRate" << endl;
    fileOut << "property float yBodyAngularRate" << endl;
    fileOut << "property float zBodyAngularRate" << endl;
    fileOut << "property float x" << endl;
    fileOut << "property float y" << endl;
    fileOut << "property float z" << endl;
    fileOut << "end_header" << endl;
    unsigned int echo_bytesize = 20*sizeof(float);
    unsigned long buffer_size = echo_bytesize * n_event;
    char * buffer = new char[buffer_size], * it = buffer;
    // compute trajecto length just for display
    XPt3D P_old = mv_georef[0].Translation();
    double trajecto_length = 0;
    for(uint i=0; i<m_n_event; i++) if(keep.size() != m_n_event || keep.at(i))
    {
        SbetEvent evt = mv_event[i];
        Write<float>(it, evt.m_time);//sec
        Write<float>(it, evt.m_latitude);//radians
        Write<float>(it, evt.m_longitude);//radians
        Write<float>(it, evt.m_altitude);//meters
        Write<float>(it, evt.m_xVelocity);//meters/second
        Write<float>(it, evt.m_yVelocity);//meters/second
        Write<float>(it, evt.m_zVelocity);//meters/second
        Write<float>(it, evt.m_roll);//radians
        Write<float>(it, evt.m_pitch);//radians
        Write<float>(it, evt.m_plateformHeading);//radians
        Write<float>(it, evt.m_wanderAngle);//radians
        Write<float>(it, evt.m_xAcceleration);//meters/second2
        Write<float>(it, evt.m_yAcceleration);//meters/second2
        Write<float>(it, evt.m_zAcceleration);//meters/second2
        Write<float>(it, evt.m_xBodyAngularRate);//radians/second
        Write<float>(it, evt.m_yBodyAngularRate);//radians/second
        Write<float>(it, evt.m_zBodyAngularRate);//radians/second
        XPt3D P = mv_georef[i].Translation();
        trajecto_length += (P-P_old).Norme();
        P_old=P;
        Write<float>(it, P.X-pivot_E);
        Write<float>(it, P.Y-pivot_N);
        Write<float>(it, P.Z);
    }
    cout << "Writing " << n_event << " sbet events of size " << echo_bytesize << "=" << 1.e-6*buffer_size << "MB to " << ply_filename << endl;
    cout << "Trajecto length (in m): " << trajecto_length << endl;
    cout << "Time span (in s): " << mv_event[m_n_event-1].m_time - mv_event[0].m_time << endl;
    fileOut.write(buffer, buffer_size); // todo: split buffer if too big
    fileOut.close();
    delete buffer;
}

//-----------------------------------------------------------------------------
XTrajecto::XTrajecto(string sbetFile, string accFile, string strProj4_out, string strProj4_in):
    XGeorefEventSeries(sbetFile, strProj4_out, strProj4_in), m_accuracy(accFile)
{}

//-----------------------------------------------------------------------------
void XTrajecto::SetParams(string sbetFile, string accFile)
{
    XTrajecto::SetFilename(sbetFile);
    m_accuracy.SetFilename(accFile);
}
bool XTrajecto::GetMetaInfo()
{
    bool success=true;
    if(HasAcc()) success=m_accuracy.GetMetaInfo();
    return success && XGeorefEventSeries::GetMetaInfo();
}

bool XTrajecto::Load(double t_min, double t_max)
{
    bool success=true;
    if(HasAcc()) success=m_accuracy.Load(t_min, t_max);
    return success && XGeorefEventSeries::Load(t_min, t_max);
}

bool XTrajecto::Unload()
{
    bool success=true;
    if(HasAcc()) success=m_accuracy.Unload();
    return success && XGeorefEventSeries::Unload();
}

bool XTrajecto::Save(std::string sbet_filename, std::string acc_filename)
{
    bool success=true;
    if(HasAcc() && acc_filename.size()>0) success=m_accuracy.Save(acc_filename);
    return success && XGeorefEventSeries::Save(sbet_filename);
}

bool XTrajecto::Interpol(SbetEvent & sol_event, AccuracyEvent & acc_event, double time)
{
    bool success=true;
    if(HasAcc()) success=m_accuracy.Interpol_event(acc_event, time);
    return success && XGeorefEventSeries::Interpol_event(sol_event, time);
}
