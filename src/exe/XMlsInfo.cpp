
#define WHAT "XMlsInfo: displays some infos on an XMls (mobile laser scan)"

#include <ctime>
#include <iostream>
#include <limits>
#include "libXMls/XMls.h"


using namespace std;

template<typename T> class interval_t
{
public:
    T m_min, m_max;
    interval_t():m_min(numeric_limits<T>::max()),m_max(numeric_limits<T>::min()){}
    void Add(const T & t)
    {
        if(t<m_min) m_min=t;
        if(t>m_max) m_max=t;
    }
};

template<typename T> ostream& operator<<(ostream& out, const interval_t<T>& itv)
{
    if(std::numeric_limits<T>::is_integer) // convert chars to ints
        return out << "[" << (int)itv.m_min << "," << (int)itv.m_max << "]";
    return out << "[" << itv.m_min << "," << itv.m_max << "]";
}

//-----------------------------------------------------------------------------
int main(int argc, char **argv)
{
    cout << WHAT << endl;
    if(argc < 4)
    {
        cout << "Usage: " << argv[0] << "  sbet laser_calib.xml ept_folder [(int)start_time (int)end_time acc]" << endl;
        cout << "sbet: path to the sbet (trajecto) file" << endl;
        cout << "laser_calib.xml: calibration file for the laser" << endl;
        cout << "ept_folder: folder containing the echo pulse tables (generated with EptExport)" << endl;
        cout << "(start|end)_time: start and end time of the laser points to export (default=everything)" << endl;
        cout << "acc: path to the accuracy file (sometimes called smrmsg)" << endl;
        return 0;
    }

    int i_arg=1;
    // required
    string sbet(argv[i_arg++]);
    string laser_calib(argv[i_arg++]);
    string ept_folder(argv[i_arg++]);

    // optional
    XSecond i_start = -1, i_end = -1;
    if(i_arg < argc) i_start = atoi(argv[i_arg++]);
    if(i_arg < argc) i_end = atoi(argv[i_arg++]);
    string acc("");
    if(i_arg < argc) acc=string(argv[i_arg++]);


    clock_t start = clock();
    // constructor and infos accessible after construction
    XTrajecto traj(sbet, acc);
    XMls mls(ept_folder, laser_calib, &traj);
    cout << mls.m_time_shift << "s time shift between ept: " << mls.m_ept_time_pivot.ToString()
         << " and sbet: " << mls.m_trajecto_time_pivot.ToString() << endl;
    cout << laser_calib << "->" << mls.m_sensor_georef.InfoTexte() << endl;
    mls.Select(i_start, i_end);
    cout << mls.NBlock() << " block(s) selected" << endl;

    // 1 indexation attribute added to each, remove from count
    cout << mls.NPulseAttrib()-1 << "/" << mls.NEchoAttrib()-1 << " pulses/echos attributes found in " << ept_folder << endl;
    XAbsoluteTime laser_start = mls.m_ept_time_pivot.PlusSecs(mls.FirstSecond());
    XAbsoluteTime laser_end = mls.m_ept_time_pivot.PlusSecs(mls.LastSecond());
    cout << "Laser is in [" << mls.FirstSecond() << "," << mls.LastSecond() << "]=[" <<
            laser_start.ToString() << "," << laser_end.ToString() << "]=" << mls.LastSecond()-mls.FirstSecond() << "s" << endl;
    XAbsoluteTime trajecto_start = mls.m_trajecto_time_pivot.PlusSecs(traj.StartTime());
    XAbsoluteTime trajecto_end = mls.m_trajecto_time_pivot.PlusSecs(traj.EndTime());
    cout << "Trajecto is in [" << traj.StartTime() << "," << traj.EndTime() << "]=[" <<
            trajecto_start.ToString() << "," << trajecto_end.ToString() << "]=" <<
            traj.EndTime() - traj.StartTime() << "s" << endl;

    // Time info and stats on loaded laser
    laser_start = mls.m_ept_time_pivot.PlusSecs(i_start);
    laser_end = mls.m_ept_time_pivot.PlusSecs(i_end);
    float time_span = i_end-i_start, micro=1e-6;
    cout.precision(4);
    cout << "--Laser--\n" << micro*mls.NTotalEcho() << " Mechos " << micro*mls.NTotalPulse() << " Mpulses in [" <<
            i_start << ", " << i_end << "]=[" << laser_start.ToString() << ", " << laser_end.ToString()
         << "]=" << time_span << "s selected in " << (double)(clock() - start) / (double) CLOCKS_PER_SEC << "s" << endl;
    float n_line = (float)mls.NTotalPulse()/(float)mls.PulsePerLine();
    cout << mls.PulsePerLine() << " pulses per line, " << n_line << " lines " << n_line/time_span << " lines per sec" << endl;

    // Time info and stats on loaded trajecto
    trajecto_start = mls.m_trajecto_time_pivot.PlusSecs(traj.StartTime());
    trajecto_end = mls.m_trajecto_time_pivot.PlusSecs(traj.EndTime());
    cout << "--Trajecto--\nLoaded " << traj.Nevent() << " sbet events and " <<
            traj.Accuracy().Nevent() << " acc events in [" <<
            traj.StartTime() << "," << traj.EndTime() << "]=["
         << trajecto_start.ToString() << ", " << trajecto_end.ToString() << "]" << endl;

    start = clock();
    // sample iteration on pulses
    interval_t<double> time_itv;
    interval_t<float> theta_itv, phi_itv;
    interval_t<int> nbOfEcho_itv;
    double d_sum=0., d_min=1.e9, d_max=0.;
    int n_int_line = 0;
    for(XBlockIndex block_idx=0; block_idx<mls.NBlock(); block_idx++)
    {
        mls.Load(block_idx);
        cout << "Block " << block_idx << ", n_pulse=" << mls.NPulse(block_idx) << endl;
        for(XPulseIndex pulse_idx=0; pulse_idx<mls.NPulse(block_idx); pulse_idx++)
        {
            time_itv.Add(mls.Time(block_idx, pulse_idx));
            theta_itv.Add(mls.Theta(block_idx, pulse_idx));
            phi_itv.Add(mls.Phi(block_idx, pulse_idx));
            nbOfEcho_itv.Add(mls.NbOfEcho(block_idx, pulse_idx));
            if(pulse_idx%mls.PulsePerLine() == 0)
            {
                XPulseIndex next_line_idx = pulse_idx + mls.PulsePerLine();
                if(next_line_idx < mls.NPulse(block_idx))
                {
                    XPt3D O = mls.Ins(block_idx, pulse_idx).Translation();
                    XPt3D O_next = mls.Ins(block_idx, next_line_idx).Translation();
                    double d = (O_next-O).Norme();
                    d_sum += d;
                    if(d<d_min) d_min = d;
                    if(d>d_max) d_max = d;
                    n_int_line++;
                }
            }
        }
        mls.Free(block_idx);
    }
    cout << "Trajecto length: " << d_sum << "m distance between lines min " << 100.*d_min << "cm avg " << 100.*d_sum/n_int_line << "cm max " << 100.*d_max << "cm" << endl;
    float dt = time_span/n_line;
    cout << "Speed: min " <<
            d_min/dt << "m/s (" << 3.6*d_min/dt << "km/h) avg " <<
            d_sum/time_span << "m/s (" << 3.6*d_sum/time_span << "km/h) max " <<
            d_max/dt << "m/s (" << 3.6*d_max/dt << "km/h)" << endl;
    cout << "Time in " << time_itv << endl;
    cout << "Theta in " << theta_itv << endl;
    cout << "Phi in " << phi_itv << endl;
    cout << "nbOfEcho in " << nbOfEcho_itv << endl;

    // sample iteration on echos
    interval_t<float> range_itv, amplitude_itv, reflectance_itv;
    interval_t<float> x_itv, y_itv, z_itv;
    interval_t<float> xw_itv, yw_itv, zw_itv;
    interval_t<int> numEcho_itv;
    interval_t<unsigned char> deviation_itv;
    // get required optional attributes
    XFloatAttrib * p_ampl = mls.GetEchoAttrib<XFloatAttrib>("amplitude");
    XFloatAttrib * p_refl = mls.GetEchoAttrib<XFloatAttrib>("reflectance");
    XUCharAttrib * p_dev = mls.GetEchoAttrib<XUCharAttrib>("deviation");
    for(XBlockIndex block_idx=0; block_idx<mls.NBlock(); block_idx++)
    {
        mls.Load(block_idx);
        cout << "Block " << block_idx << ", n_echo=" << mls.NEcho(block_idx) << endl;
        for(XEchoIndex echo_idx=0; echo_idx<mls.NEcho(block_idx); echo_idx++)
        {
            range_itv.Add(mls.Range(block_idx, echo_idx));
            amplitude_itv.Add(p_ampl->at(block_idx)[echo_idx]);
            reflectance_itv.Add(p_refl->at(block_idx)[echo_idx]);
            deviation_itv.Add(p_dev->at(block_idx)[echo_idx]);
            numEcho_itv.Add(mls.NumEcho(block_idx, echo_idx));
            XPt3D P = mls.P(block_idx, echo_idx), Pw=mls.Pworld(block_idx, echo_idx);
            x_itv.Add(P.X);
            y_itv.Add(P.Y);
            z_itv.Add(P.Z);
            xw_itv.Add(Pw.X);
            yw_itv.Add(Pw.Y);
            zw_itv.Add(Pw.Z);
        }
        mls.Free(block_idx);
    }
    cout << "Range in " << range_itv << endl;
    cout << "Amplitude in " << amplitude_itv << endl;
    cout << "Reflectance in " << reflectance_itv << endl;
    cout << "Deviation in " << deviation_itv << endl;
    cout << "numEcho in " << numEcho_itv << endl;
    cout << "X in " << x_itv << endl;
    cout << "Y in " << y_itv << endl;
    cout << "Z in " << z_itv << endl;
    cout << "Xw in " << xw_itv << endl;
    cout << "Yw in " << yw_itv << endl;
    cout << "Zw in " << zw_itv << endl;
    cout << "Done in " << (double)(clock() - start) / (double) CLOCKS_PER_SEC << "s" << endl;
    return 0;
}

