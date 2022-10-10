
#define WHAT "XMlsAttrib: example of computing and adding a new attribute to a XMls."

#include "libXMls/XMls.h"

using namespace std;

//-----------------------------------------------------------------------------
int main(int argc, char **argv)
{
    cout << WHAT << endl;
    if(argc < 5)
    {
        cout << "Usage: " << argv[0] << "  sbet acc laser_calib.xml ept_folder [(int)start_time (int)end_time]" << endl;
        cout << "sbet: path to the sbet (trajecto) file" << endl;
        cout << "acc: path to the accuracy file (sometimes called smrmsg)" << endl;
        cout << "laser_calib.xml: calibration file for the laser" << endl;
        cout << "ept_folder: folder containing the echo pulse tables (generated with EptExport)" << endl;
        return 0;
    }

    int i_arg=1;
    // required
    string sbet(argv[i_arg++]);
    string acc(argv[i_arg++]);
    string laser_calib(argv[i_arg++]);
    string ept_folder(argv[i_arg++]);

    // optional
    XSecond i_start = -1, i_end = -1;
    if(i_arg < argc) i_start = atoi(argv[i_arg++]);
    if(i_arg < argc) i_end = atoi(argv[i_arg++]);

    // read all echo/pulse tables in time interval
    clock_t start = clock();
    XTrajecto traj(sbet, acc);
    XMls mls(ept_folder, laser_calib, &traj);
    cout << mls.NPulseAttrib() << "/" << mls.NEchoAttrib() << " echos/pulses attributes found in " << ept_folder << endl;
    mls.Select(i_start, i_end);
    cout << "Selected in " << (double)(clock() - start) / (double) CLOCKS_PER_SEC << "s" << endl;
    XBlockIndex last_block_idx=mls.NBlock()-1;
    cout << "Time " << mls.Time(0, (XPulseIndex)0) << "-" << mls.Time(last_block_idx, (XPulseIndex)(mls.NPulse(last_block_idx)-1)) << endl;
    cout << "Pulses per line: " << mls.PulsePerLine() << endl;

    start = clock();
    // iteration on echos
    XFloatAttrib * xw = mls.AddEchoAttrib<XFloatAttrib>("xw");
    XFloatAttrib * yw = dynamic_cast<XFloatAttrib *>(mls.AddEchoAttrib("yw", "float32"));
    XFloatAttrib * zw = dynamic_cast<XFloatAttrib *>(mls.AddEchoAttrib("zw", "float32"));
    for(XBlockIndex block_idx=0; block_idx<mls.NBlock(); block_idx++)
    {
        mls.Load(block_idx);
        cout << "Block " << block_idx << ", n_echo=" << mls.NEcho(block_idx) << endl;
        for(XEchoIndex echo_idx=0; echo_idx<mls.NEcho(block_idx); echo_idx++)
        {
            XPt3D Pw=mls.Pworld(block_idx, echo_idx);
            xw->at(block_idx).at(echo_idx) = Pw.X;
            yw->at(block_idx).at(echo_idx) = Pw.Y;
            zw->at(block_idx).at(echo_idx) = Pw.Z;
        }
        xw->Save(block_idx);
        yw->Save(block_idx);
        zw->Save(block_idx);
        mls.Free(block_idx);
    }
    cout << "Done in " << (double)(clock() - start) / (double) CLOCKS_PER_SEC << "s" << endl;
    return 0;
}

