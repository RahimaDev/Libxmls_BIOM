
#define WHAT "SensorEdges: create a single .ply edge mesh for all the specified time interval using sensor topology"

#include <ctime>
#include <iostream>
#include <limits>
#include "libXMls/XMls.h"
#include "libXBase/XPt2D.h"


using namespace std;

struct param
{
    string sbet;
    string acc;
    string laser_calib;
    string ept_folder;
    string output, format;

    // optional
    float threshold, lambda;
    XSecond i_start, i_end;
    int pivot_E, pivot_N, pivot_H;
    int add_rgb;

    // header
    string datetime_;
    param():threshold(0.5), lambda(1), i_start(-1), i_end(-1), pivot_E(-1), pivot_N(-1), pivot_H(0), add_rgb(false){}
    string Format()
    {
        string format = output.substr(output.size()-3,3);
        cout << "format: " << format << endl;
        return format;
    }
};

struct Edge
{
    XBlockIndex b1, b2;
    XEchoIndex e1, e2;
    float C0, C1;
    Edge(XBlockIndex b1_=0, XBlockIndex b2_=0, XEchoIndex e1_=0, XEchoIndex e2_=0):
        b1(b1_), b2(b2_), e1(e1_), e2(e2_) {}
};

template <typename T> void Write(char * & it, T data)
{
    *reinterpret_cast<T*>(it) = data;
    it += sizeof(T);
}

void WritePly(XMls & mls, vector<Edge> v_edge, param params)
{
    ofstream fileOut(params.output);
    if(!fileOut.good())
    {
        cout << "Cannot open " + params.output + " for writing\n";
        return;
    }

    // write text header
    fileOut << "ply\nformat binary_little_endian 1.0" << endl;
    fileOut << "comment Generated from " << params.ept_folder << " secs " << mls.FirstSecond() << " to " << mls.LastSecond() << endl;
    fileOut << "comment IGN offset Pos " << params.pivot_E << " " << params.pivot_N << " " << params.pivot_H << endl;
    fileOut << "element vertex " << v_edge.size() << endl;
    fileOut << "property float x" << endl;
    fileOut << "property float y" << endl;
    fileOut << "property float z" << endl;
    fileOut << "property float x2" << endl;
    fileOut << "property float y2" << endl;
    fileOut << "property float z2" << endl;
    fileOut << "property float C0" << endl;
    fileOut << "property float C1" << endl;
    if(params.add_rgb == 1)
    {
        fileOut << "property uchar red" << endl;
        fileOut << "property uchar green" << endl;
        fileOut << "property uchar blue" << endl;
    } else if(params.add_rgb == 2)
    {
        fileOut << "property float quality" << endl;
    }
    fileOut << "end_header" << endl;
    // vertex list
    unsigned int edge_bytesize = 8*sizeof(float);
    unsigned long edge_buffer_size = edge_bytesize * v_edge.size();
    XFloatAttrib * p_refl=NULL;
    if(params.add_rgb > 0) p_refl = mls.GetEchoAttrib<XFloatAttrib>("reflectance");
    if(params.add_rgb == 1) edge_buffer_size += 3*sizeof(unsigned char) * v_edge.size();
    else if(params.add_rgb == 2) edge_buffer_size += sizeof(float) * v_edge.size();

    char * buffer = new char[edge_buffer_size], * it = buffer;
    XBlockIndex last_loaded_block(0);
    cout << "Loading block 0" << endl;
    mls.Load(0);
    XPt3D pivot(params.pivot_E, params.pivot_N, params.pivot_H);
    for(vector<Edge>::iterator edge_it=v_edge.begin(); edge_it<v_edge.end(); edge_it++)
    {
        // ensure b1 and b2 are loaded, and free before
        while(edge_it->b2 > last_loaded_block)
        {
            last_loaded_block++;
            cout << "Loading block" << last_loaded_block << endl;
            mls.Load(last_loaded_block);
            if(last_loaded_block>1)
            {
                cout << "Freeing block" << last_loaded_block-2 << endl;
                mls.Free(last_loaded_block-2);
            }
        }
        XPt3D P1 = mls.Pworld(edge_it->b1, edge_it->e1)-pivot, P2 = mls.Pworld(edge_it->b2, edge_it->e2)-pivot;
        Write<float>(it, P1.X);
        Write<float>(it, P1.Y);
        Write<float>(it, P1.Z);
        Write<float>(it, P2.X);
        Write<float>(it, P2.Y);
        Write<float>(it, P2.Z);
        Write<float>(it, edge_it->C0);
        Write<float>(it, edge_it->C1);
        if(params.add_rgb == 1) // r g b mode
        {
            float g = 12.75f*(0.5*(p_refl->at(edge_it->b1)[edge_it->e1]+p_refl->at(edge_it->b2)[edge_it->e2])+20.f); // rescale [-20,0] to [0,256] TODO: parameters
            if(g<0.f) g=0.f; else if(g>255.f) g=255.f;
            unsigned char ug = g;
            Write<unsigned char>(it, ug);
            Write<unsigned char>(it, ug);
            Write<unsigned char>(it, ug);
        }
        else if(params.add_rgb == 2) // quality mode
        {
            Write<float>(it, 0.5*(p_refl->at(edge_it->b1)[edge_it->e1]+p_refl->at(edge_it->b2)[edge_it->e2]));
        }
    }

    cout << "Writing " << v_edge.size() << " edges of size " << edge_buffer_size << "=" << 1.e-6*edge_buffer_size << "MB" << endl;
    fileOut.write(buffer, edge_buffer_size);
    delete buffer;
    fileOut.close();
}

inline XPt3D V12(XMls & mls, Edge & e, bool normalise=true)
{
    XPt3D V = mls.Pworld(e.b1, e.e1)-mls.Pworld(e.b2, e.e2);
    if(normalise) V.Normalise();
    return V;
}

//-----------------------------------------------------------------------------
int main(int argc, char **argv)
{
    cout << WHAT << endl;
    if(argc < 6)
    {
        cout << "Usage: " << argv[0] << "  sbet acc laser_calib.xml ept_folder output [threshold=0.5 lambda=1 (int)start_time (int)end_time pivot_E=-1 pivot_N=-1 pivot_H=0 add_rgb=0]" << endl;
        cout << "sbet: path to the sbet (trajecto) file" << endl;
        cout << "acc: path to the accuracy file (sometimes called smrmsg)" << endl;
        cout << "laser_calib.xml: calibration file for the laser" << endl;
        cout << "ept_folder: folder containing the echo pulse tables (generated with EptExport)" << endl;
        cout << "output: name of the output mesh file (.ply or .off only)" << endl;
        cout << "threshold: threshold on regularity to keep an edge (default=0.5)" << endl;
        cout << "lambda: weight on C1 regularity from 0 (non) to infinity (no C0 reg) (default=1)" << endl;
        //cout << "max_DP_error: maximum Douglas-Peucker error (default=0m=do not decimate)" << endl;
        cout << "(start|end)_time: start and end time of the laser points to export (default=everything)" << endl;
        cout << "pivot_(E|N|H): optionally set manually the pivot point (default=-1=first point of the trajectory rounded at 100m for EN, 0 for H)" << endl;
        cout << "add_rgb: if not 0, add r,g, b attributes in .ply mode based on reflectance (default=0=don't)" << endl;
        return 0;
    }

    int i_arg=1;
    param params;
    // required
    params.sbet = string((argv[i_arg++]));
    params.acc = string(argv[i_arg++]);
    params.laser_calib = string(argv[i_arg++]);
    params.ept_folder = string(argv[i_arg++]);
    params.output = string(argv[i_arg++]);

    // optional
    if(i_arg < argc) params.threshold = atof(argv[i_arg++]);
    if(i_arg < argc) params.lambda = atof(argv[i_arg++]);
    if(i_arg < argc) params.i_start = atoi(argv[i_arg++]);
    if(i_arg < argc) params.i_end = atoi(argv[i_arg++]);
    if(i_arg < argc) params.pivot_E = atoi(argv[i_arg++]);
    if(i_arg < argc) params.pivot_N = atoi(argv[i_arg++]);
    if(i_arg < argc) params.pivot_H = atoi(argv[i_arg++]);
    if(i_arg < argc) params.add_rgb = atoi(argv[i_arg++]);

    clock_t start = clock();
    // constructor and infos accessible after construction
    XTrajecto traj(params.sbet);
    XMls mls(params.ept_folder, params.laser_calib, &traj);
    mls.Select(params.i_start, params.i_end);
    cout << mls.NBlock() << " block(s) selected" << endl;

    if(params.pivot_E == -1)
    {
        XPt3D Pivot = traj.GetGeoref(0).Translation();
        params.pivot_E = 100*(int)(Pivot.X/100);
        params.pivot_N = 100*(int)(Pivot.Y/100);
    }

    start = clock();
    int PPL = mls.PulsePerLine();

    // create edge list
    vector<Edge> v_edge;
    for(XBlockIndex block_idx=0; block_idx<mls.NBlock(); block_idx++)
    {
        mls.Load(block_idx);
        cout << "Time " << block_idx << "=" << mls.Time(block_idx, XPulseIndex(0)) << endl;
        XPulseIndex n_pulse = mls.NPulse(block_idx);
        for(XPulseIndex pulse_idx=0; pulse_idx<n_pulse-3*PPL-3; pulse_idx++)
        {
            // regularity criteria involes 4 consecutive pulses along a line and applies to the middle edge
            // along sensor line
            XPt3D l1 = mls.RayWorld(block_idx, pulse_idx+1);
            XEchoIndex fe0 = mls.IdxFirstEcho(block_idx, pulse_idx); // first echo in pulse 0
            XEchoIndex fe1 = mls.IdxFirstEcho(block_idx, pulse_idx.Plus(1)); // first echo in pulse 1
            XEchoIndex fe2 = mls.IdxFirstEcho(block_idx, pulse_idx.Plus(2)); // first echo in pulse 2
            XEchoIndex fe3 = mls.IdxFirstEcho(block_idx, pulse_idx.Plus(3)); // first echo in pulse 3
            for(int e1=0; e1<mls.NbOfEcho(block_idx, pulse_idx.Plus(1)); e1++)
            {
                for(int e2=0; e2<mls.NbOfEcho(block_idx, pulse_idx.Plus(2)); e2++)
                {
                    Edge e12(block_idx, block_idx, fe1+e1, fe2+e2);
                    XPt3D v12 = V12(mls, e12);
                    float C0 = acos(fabs(prodScal(l1,v12))), C1 = 0;
                    bool pass=true;
                    if(C0 < params.threshold)
                    {
                        for(int e0=0; e0<mls.NbOfEcho(block_idx, pulse_idx); e0++)
                        {
                            Edge e01(block_idx, block_idx, fe0+e0, fe1+e1);
                            XPt3D v01 = V12(mls, e01);
                            float cur_C1 = fabs(prodScal(v01,v12));
                            if(cur_C1>C1) C1=cur_C1;
                        }
                        for(int e3=0; e3<mls.NbOfEcho(block_idx, pulse_idx.Plus(3)); e3++)
                        {
                            Edge e23(block_idx, block_idx, fe2+e2, fe3+e3);
                            XPt3D v23 = V12(mls, e23);
                            float cur_C1 = fabs(prodScal(v12,v23));
                            if(cur_C1>C1) C1=cur_C1;
                        }
                        C1 = acos(C1);
                        // remove edge if C0 is too low but above C1
                        if(C0 > params.lambda * C1) pass = false;
                    }
                    if(pass)
                    {
                        e12.C0 = C0;
                        e12.C1 = C1;
                        v_edge.push_back(e12);
                    }
                }
            }
            // TODO: 2 diagonals
        }
        mls.Free(block_idx);
    }
    cout << "Creating " << params.output << endl;
    //if(params.Format() == "off") WriteOff(mls, v_tri, idx, params);
    WritePly(mls, v_edge, params);
    cout << "Done in " << (double)(clock() - start) / (double) CLOCKS_PER_SEC << "s" << endl;
    return 0;
}

