
#define WHAT "SensorMesh3: create a single .ply mesh for all the specified time interval using sensor topology and filtering out accumulated lines"

#include <ctime>
#include <iostream>
#include <limits>
#include "libXMls/XMls.h"
#include "libXBase/XPt2D.h"


using namespace std;

struct param
{
    XMls::Param xmls_param;
    string output, format;
    double block_length;

    // optional
    float dist_threshold;
    float tri_threshold;
    float max_DP_error;
    int pivot_E, pivot_N, pivot_H;
    int add_rgb;

    // header
    string datetime_;
    param():block_length(50.), dist_threshold(0.02), tri_threshold(0.5), max_DP_error(0), i_start(-1), i_end(-1), pivot_E(-1), pivot_N(-1), pivot_H(0), add_rgb(false){}
    string Format()
    {
        string format = output.substr(output.size()-3,3);
        cout << "format: " << format << endl;
        return format;
    }
};

struct Triangle
{
    unsigned int i,j,k; // pulse with return idx
    Triangle(unsigned int i_=0, unsigned int j_=0, unsigned int k_=0):
        i(i_), j(j_), k(k_) {}
};

template <typename T> void Write(char * & it, T data)
{
    *reinterpret_cast<T*>(it) = data;
    it += sizeof(T);
}

void WritePly(XMls & mls, vector<Triangle> v_tri,
              int n_vertex, param params)
{
    ofstream fileOut(params.output);
    if(!fileOut.good())
    {
        cout << "Cannot open " + params.output + " for writing\n";
        return;
    }

    // write text header
    fileOut << "ply\nformat binary_little_endian 1.0" << endl;
    fileOut << "comment Generated from " << params.xmls_param.ept_folder << " secs " << mls.FirstSecond() << " to " << mls.LastSecond() << endl;
    fileOut << "comment IGN offset Pos " << params.pivot_E << " " << params.pivot_N << " " << params.pivot_H << endl;
    fileOut << "element vertex " << n_vertex << endl;
    fileOut << "property float x" << endl;
    fileOut << "property float y" << endl;
    fileOut << "property float z" << endl;
    if(params.add_rgb == 1)
    {
        fileOut << "property uchar red" << endl;
        fileOut << "property uchar green" << endl;
        fileOut << "property uchar blue" << endl;
    } else if(params.add_rgb == 2)
    {
        fileOut << "property float quality" << endl;
    }
    fileOut << "element face " << v_tri.size() << endl;
    fileOut << "property list uchar int vertex_indices" << endl;
    fileOut << "end_header" << endl;
    // vertex list
    unsigned int vertex_bytesize = 3*sizeof(float);
    unsigned long vertex_buffer_size = vertex_bytesize * n_vertex;
    XFloatAttrib * p_refl=NULL;
    if(params.add_rgb > 0) p_refl = mls.GetEchoAttrib<XFloatAttrib>("reflectance");
    if(params.add_rgb == 1) vertex_buffer_size += 3*sizeof(unsigned char) * n_vertex;
    else if(params.add_rgb == 2) vertex_buffer_size += sizeof(float) * n_vertex;

    char * buffer = new char[vertex_buffer_size], * it = buffer;
    for(XBlockIndex block_idx=0; block_idx<mls.NBlock(); block_idx++)
    {
        mls.Load(block_idx);
        for(XPulseIndex pulse_idx=0; pulse_idx<mls.NPulse(block_idx); pulse_idx++) if(mls.NbOfEcho(block_idx, pulse_idx)>0)
        {
            XEchoIndex last_echo_idx = mls.IdxLastEcho(block_idx, pulse_idx);
            XPt3D Pw = mls.Pworld(block_idx, last_echo_idx);
            float e=Pw.X-params.pivot_E, n=Pw.Y-params.pivot_N, h=Pw.Z-params.pivot_H;
            Write<float>(it, e);
            Write<float>(it, n);
            Write<float>(it, h);
            if(params.add_rgb == 1) // r g b mode
            {
                float g = 12.75f*(p_refl->at(block_idx)[last_echo_idx]+20.f); // rescale [-20,0] to [0,256] TODO: parameters
                if(g<0.f) g=0.f; else if(g>255.f) g=255.f;
                unsigned char ug = g;
                Write<unsigned char>(it, ug);
                Write<unsigned char>(it, ug);
                Write<unsigned char>(it, ug);
            }
            else if(params.add_rgb == 2) // quality mode
            {
                Write<float>(it, p_refl->at(block_idx)[last_echo_idx]);
            }
        }
        mls.Free(block_idx);
    }
    cout << "Writing " << n_vertex << " vertices of size " << vertex_buffer_size << "=" << 1.e-6*vertex_buffer_size << "MB" << endl;
    fileOut.write(buffer, vertex_buffer_size);
    delete buffer;
    // triangle list
    unsigned int tri_bytesize = sizeof(unsigned char)+3*sizeof(int);
    unsigned long tri_buffer_size = tri_bytesize*v_tri.size();
    buffer = new char[tri_buffer_size]; it = buffer;
    for(auto & tri:v_tri)
    {
        Write<unsigned char>(it, 3);
        Write<int>(it, tri.i);
        Write<int>(it, tri.j);
        Write<int>(it, tri.k);
    }
    cout << "Writing " << v_tri.size() << " triangles of size " << tri_buffer_size << "=" << 1.e-6*tri_buffer_size << "MB" << endl;
    cout << "Total " << 1.e-6*(vertex_buffer_size+tri_buffer_size) << "MB" << endl;
    fileOut.write(buffer, tri_buffer_size);
    fileOut.close();
}

void WriteOff(XMls & mls, vector<Triangle> v_tri,
              int n_vertex, param params)
{
    ofstream fileOut(params.output);
    if(!fileOut.good())
    {
        cout << "Cannot open " + params.output + " for writing\n";
        return;
    }
    // write text header
    fileOut << "OFF" << endl;
    fileOut << "# Generated from " << params.xmls_param.ept_folder << " secs " << params.i_start << " to " << params.i_end << endl;
    fileOut << "# IGN offset Pos " << params.pivot_E << " " << params.pivot_N << " " << params.pivot_H << endl << endl;
    // number of edges is not used in OFF but Meshlab fails if absent. It is long to compute so we provide an estimate (exact for watertight meshes)
    fileOut << n_vertex << " " << v_tri.size() << " " << v_tri.size() << endl;
    // vertex list
    for(XBlockIndex block_idx=0; block_idx<mls.NBlock(); block_idx++)
    {
        mls.Load(block_idx);
        for(XPulseIndex pulse_idx=0; pulse_idx<mls.NPulse(block_idx); pulse_idx++) if(mls.NbOfEcho(block_idx, pulse_idx)>0)
        {
            XPt3D Pw = mls.Pworld(block_idx, mls.IdxLastEcho(block_idx, pulse_idx));
            float e=Pw.X-params.pivot_E, n=Pw.Y-params.pivot_N, h=Pw.Z;
            fileOut << e << " " << n << " " << h << endl;
        }
        mls.Free(block_idx);
    }

    // triangle list
    for(auto & tri:v_tri) fileOut << "3 " << tri.i << " " << tri.j << " " << tri.k << endl;
    fileOut.close();
}

float MaxEdgeSize(XMls & mls, XBlockIndex b1, XPulseIndex id1, XBlockIndex b2, XPulseIndex id2, XBlockIndex b3, XPulseIndex id3)
{
    XPt3D P1 = mls.Pworld(b1, mls.IdxLastEcho(b1, id1));
    XPt3D P2 = mls.Pworld(b2, mls.IdxLastEcho(b2, id2));
    XPt3D P3 = mls.Pworld(b3, mls.IdxLastEcho(b3, id3));
    double d12 = (P1-P2).Norme();
    double d23 = (P2-P3).Norme();
    double d31 = (P3-P1).Norme();
    return max(d12,max(d23,d31));
}

float MaxEdgeSize(XMls & mls, XBlockIndex b, XPulseIndex id1, XPulseIndex id2, XPulseIndex id3)
{
    return MaxEdgeSize(mls, b, id1, b, id2, b, id3);
}

// make a wedge based on the four vertices ABXY = (AXB)+(ABY) triangles
void MakeWedge(XMls & mls, vector<Triangle> & v_tri, vector<vector<int> > & vv_pulse_with_echo_idx, float tri_threshold,
               XBlockIndex ba, XPulseIndex pa, XBlockIndex bb, XPulseIndex pb,
               XBlockIndex bx, XPulseIndex px, XBlockIndex by, XPulseIndex py)
{
    if(mls.NbOfEcho(ba, pa)==0 || mls.NbOfEcho(bb, pb) == 0) return;
    //cout << "ABXY (" << ba << "," << pa << ")("<< bb << "," << pb << ")("<< bx << "," << px << ")("<< by << "," << py << ")" << endl;
    int A=vv_pulse_with_echo_idx[ba][pa], B=vv_pulse_with_echo_idx[bb][pb];
    int X=vv_pulse_with_echo_idx[bx][px], Y=vv_pulse_with_echo_idx[by][py];
    if(mls.NbOfEcho(bx, px)>0)
    {
        if(MaxEdgeSize(mls, ba, pa, bb, pb, bx, px) < tri_threshold)
            v_tri.push_back(Triangle(A, B, X));
    }
    if(mls.NbOfEcho(by, py)>0)
    {
        if(MaxEdgeSize(mls, ba, pa, bb, pb, by, py) < tri_threshold)
            v_tri.push_back(Triangle(A, Y, B));
    }
}

// mono block version
inline void MakeWedge(XMls & mls, vector<Triangle> & v_tri, vector<vector<int> > & vv_pulse_with_echo_idx, float tri_threshold,
                      XBlockIndex b, XPulseIndex pa, XPulseIndex pb, XPulseIndex px, XPulseIndex py)
{
    MakeWedge(mls, v_tri, vv_pulse_with_echo_idx, tri_threshold, b, pa, b, pb, b, px, b, py);
}

//-----------------------------------------------------------------------------
int main(int argc, char **argv)
{
    cout << WHAT << endl;
    param params;
    if(argc < 6)
    {
        cout << "Usage: " << argv[0] << "  sbet acc laser_calib.xml ept_folder output_folder [block_length dist_threshold tri_threshold pivot_E pivot_N pivot_H add_rgb]" << endl;
        cout << "sbet: sbet file" << endl;
        cout << "acc: accuracy file (smrmsg)" << endl;
        cout << "laser_calib.xml: calibration file for the laser" << endl;
        cout << "ept_folder: folder containing the echo pulse tables (generated with EptExport)" << endl;
        cout << "output: name of the output mesh file (.ply or .off only)" << endl;
        cout << "block_length: expected length of an exported block on min distance between lines to keep it (default=" << params.block_length << "m)" << endl;
        cout << "dist_threshold: threshold on min distance between lines to keep it (default=" << params.dist_threshold << "m)" << endl;
        cout << "tri_threshold: threshold on max triangle edge size to add the triangle (default=" << params.tri_threshold << "m)" << endl;
        cout << "pivot_(E|N|H): optionally set manually the pivot point (default=-1=first point of the trajectory rounded at 100m for EN, 0 for H)" << endl;
        cout << "add_rgb: if not 0, add r,g, b attributes in .ply mode based on reflectance (default=0=don't)" << endl;
        return 0;
    }
    int i_arg=1;

    // required
    params.xmls_param.sbet = string((argv[i_arg++]));
    params.xmls_param.acc = string(argv[i_arg++]);
    params.xmls_param.laser_calib = string(argv[i_arg++]);
    params.xmls_param.ept_folder = string(argv[i_arg++]);
    params.output = string(argv[i_arg++]);
    params.block_length = atof(argv[i_arg++]);
    params.dist_threshold = atof(argv[i_arg++]);

    // optional
    if(i_arg < argc) params.tri_threshold = atof(argv[i_arg++]);
    //if(i_arg < argc) params.max_DP_error = atof(argv[i_arg++]);
    if(i_arg < argc) params.pivot_E = atoi(argv[i_arg++]);
    if(i_arg < argc) params.pivot_N = atoi(argv[i_arg++]);
    if(i_arg < argc) params.pivot_H = atoi(argv[i_arg++]);
    if(i_arg < argc) params.add_rgb = atoi(argv[i_arg++]);


    clock_t start = clock();
    // constructor and infos accessible after construction
    XMls mls(params.xmls_param);
    mls.Select(params.i_start, params.i_end);
    cout << mls.NBlock() << " block(s) selected" << endl;

    if(params.pivot_E == -1)
    {
        XPt3D Pivot = mls.m_trajecto.SbetSeries().GetGeoref(0).Translation();
        params.pivot_E = 100*(int)(Pivot.X/100);
        params.pivot_N = 100*(int)(Pivot.Y/100);
    }

    start = clock();
    // lines are defined as sets of successive pulses with increasing theta
    // a new line starts when theta finishes a 2\pi rotation (goes back to 0)
    double d_sum=0., d_min=1.e9, d_max=0., d_cur=0; // distance between line, d_cur is distance to the last kept line
    vector< vector<int> > vv_next_line(mls.NBlock()); // number of lines to increment, -1 if on a skipped line
    for(XBlockIndex block_idx=0; block_idx<mls.NBlock(); block_idx++)
    {
        int next_line=1;
        mls.Load(block_idx);
        cout << "Block " << block_idx << ", n_pulse=" << mls.NPulse(block_idx) << endl;
        vv_next_line[block_idx] = vector<int>(mls.NPulse(block_idx), -1);
        for(XPulseIndex pulse_idx=0; pulse_idx<mls.NPulse(block_idx)-1; pulse_idx++)
        {
            if(mls.Theta(block_idx, pulse_idx+1)<mls.Theta(block_idx, pulse_idx)) // reached the next line
            {
                next_line=1;
                do
                {
                    XPulseIndex next_line_idx = pulse_idx + next_line*mls.PulsePerLine();
                    bool ok = (next_line_idx < mls.NPulse(block_idx));
                    if(ok)
                    {
                        XPt3D O = mls.Ins(block_idx, pulse_idx).Translation();
                        XPt3D O_next = mls.Ins(block_idx, next_line_idx).Translation();
                        double d = (O_next-O).Norme();
                        d_cur += d;
                        next_line++;
                        // stats
                        d_sum += d;
                        if(d<d_min) d_min = d;
                        if(d>d_max) d_max = d;
                    }
                } while(ok && d_cur < params.dist_threshold);
            }
            vv_next_line[block_idx][pulse_idx] = next_line;

        }
        mls.Free(block_idx);
    }

    // indexation of kept echoes (non last echoes and lines too closed are removed)
    int PPL = mls.PulsePerLine();
    int idx=0;
    vector< vector<int> > vv_pulse_with_echo_idx(mls.NBlock()); // unique indexation of pulses with at least one echo across blocks
    for(XBlockIndex block_idx=0; block_idx<mls.NBlock(); block_idx++)
    {
        mls.Load(block_idx);
        // index pulses with at leat one echo
        vv_pulse_with_echo_idx[block_idx] = vector<int>(mls.NPulse(block_idx), -1);
        for(XPulseIndex pulse_idx=0; pulse_idx<mls.NPulse(block_idx); pulse_idx++)
            if(mls.NbOfEcho(block_idx, pulse_idx)>0)
                vv_pulse_with_echo_idx[block_idx][pulse_idx] = idx++;
        mls.Free(block_idx);
    }

    // create triangles list
    vector<Triangle> v_tri;
    mls.Load(0);
    for(XBlockIndex block_idx=0; block_idx<mls.NBlock(); block_idx++)
    {
        XPulseIndex n_pulse = mls.NPulse(block_idx);
        for(XPulseIndex pulse_idx=0; pulse_idx<n_pulse-PPL-1; pulse_idx++)
        {
            MakeWedge(mls, v_tri, vv_pulse_with_echo_idx, params.tri_threshold,
                      block_idx, pulse_idx, pulse_idx+PPL+1, pulse_idx+PPL, pulse_idx+1);
        }
        // triangles between current block and next one
        if(block_idx<mls.NBlock()-1)
        {
            mls.Load(block_idx+1);
            cout << "PPL " << PPL << " n_pulse " << n_pulse << " block " << block_idx << "/" << mls.NBlock() << endl;
            // first wedge (3 vertices in block, 1 in block+1)
            MakeWedge(mls, v_tri, vv_pulse_with_echo_idx, params.tri_threshold,
                      block_idx, n_pulse-PPL-1, block_idx+1, 0,
                      block_idx, n_pulse-1, block_idx, n_pulse-PPL);
            for(XPulseIndex pulse_idx=n_pulse-PPL; pulse_idx<n_pulse-1; pulse_idx++)
            {
                // strip wedges (2 vertices in block, 2 in block+1)
                MakeWedge(mls, v_tri, vv_pulse_with_echo_idx, params.tri_threshold,
                          block_idx, pulse_idx, block_idx+1, pulse_idx+PPL+1-n_pulse,
                          block_idx+1, pulse_idx+PPL-n_pulse, block_idx, pulse_idx+1);
            }
            // last wedge (1 vertex in block, 3 in block+1)
            MakeWedge(mls, v_tri, vv_pulse_with_echo_idx, params.tri_threshold,
                      block_idx, n_pulse-1, block_idx+1, PPL,
                      block_idx+1, PPL-1, block_idx+1, 0);
        }
        mls.Free(block_idx);
    }
    cout << "Creating " << params.output << endl;
    if(params.Format() == "off") WriteOff(mls, v_tri, idx, params);
    else WritePly(mls, v_tri, idx, params);
    cout << "Done in " << (double)(clock() - start) / (double) CLOCKS_PER_SEC << "s" << endl;
    return 0;
}

