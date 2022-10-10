
#define WHAT "SensorMesh3: create a single .ply mesh for all the specified time interval using sensor topology and filtering out accumulated lines"

#include <ctime>
#include <iostream>
#include <limits>
#include "libXMls/XMls.h"
#include "libXBase/XPt2D.h"


using namespace std;

struct param
{
    //string sbet, acc;
    XMls::Param xmls_param;
    string output, format;
    double block_length;

    // optional
    float dist_threshold;
    float tri_threshold;
    float max_DP_error;
    XSecond i_start, i_end;
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

// A (x,y,z) point with an index
struct PointIdx
{
    XPt3D P;
    int idx;
    PointIdx(float x=0, float y=0, float z=0, unsigned int idx_=-1):
        P(x, y, z), idx(idx_) {}
    PointIdx(XPt3D P_, unsigned int idx_=-1):
        P(P_), idx(idx_) {}
    inline bool Exists(){return (idx > -1);}
};

// A (x,y,z) point and a center point with an index
struct Point0Idx
{
    XPt3D P, P0;
    int idx;
    Point0Idx(float x=0, float y=0, float z=0,
             float x0=0, float y0=0, float z0=0, unsigned int idx_=-1):
        P(x, y, z), P0(x0, y0, z0), idx(idx_) {}
    Point0Idx(XPt3D P_, XPt3D P0_, unsigned int idx_=-1):
        P(P_), P0(P0_), idx(idx_) {}
    inline bool Exists(){return (idx > -1);}
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
              int n_vertex, param params, vector<vector<int>> vv_pulse_with_echo_idx, XBlockIndex debut, XBlockIndex fin)
{
    cout << "n_vertex = " << n_vertex << endl;
    cout << "v_tri.size() = " << v_tri.size() << endl;
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
    for(XBlockIndex block_idx=debut; block_idx<fin; block_idx++)
    {
        mls.Load(block_idx);
        cout << "block_idx_ply = " <<block_idx << endl;
        for(XPulseIndex pulse_idx=0; pulse_idx<mls.NPulse(block_idx); pulse_idx++)
        {
            if(vv_pulse_with_echo_idx[block_idx][pulse_idx] != -1)
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
        }
        mls.Free(block_idx);
        //cout << "\n i am here now!!" << endl;
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

void WritePly(vector<Point0Idx> & v_ptidx, vector<Triangle> & v_tri,
              int n_vertex, param params, string filename)
{
    ofstream fileOut(filename);
    if(!fileOut.good())
    {
        cout << "Cannot open " + filename + " for writing\n";
        return;
    }

    // write text header
    fileOut << "ply\nformat binary_little_endian 1.0" << endl;
    fileOut << "comment Generated from " << params.xmls_param.ept_folder << " secs " << 0 << " to " << 0 << endl; // TODO
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
    } else if(params.add_rgb == 3)
    {
        fileOut << "property float x0" << endl;
        fileOut << "property float y0" << endl;
        fileOut << "property float z0" << endl;
    }
    fileOut << "element face " << v_tri.size() << endl;
    fileOut << "property list uchar int vertex_indices" << endl;
    fileOut << "end_header" << endl;
    // vertex list
    unsigned int vertex_bytesize = 3*sizeof(float);
    unsigned long vertex_buffer_size = vertex_bytesize * n_vertex;
    //XFloatAttrib * p_refl=NULL;
    //if(params.add_rgb > 0) p_refl = mls.GetEchoAttrib<XFloatAttrib>("reflectance");
    if(params.add_rgb == 1) vertex_buffer_size += 3*sizeof(unsigned char) * n_vertex;
    else if(params.add_rgb == 2) vertex_buffer_size += sizeof(float) * n_vertex;
    else if(params.add_rgb == 3) vertex_buffer_size += 3*sizeof(float) * n_vertex;

    char * buffer = new char[vertex_buffer_size], * it = buffer;
    for(auto & ptidx:v_ptidx) if(ptidx.Exists())
    {
        Write<float>(it, ptidx.P.X);
        Write<float>(it, ptidx.P.Y);
        Write<float>(it, ptidx.P.Z);
        if(params.add_rgb == 1) // r g b mode
        {
            /*float g = 12.75f*(p_refl->at(block_idx)[last_echo_idx]+20.f); // rescale [-20,0] to [0,256] TODO: parameters
            if(g<0.f) g=0.f; else if(g>255.f) g=255.f;
            unsigned char ug = g;
            Write<unsigned char>(it, ug);
            Write<unsigned char>(it, ug);
            Write<unsigned char>(it, ug);*/
        }
        else if(params.add_rgb == 2) // quality mode
        {
            //Write<float>(it, p_refl->at(block_idx)[last_echo_idx]);
        }
        else if(params.add_rgb == 3) // center mode
        {
            Write<float>(it, ptidx.P0.X);
            Write<float>(it, ptidx.P0.Y);
            Write<float>(it, ptidx.P0.Z);
        }
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
              int n_vertex, param params, vector<vector<int>> vv_pulse_with_echo_idx, XBlockIndex debut, XBlockIndex fin)
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
    for(XBlockIndex block_idx=debut; block_idx<fin; block_idx++)
    {
        mls.Load(block_idx);
        for(XPulseIndex pulse_idx=0; pulse_idx<mls.NPulse(block_idx); pulse_idx++) if(vv_pulse_with_echo_idx[block_idx][pulse_idx] != -1)
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

float MaxEdgeSize(XPt3D P1, XPt3D P2, XPt3D P3)
{
    double d12 = (P1-P2).Norme();
    double d23 = (P2-P3).Norme();
    double d31 = (P3-P1).Norme();
    return max(d12,max(d23,d31));
}

float MaxEdgeSize(XMls & mls, XBlockIndex b1, XPulseIndex id1, XBlockIndex b2, XPulseIndex id2, XBlockIndex b3, XPulseIndex id3)
{
    XPt3D P1 = mls.Pworld(b1, mls.IdxLastEcho(b1, id1));
    XPt3D P2 = mls.Pworld(b2, mls.IdxLastEcho(b2, id2));
    XPt3D P3 = mls.Pworld(b3, mls.IdxLastEcho(b3, id3));
    return MaxEdgeSize(P1, P2, P3);
}

float MaxEdgeSize(XMls & mls, XBlockIndex b, XPulseIndex id1, XPulseIndex id2, XPulseIndex id3)
{
    return MaxEdgeSize(mls, b, id1, b, id2, b, id3);
}

float MaxEdgeSize(vector<Point0Idx> & v_ptidx, int i1, int i2, int i3)
{
    XPt3D P1 = v_ptidx[i1].P;
    XPt3D P2 = v_ptidx[i2].P;
    XPt3D P3 = v_ptidx[i3].P;
    return MaxEdgeSize(P1, P2, P3);
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

vector<Triangle> MakeTriangles(vector<Point0Idx> & v_ptidx, int PPL, float tri_threshold)
{
    vector<Triangle> v_tri;
    for(int i=0; i<(int)v_ptidx.size()-PPL-1; i++)
    {
        // make the 2 triangles based on i, i+PPL, i+1, I+PPL+1 = (i, i+1, i+PPL) and (i+1, i+PPL+1, i+PPL)
        if(v_ptidx[i+1].Exists() && v_ptidx[i+PPL].Exists())
        {
            if(v_ptidx[i].Exists())
            {
                if(MaxEdgeSize(v_ptidx, i, i+1, i+PPL) < tri_threshold)
                    v_tri.push_back(Triangle(v_ptidx[i].idx, v_ptidx[i+1].idx, v_ptidx[i+PPL].idx));
            }
            if(v_ptidx[i+PPL+1].Exists())
            {
                if(MaxEdgeSize(v_ptidx, i+1, i+PPL+1, i+PPL) < tri_threshold)
                    v_tri.push_back(Triangle(v_ptidx[i+1].idx, v_ptidx[i+PPL+1].idx, v_ptidx[i+PPL].idx));
            }
        }
    }
    return v_tri;
}

//-----------------------------------------------------------------------------
int main(int argc, char **argv)
{
    cout << WHAT << endl;
    param params;
    if(argc < 9)
    {
        cout << "Usage: " << argv[0] << "  sbet laser_calib.xml ept_folder output i_start i_end ";
        cout << "[block_length dist_threshold tri_threshold pivot_E pivot_N pivot_H add_rgb]" << endl;
        cout << "sbet: sbet file" << endl;
        cout << "acc: accuracy file (smrmsg)" << endl;
        cout << "laser_calib.xml: calibration file for the laser" << endl;
        cout << "ept_folder: folder containing the echo pulse tables (generated with EptExport)" << endl;
        cout << "output: name of the output folder" << endl;
        cout << "block_length: expected length of an exported block on min distance between lines to keep it (default="
             << params.block_length << "m)" << endl;
        cout << "dist_threshold: threshold on min distance between lines to keep it (default="
             << params.dist_threshold << "m)" << endl;
        cout << "tri_threshold: threshold on max triangle edge size to add the triangle (default="
             << params.tri_threshold << "m)" << endl;
        cout << "pivot_(E|N|H): optionally set manually the pivot point "
             << "(default=-1=first point of the trajectory rounded at 100m for EN, 0 for H)" << endl;
        cout << "add_rgb: (0:add nothing, 1: add reflectance in r,g,b attributes, 2:add reflectance in quality, "
             << "3:add center coords) in .ply (default=0)" << endl;
        return 0;
    }
    int i_arg=1;

    // required
    string sbet(argv[i_arg++]);
    string laser_calib(argv[i_arg++]);
    string ept_folder(argv[i_arg++]);
    params.output = string(argv[i_arg++]);
    params.i_start = atof(argv[i_arg++]);
    params.i_end = atof(argv[i_arg++]);

    // optional
    if(i_arg < argc) params.block_length = atof(argv[i_arg++]);
    if(i_arg < argc) params.dist_threshold = atof(argv[i_arg++]);
    if(i_arg < argc) params.tri_threshold = atof(argv[i_arg++]);
    //if(i_arg < argc) params.max_DP_error = atof(argv[i_arg++]);
    if(i_arg < argc) params.pivot_E = atoi(argv[i_arg++]);
    if(i_arg < argc) params.pivot_N = atoi(argv[i_arg++]);
    if(i_arg < argc) params.pivot_H = atoi(argv[i_arg++]);
    if(i_arg < argc) params.add_rgb = atoi(argv[i_arg++]);


    clock_t start = clock();
    // constructor and infos accessible after construction
    XTrajecto traj(sbet);
    XMls mls(ept_folder, laser_calib, &traj);
    mls.Select(params.i_start, params.i_end);
    cout << mls.NBlock() << " block(s) selected" << endl;

    if(params.pivot_E == -1)
    {
        XPt3D Pivot = traj.GetGeoref(0).Translation();
        params.pivot_E = 100*(int)(Pivot.X/100);
        params.pivot_N = 100*(int)(Pivot.Y/100);
    }

    start = clock();

    // lines are defined as sets of successive pulses with increasing theta
    // a new line starts when theta finishes a 2\pi rotation (goes back to 0)
    // note: this definition is arbitrary because in fact there is a single helicoidal line
    double d_sum=0., d_min=1.e9, d_max=0., d_mesh=0; // distance between line, d is distance to the last kept line
    XPt3D O=mls.Ins(0, 0).Translation(); // last trajectory point corresponding to a kept line
    vector<Point0Idx> v_ptidx; // stores the kept point (not on removed lines)
    bool keeping = true, first=true; // flag to know if we are currently on a skipped line or not, and if we are on the first line or not (we always keep first line)
    const int PPL = mls.PulsePerLine();
    int mesh_start_time=params.i_start, mesh_end_time=0; // integer time of the start and end of mesh
    int idx=0, n_line=0, n_skip=0, n_line_mesh=0, n_skip_mesh=0; // existing point index, number of (skipped) lines (per mesh)
    for(XBlockIndex block_idx=0; block_idx<mls.NBlock(); block_idx++)
    {
        mls.Load(block_idx);
        cout << "Extracting points from Block " << block_idx << ", n_pulse=" << mls.NPulse(block_idx) << endl;
        for(XPulseIndex pulse_idx=0; pulse_idx<mls.NPulse(block_idx)-1; pulse_idx++)
        {
            if(keeping)
            {
                if(mls.NbOfEcho(block_idx, pulse_idx)>0)
                {
                    XEchoIndex last_echo_idx = mls.IdxLastEcho(block_idx, pulse_idx);
                    XPt3D Pw = mls.Pworld(block_idx, last_echo_idx);
                    XPt3D Cw = mls.Cworld(block_idx, last_echo_idx);
                    float e=Pw.X-params.pivot_E, n=Pw.Y-params.pivot_N, h=Pw.Z-params.pivot_H;
                    float e0=Cw.X-params.pivot_E, n0=Cw.Y-params.pivot_N, h0=Cw.Z-params.pivot_H;
                    v_ptidx.push_back(Point0Idx(e, n, h, e0, n0, h0, idx++));
                }
                else v_ptidx.push_back(Point0Idx());
            }
            if(mls.Theta(block_idx, (XPulseIndex)(pulse_idx+1))<mls.Theta(block_idx, (XPulseIndex)pulse_idx)) // reached the next line
            {
                n_line_mesh++;
                XPt3D O_next = mls.Ins(block_idx, pulse_idx).Translation();
                double d = (O_next-O).Norme();
                // stats
                d_sum += d;
                d_mesh += d;
                //cout << "Pulse " << pulse_idx << " d=" << d << " d_sum=" << d_sum << " d_mesh=" << d_mesh << endl;
                if(d<d_min) d_min = d;
                if(d>d_max) d_max = d;
                O = O_next;
                if(first || d_sum > params.dist_threshold)
                {
                    keeping=true;
                    d_sum=0;
                    first = false;
                } else
                {
                    keeping=false;
                    n_skip_mesh++;
                }
                if (d_mesh > params.block_length)
                {
                    mesh_end_time=(int)mls.Time(block_idx, pulse_idx);
                    ostringstream oss;
                    oss << params.output << "/" << mesh_start_time << "-" << mesh_end_time << ".ply";
                    cout << "Creating mesh " << oss.str() << " with " <<
                            n_skip_mesh << "/" << n_line_mesh << " lines skipped" << endl;
                    vector<Triangle> v_tri = MakeTriangles(v_ptidx, PPL, params.tri_threshold);
                    WritePly(v_ptidx, v_tri, idx, params, oss.str());
                    n_line+=n_line_mesh;
                    n_skip+=n_skip_mesh;
                    d_mesh=0;
                    n_line_mesh=0;
                    n_skip_mesh=0;
                    v_ptidx.clear();
                    idx = 0;
                    mesh_start_time = mesh_end_time;
                }
            }
        }
        mls.Free(block_idx);
    }
    // create the last mesh
    mesh_end_time=params.i_end;
    ostringstream oss;
    oss << params.output << "/" << mesh_start_time << "-" << mesh_end_time << ".ply";
    cout << "Creating mesh " << oss.str() << " with " <<
            n_skip_mesh << "/" << n_line_mesh << " lines skipped" << endl;
    vector<Triangle> v_tri = MakeTriangles(v_ptidx, PPL, params.tri_threshold);
    WritePly(v_ptidx, v_tri, idx, params, oss.str());

    cout << "Done in " << (double)(clock() - start) / (double) CLOCKS_PER_SEC << "s" << endl;

    // Stats
    cout << "d_min=" << d_min << " d_max=" << d_max << endl;
    cout << n_skip << "/" << n_line << " lines skipped" << endl;

    return 0;

}


