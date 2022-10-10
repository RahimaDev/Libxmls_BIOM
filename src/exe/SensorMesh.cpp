
#define WHAT "SensorMesh: create a .ply mesh per block using sensor topology"

#include <ctime>
#include <iostream>
#include <limits>
#include "libXMls/XMls.h"
#include "libXBase/XPt2D.h"


using namespace std;

struct Triangle
{
    unsigned int i,j,k;
    Triangle(unsigned int i_=0, unsigned int j_=0, unsigned int k_=0):i(i_),j(j_),k(k_){}
};

struct DP_point
{
    XPt3D P;
    XPulseIndex pulse_idx;
    XEchoIndex echo_idx;
    int dec_idx; // indexation of decimated points
    float theta;
    bool keep;
};

template <typename T> void Write(char * & it, T data)
{
    *reinterpret_cast<T*>(it) = data;
    it += sizeof(T);
}

float MaxEdgeSize(XMls & mls, XBlockIndex block_idx, XPulseIndex id1, XPulseIndex id2, XPulseIndex id3)
{
    XPt3D P1 = mls.Pworld(block_idx, mls.IdxLastEcho(block_idx, id1));
    XPt3D P2 = mls.Pworld(block_idx, mls.IdxLastEcho(block_idx, id2));
    XPt3D P3 = mls.Pworld(block_idx, mls.IdxLastEcho(block_idx, id3));
    double d12 = (P1-P2).Norme();
    double d23 = (P2-P3).Norme();
    double d31 = (P3-P1).Norme();
    return max(d12,max(d23,d31));
}

void WritePly(XMls & mls, string ept_folder, XBlockIndex block_idx,
              vector<bool> has_echo, vector<Triangle> v_tri,
              string filename, int n_vertex,
              int pivot_E, int pivot_N)
{
    ofstream fileOut(filename);
    if(!fileOut.good())
    {
        cout << "Cannot open " + filename + " for writing\n";
        return;
    }

    // write text header
    fileOut << "ply\nformat binary_little_endian 1.0" << endl;
    fileOut << "comment Generated from " << ept_folder << " block " << block_idx << endl;
    fileOut << "comment IGN offset Pos " << pivot_E << " " << pivot_N << " 0" << endl;
    fileOut << "element vertex " << n_vertex << endl;
    fileOut << "property float x" << endl;
    fileOut << "property float y" << endl;
    fileOut << "property float z" << endl;
    fileOut << "element face " << v_tri.size() << endl;
    fileOut << "property list uchar int vertex_indices" << endl;
    fileOut << "end_header" << endl;
    // vertex list
    unsigned int vertex_bytesize = 3*sizeof(float);
    unsigned long vertex_buffer_size = vertex_bytesize * n_vertex;
    char * buffer = new char[vertex_buffer_size], * it = buffer;
    for(XPulseIndex pulse_idx=0; pulse_idx<mls.NPulse(block_idx); pulse_idx++) if(has_echo[pulse_idx])
    {
        XPt3D Pw = mls.Pworld(block_idx, mls.IdxLastEcho(block_idx, pulse_idx));
        float e=Pw.X-pivot_E, n=Pw.Y-pivot_N, h=Pw.Z;
        Write<float>(it, e);
        Write<float>(it, n);
        Write<float>(it, h);
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

void WritePly(XMls & mls, string ept_folder, XBlockIndex block_idx,
              vector<DP_point> v_pts, vector<Triangle> v_tri, string filename,
              int pivot_E, int pivot_N)
{
    ofstream fileOut(filename);
    if(!fileOut.good())
    {
        cout << "Cannot open " + filename + " for writing\n";
        return;
    }

    // write text header
    fileOut << "ply\nformat binary_little_endian 1.0" << endl;
    fileOut << "comment Generated from " << ept_folder << " block " << block_idx << endl;
    fileOut << "comment IGN offset Pos " << pivot_E << " " << pivot_N << " 0" << endl;
    fileOut << "element vertex " << v_pts.size() << endl;
    fileOut << "property float x" << endl;
    fileOut << "property float y" << endl;
    fileOut << "property float z" << endl;
    fileOut << "element face " << v_tri.size() << endl;
    fileOut << "property list uchar int vertex_indices" << endl;
    fileOut << "end_header" << endl;
    // vertex list
    unsigned int vertex_bytesize = 3*sizeof(float);
    unsigned long vertex_buffer_size = vertex_bytesize * v_pts.size();
    char * buffer = new char[vertex_buffer_size], * it = buffer;
    for(auto & pt:v_pts)
    {
        XPt3D Pw = mls.Pworld(block_idx, pt.echo_idx);
        float e=Pw.X-pivot_E, n=Pw.Y-pivot_N, h=Pw.Z;
        Write<float>(it, e);
        Write<float>(it, n);
        Write<float>(it, h);
    }
    cout << "Writing " << v_pts.size() << " vertices of size " << vertex_buffer_size << "=" << 1.e-6*vertex_buffer_size << "MB" << endl;
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

void WriteOff(XMls & mls, string ept_folder, XBlockIndex block_idx,
              vector<bool> has_echo, vector<Triangle> v_tri,
              string filename, int n_vertex,
              int pivot_E, int pivot_N)
{
    ofstream fileOut(filename);
    if(!fileOut.good())
    {
        cout << "Cannot open " + filename + " for writing\n";
        return;
    }
    // write text header
    fileOut << "OFF" << endl;
    fileOut << "# Generated from " << ept_folder << " block " << block_idx << endl;
    fileOut << "# IGN offset Pos " << pivot_E << " " << pivot_N << " 0" << endl << endl;
    // number of edges is not used in OFF but Meshlab fails if absent. It is long to compute so we provide an estimate (exact for watertight meshes)
    fileOut << n_vertex << " " << v_tri.size() << " " << v_tri.size() << endl;
    // vertex list
    for(XPulseIndex pulse_idx=0; pulse_idx<mls.NPulse(block_idx); pulse_idx++) if(has_echo[pulse_idx])
    {
        XPt3D Pw = mls.Pworld(block_idx, mls.IdxLastEcho(block_idx, pulse_idx));
        float e=Pw.X-pivot_E, n=Pw.Y-pivot_N, h=Pw.Z;
        fileOut << e << " " << n << " " << h << endl;
    }
    // triangle list
    for(auto & tri:v_tri) fileOut << "3 " << tri.i << " " << tri.j << " " << tri.k << endl;
    fileOut.close();
}

void WriteOff(XMls & mls, string ept_folder, XBlockIndex block_idx,
              vector<DP_point> v_pts, vector<Triangle> v_tri,
              string filename,
              int pivot_E, int pivot_N)
{
    ofstream fileOut(filename);
    if(!fileOut.good())
    {
        cout << "Cannot open " + filename + " for writing\n";
        return;
    }
    // write text header
    fileOut << "OFF" << endl;
    fileOut << "# Generated from " << ept_folder << " block " << block_idx << endl;
    fileOut << "# IGN offset Pos " << pivot_E << " " << pivot_N << " 0" << endl << endl;
    // number of edges is not used in OFF but Meshlab fails if absent. It is long to compute so we provide an estimate (exact for watertight meshes)
    fileOut << v_pts.size() << " " << v_tri.size() << " " << v_tri.size() << endl;
    // vertex list
    for(auto & pt:v_pts)
    {
        XPt3D Pw = mls.Pworld(block_idx, pt.echo_idx);
        float e=Pw.X-pivot_E, n=Pw.Y-pivot_N, h=Pw.Z;
        fileOut << e << " " << n << " " << h << endl;
    }
    // triangle list
    for(auto & tri:v_tri) fileOut << "3 " << tri.i << " " << tri.j << " " << tri.k << endl;
    fileOut.close();
}

// recursive Douglas-Peucker function decimating (P_i_min,...,P_i_max).
// Computes the approximation of input in the given range as a subset of input.
// decimation if performed by setting the keep member of the points for points that should be kept in the decimated result
void DPDecimate(vector<DP_point> & input, float max_DP_error, int i_min=0, int i_max=0)
{
    if(0 == i_max) i_max=input.size()-1;
    input.at(i_min).keep=true;
    input.at(i_max).keep=true;
    float d_max=0.;
    int i_star=i_min;
    XPt3D AB = input.at(i_max).P-input.at(i_min).P;
    AB.Normalise();
    for(int i=i_min+1; i<i_max; i++)
    {
        XPt3D AC = input.at(i).P-input.at(i_min).P;
        float d = fabs(AB.Y*AC.X-AB.X*AC.Y);
        if(d > d_max)
        {
            d_max=d;
            i_star=i;
        }
    }
    if(d_max > max_DP_error)
    {
        DPDecimate(input, max_DP_error, i_min, i_star);
        DPDecimate(input, max_DP_error, i_star, i_max);
    }
}

//-----------------------------------------------------------------------------
int main(int argc, char **argv)
{
    cout << WHAT << endl;
    if(argc < 5)
    {
        cout << "Usage: " << argv[0] << "  sbet acc laser_calib.xml ept_folder output_folder [tri_threshold=0.5 max_DP_error=0 format=ply (int)start_time (int)end_time]" << endl;
        cout << "sbet: path to the sbet (trajecto) file" << endl;
        cout << "laser_calib.xml: calibration file for the laser" << endl;
        cout << "ept_folder: folder containing the echo pulse tables (generated with EptExport)" << endl;
        cout << "tri_threshold: threshold on max triangle edge size to add the triangle (default=0.5m)" << endl;
        cout << "max_DP_error: maximum Douglas-Peucker error (default=0m=do not decimate)" << endl;
        cout << "format: choose between: ply (default), off" << endl;
        cout << "(start|end)_time: start and end time of the laser points to export (default=everything)" << endl;
        return 0;
    }

    int i_arg=1;
    // required
    string sbet(argv[i_arg++]);
    string laser_calib(argv[i_arg++]);
    string ept_folder(argv[i_arg++]);
    string output_folder(argv[i_arg++]);

    // optional
    float tri_threshold=0.5;
    if(i_arg < argc) tri_threshold = atof(argv[i_arg++]);
    float max_DP_error=0;
    if(i_arg < argc) max_DP_error = atof(argv[i_arg++]);
    string format = "ply";
    if(i_arg < argc) format = string(argv[i_arg++]);
    XSecond i_start = -1, i_end = -1;
    if(i_arg < argc) i_start = atoi(argv[i_arg++]);
    if(i_arg < argc) i_end = atoi(argv[i_arg++]);


    clock_t start = clock();
    // constructor and infos accessible after construction
    XTrajecto traj(sbet);
    XMls mls(ept_folder, laser_calib, &traj);
    mls.Select(i_start, i_end);
    cout << mls.NBlock() << " block(s) selected" << endl;

    XPt3D Pivot = traj.GetGeoref(0).Translation();
    int pivot_E = 100*(int)(Pivot.X/100);
    int pivot_N = 100*(int)(Pivot.Y/100);

    start = clock();
    int PPL = mls.PulsePerLine();
    for(XBlockIndex block_idx=0; block_idx<mls.NBlock(); block_idx++)
    {
        mls.Load(block_idx);
        ostringstream oss;
        oss << output_folder << "/" << mls.Second(block_idx) << "." << format;
        cout << "Creating " << oss.str() << ", n_pulse=" << mls.NPulse(block_idx) << endl;
        if(0==max_DP_error)
        {
            // index pulses with at leat one echo
            vector<bool> has_echo(mls.NPulse(block_idx),false);
            vector<int> v_pulse_with_echo_idx(mls.NPulse(block_idx), -1);
            int idx=0;
            for(XPulseIndex pulse_idx=0; pulse_idx<mls.NPulse(block_idx); pulse_idx++)
            {
                has_echo[pulse_idx] = (mls.NbOfEcho(block_idx, pulse_idx)>0);
                if(has_echo[pulse_idx]) v_pulse_with_echo_idx[pulse_idx] = idx++;
            }
            // create triangles list TODO remove large triangles
            vector<Triangle> v_tri;
            for(XPulseIndex pulse_idx=0; pulse_idx<mls.NPulse(block_idx)-PPL-1; pulse_idx++)
            {
                // build triangles if the 3 pulses echoed
                if(mls.NbOfEcho(block_idx, pulse_idx)>0 && mls.NbOfEcho(block_idx, XPulseIndex(pulse_idx+PPL+1))>0)
                {
                    if(mls.NbOfEcho(block_idx, XPulseIndex(pulse_idx+PPL))>0)
                    {
                        if(MaxEdgeSize(mls, block_idx, pulse_idx, pulse_idx+PPL+1, pulse_idx+PPL) < tri_threshold)
                        {
                            v_tri.push_back(
                                        Triangle(v_pulse_with_echo_idx[pulse_idx],
                                                 v_pulse_with_echo_idx[pulse_idx+PPL+1],
                                        v_pulse_with_echo_idx[pulse_idx+PPL]));
                        }
                    }
                    if(mls.NbOfEcho(block_idx, XPulseIndex(pulse_idx+1))>0)
                    {
                        if(MaxEdgeSize(mls, block_idx, pulse_idx, pulse_idx+1, pulse_idx+PPL+1) < tri_threshold)
                        {
                            v_tri.push_back(
                                        Triangle(v_pulse_with_echo_idx[pulse_idx],
                                                 v_pulse_with_echo_idx[pulse_idx+1],
                                        v_pulse_with_echo_idx[pulse_idx+PPL+1]));
                        }
                    }
                }
            }
            if("off" == format) WriteOff(mls, ept_folder, block_idx, has_echo, v_tri, oss.str().c_str(), idx, pivot_E, pivot_N);
            else WritePly(mls, ept_folder, block_idx, has_echo, v_tri, oss.str().c_str(), idx, pivot_E, pivot_N);
        }
        else // Douglas-Peucker
        {
            vector<DP_point> v_dppts;
            // start of the interval for DP algo
            unsigned int idx_start=0;
            float theta_prev=mls.Theta(block_idx, (XPulseIndex)0);
            float half_pi = 1.57079632679;
            for(XPulseIndex pulse_idx=0; pulse_idx<mls.NPulse(block_idx); pulse_idx++)
            {
                if(mls.NbOfEcho(block_idx, pulse_idx)>0)
                {
                    DP_point pt;
                    pt.pulse_idx = pulse_idx;
                    pt.echo_idx = mls.IdxLastEcho(block_idx, pulse_idx);
                    pt.P = mls.P(block_idx, pt.echo_idx);
                    pt.keep = false; // init
                    pt.dec_idx = -1;
                    pt.theta = mls.Theta(block_idx, pulse_idx);
                    v_dppts.push_back(pt);

                    if(fmod(pt.theta, half_pi) < fmod(theta_prev, half_pi)) // run DP every pi/2 lidar rotation
                    {
                        //cout << "Quadrant: " << v_dppts[idx_start].theta << "," << pt.theta << endl;
                        //cout << "indices: " << idx_start << "," << v_dppts.size()-1 << endl;
                        DPDecimate(v_dppts, max_DP_error, idx_start, v_dppts.size()-1);
                        idx_start=v_dppts.size()-1;
                    }
                    theta_prev = pt.theta;
                }
            }
            // create decimated points vector
            vector<DP_point> v_dec;
            for(unsigned int i=0; i<v_dppts.size(); i++) if(v_dppts[i].keep) v_dec.push_back(v_dppts[i]);
            cout << "keeping: " << v_dec.size() << "/" << v_dppts.size() <<
                    " points=%" << 100.*(float)v_dec.size()/(float)v_dppts.size() << endl;
            vector<DP_point>().swap(v_dppts); // clear original points
            // triangulate
            unsigned int i1=0, i2=0; // (decimated) indices of points in 2 successive lines
            while(v_dec[i2].pulse_idx < mls.PulsePerLine()) i2++;
            i2--; // index of first point of second line
            vector<Triangle> v_tri;
            do
            {
                float delta = v_dec[i2].theta - v_dec[i1].theta;
                while(delta > M_PI) delta-=2*M_PI;
                while(delta < -M_PI) delta+=2*M_PI;
                if(delta>0) {v_tri.push_back(Triangle(i1, i1+1, i2)); i1++;}
                else {v_tri.push_back(Triangle(i1, i2+1, i2)); i2++;}
                if(i1 > i2) cout << "ERR " << i1 << ">" << i2 << endl;
                // else cout << v_dec[i1].theta << " " << v_dec[i2].theta << endl;
            } while(i2<v_dec.size());
            if("off" == format) WriteOff(mls, ept_folder, block_idx, v_dec, v_tri, oss.str().c_str(), pivot_E, pivot_N);
            else WritePly(mls, ept_folder, block_idx, v_dec, v_tri, oss.str().c_str(), pivot_E, pivot_N);
        }
        mls.Free(block_idx);
    }
    cout << "Done in " << (double)(clock() - start) / (double) CLOCKS_PER_SEC << "s" << endl;
    return 0;
}

