
#define WHAT "RieglFlexport: flexible Riegl exporter from echo/pulse tables (ept)"

#include <ctime>
#include <iostream>
#include <fstream>
#include <set>
#include "libXMls/XMls.h"
//#include "libXBase/XArchiXMLException.h"

#define COUT_SUB 1000000

using namespace std;

int Typesize(string type_name)
{
    if(type_name == "float32" || type_name == "int32" || type_name == "uint32") return 4;
    if(type_name == "float64" || type_name == "int64" || type_name == "uint64") return 8;
    if(type_name == "int8" || type_name == "uint8") return 1;
    if(type_name == "int16" || type_name == "uint16") return 2;
    return 0;
}

struct attrib_meta_info_t
{
    string name, type;
    unsigned int bytesize;
    bool required;
    attrib_meta_info_t(string name_, string type_):name(name_), type(type_), bytesize(Typesize(type_)), required(false){}
};

typedef vector<attrib_meta_info_t> v_attrib_info_t;

inline void AddAttribInfo(v_attrib_info_t & v_info, string attrib_name, string type_name)
{
    v_info.push_back(attrib_meta_info_t(attrib_name, type_name));
}

v_attrib_info_t ExportableAttributes()
{
    // order matters ! if you add new attribs, compute them in the same order
    v_attrib_info_t v_attrib_info;
    // preffered but non mandatory in decreasing type size order
    // if you add new attribs not at end update indices and boolean computation before main loop
    // 0
    AddAttribInfo(v_attrib_info, "GPS_time", "float64");
    AddAttribInfo(v_attrib_info, "range", "float32");
    AddAttribInfo(v_attrib_info, "theta", "float32");
    AddAttribInfo(v_attrib_info, "phi", "float32");

    // 4
    AddAttribInfo(v_attrib_info, "x_sensor", "float32");
    AddAttribInfo(v_attrib_info, "y_sensor", "float32");
    AddAttribInfo(v_attrib_info, "z_sensor", "float32");

    // 7
    AddAttribInfo(v_attrib_info, "x_origin", "float32");
    AddAttribInfo(v_attrib_info, "y_origin", "float32");
    AddAttribInfo(v_attrib_info, "z_origin", "float32");

    // 10
    AddAttribInfo(v_attrib_info, "x_ins", "float32");
    AddAttribInfo(v_attrib_info, "y_ins", "float32");
    AddAttribInfo(v_attrib_info, "z_ins", "float32");

    // 13
    AddAttribInfo(v_attrib_info, "x", "float32");
    AddAttribInfo(v_attrib_info, "y", "float32");
    AddAttribInfo(v_attrib_info, "z", "float32");

    // 16
    AddAttribInfo(v_attrib_info, "xVelocity", "float32");
    AddAttribInfo(v_attrib_info, "yVelocity", "float32");
    AddAttribInfo(v_attrib_info, "zVelocity", "float32");

    // 19
    AddAttribInfo(v_attrib_info, "roll", "float32");
    AddAttribInfo(v_attrib_info, "pitch", "float32");
    AddAttribInfo(v_attrib_info, "plateformHeading", "float32");
    AddAttribInfo(v_attrib_info, "wanderAngle", "float32");

    // 23
    AddAttribInfo(v_attrib_info, "xAcceleration", "float32");
    AddAttribInfo(v_attrib_info, "yAcceleration", "float32");
    AddAttribInfo(v_attrib_info, "zAcceleration", "float32");

    // 26
    AddAttribInfo(v_attrib_info, "xBodyAngularRate", "float32");
    AddAttribInfo(v_attrib_info, "yBodyAngularRate", "float32");
    AddAttribInfo(v_attrib_info, "zBodyAngularRate", "float32");

    // 29
    AddAttribInfo(v_attrib_info, "northPositionRMSError", "float32");
    AddAttribInfo(v_attrib_info, "eastPositionRMSError", "float32");
    AddAttribInfo(v_attrib_info, "downPositionRMSError", "float32");

    // 32
    AddAttribInfo(v_attrib_info, "northVelocityRMSError", "float32");
    AddAttribInfo(v_attrib_info, "eastVelocityRMSError", "float32");
    AddAttribInfo(v_attrib_info, "downVelocityRMSError", "float32");

    // 35
    AddAttribInfo(v_attrib_info, "RollRMSError", "float32");
    AddAttribInfo(v_attrib_info, "PitchRMSError", "float32");
    AddAttribInfo(v_attrib_info, "headingRMSError", "float32");

    // 38
    AddAttribInfo(v_attrib_info, "amplitude", "float32");
    AddAttribInfo(v_attrib_info, "reflectance", "float32");

    // 40
    AddAttribInfo(v_attrib_info, "deviation", "uint8");
    AddAttribInfo(v_attrib_info, "nb_of_echo", "uint8");
    AddAttribInfo(v_attrib_info, "num_echo", "uint8");

    return v_attrib_info;
}

// not implemented as a map for 2 reasons:
// - We need to ensure the attributes ordering
// - efficiency in critical loop
inline bool Require(v_attrib_info_t & v_attrib_info, string name)
{
    if(name.empty()) return false;
    for(v_attrib_info_t::iterator it=v_attrib_info.begin(); it!=v_attrib_info.end(); it++)
    {
        if(it->name == name)
        {
            it->required = true;
            return true;
        }
    }
    if(name[0]!='(') cout << "Warning: Required unknown attribute " << name << endl;
    return false;
}

// return true if at least one attribute in the range is required
inline bool Requires(v_attrib_info_t & v_attrib_info, int i_start, int i_end)
{
    for(int i=i_start; i<i_end; i++) if(v_attrib_info[i].required) return true;
    return false;
}

template <typename T> void Write(char * & it, T data)
{
    *reinterpret_cast<T*>(it) = data;
    it += sizeof(T);
}

//-----------------------------------------------------------------------------
int main(int argc, char **argv)
{
    cout << WHAT << endl;
    v_attrib_info_t v_attrib_info = ExportableAttributes();
    if(argc < 7)
    {
        // TODO: xml to store additionally the coord system and list of attribs to store
        cout << "Usage: " << argv[0] << "  sbet acc laser_calib.xml ept_folder attrib.txt out_name.ply [(int)start_time (int)end_time pivot_E=0 pivot_N=0 z_shift=0]" << endl;
        cout << "sbet: path to the sbet (trajecto) file" << endl;
        cout << "acc: path to the accuracy file (sometimes called smrmsg)" << endl;
        cout << "laser_calib.xml: calibration file for the laser" << endl;
        cout << "ept_folder: folder containing the echo pulse tables (generated with EptExport)" << endl;
        cout << "out_name.ply: name of the ply file to write (cf attrib list)" << endl;
        cout << "(start|end)_time: start and end time of the laser points to export (default=everything)" << endl;
        cout << "pivot_(E|N): easting/northing of the pivot (integer, default: use first trajectory point, round at 100m)" << endl;
        cout << "attrib.txt: text file with the names of all attributes to export among:" << endl;
        for(v_attrib_info_t::iterator it=v_attrib_info.begin(); it!=v_attrib_info.end(); it++)
        {
            cout << "(" << it->type << ") " << it->name << endl;
        }
        cout << "The most simple is to copy paste the lines above in attrib.txt file and keep only the lines you need" << endl;
        return 0;
    }

    int i_arg=1;
    // required
    string sbet(argv[i_arg++]);
    string acc(argv[i_arg++]);
    string laser_calib(argv[i_arg++]);
    string ept_folder(argv[i_arg++]);
    string attrib_filename(argv[i_arg++]);
    string out_name(argv[i_arg++]);

    // optional
    XSecond i_start = -1, i_end = -1;
    if(i_arg < argc) i_start = atoi(argv[i_arg++]);
    if(i_arg < argc) i_end = atoi(argv[i_arg++]);
    int pivot_E = 0, pivot_N = 0; double z_shift=0.;
    if(i_arg < argc) pivot_E = atof(argv[i_arg++]);
    if(i_arg < argc) pivot_N = atof(argv[i_arg++]);
    if(i_arg < argc) z_shift = atof(argv[i_arg++]);

    // read attrib file
    ifstream attrib_file(attrib_filename);
    int n_attrib = 0;
    do
    {
        string attrib; attrib_file >> attrib;
        if(Require(v_attrib_info, attrib)) n_attrib++;
    } while(!attrib_file.eof());
    cout << attrib_filename << " selects " << n_attrib << " known attributes" << endl;

    // read all echo/pulse tables in time interval
    clock_t start = clock();
    XTrajecto traj(sbet, acc);
    XMls mls(ept_folder, laser_calib, &traj);
    cout << mls.NPulseAttrib() << "/" << mls.NEchoAttrib() << " echos/pulses attributes found in " << ept_folder << endl;
    mls.Select(i_start, i_end);

    if(pivot_E == 0 && pivot_N == 0)
    {
        XArchiGeoref G = traj.GetGeoref(0);
        pivot_E = 100*(int)(G.Translation().X/100);
        pivot_N = 100*(int)(G.Translation().Y/100);
    }

    // iteration on echos
    start = clock();
    ofstream fileOut(out_name.c_str());
    if(!fileOut.good())
    {
        cout << "Cannot open " + out_name + " for writing\n";
        return 3;
    }

    // write text header
    fileOut << "ply\nformat binary_little_endian 1.0" << endl;
    fileOut << "comment Generated with RieglFlexport" << endl;
    fileOut << "comment IGN offset Pos " << pivot_E << " " << pivot_N << " 0" << endl; // z not handled by lidarformat, apply it to coordinates
    fileOut << "comment IGN offset Time " << mls.m_ept_time_pivot.ToString() << endl;
    fileOut << "comment IGN EPT " << ept_folder << endl;
    fileOut << "comment IGN sbet " << sbet << endl;
    fileOut << "comment IGN acc " << acc << endl;
    fileOut << "comment IGN calib " << laser_calib << endl;
    fileOut << "element vertex " << mls.NTotalEcho() << endl;
    int echo_bytesize = 0;
    for(v_attrib_info_t::iterator it=v_attrib_info.begin(); it!=v_attrib_info.end(); it++) if(it->required)
    {
        fileOut << "property " << it->type << " " << it->name << endl;
        echo_bytesize += it->bytesize;
    }
    fileOut << "end_header" << endl;
    // precompute booleans to know what computation blocks are needed (for efficiency)
    bool requires_xyz_sensor = Requires(v_attrib_info, 4, 7);
    bool requires_xyz_origin = Requires(v_attrib_info, 7, 10);
    bool requires_xyz_ins = Requires(v_attrib_info, 10, 13);
    bool requires_xyz = Requires(v_attrib_info, 13, 16);
    bool requires_sbet = Requires(v_attrib_info, 16, 29);
    bool requires_accuracy = Requires(v_attrib_info, 29, 38);

    // gather the required attribs
    XFloatAttrib * p_ampl = NULL, * p_refl = NULL;
    XUCharAttrib * p_dev = NULL;
    if(v_attrib_info[38].required) p_ampl = mls.GetEchoAttrib<XFloatAttrib>("amplitude");
    if(v_attrib_info[39].required) p_refl = mls.GetEchoAttrib<XFloatAttrib>("reflectance");
    if(v_attrib_info[40].required) p_dev = mls.GetEchoAttrib<XUCharAttrib>("deviation");
    // TODO: authorize export of any attribute in ept_folder

    unsigned long buffer_size = echo_bytesize * mls.NTotalEcho();
    char * buffer = new char[buffer_size], * it = buffer;
    cout.precision(16);
    bool compare = false;
    double max_error = 0.;
    for(XBlockIndex block_idx=0; block_idx<mls.NBlock(); block_idx++)
    {
        mls.Load(block_idx);
        for(XEchoIndex i_echo=0; i_echo<mls.NEcho(block_idx); i_echo++)
        {
            bool display_echo = (i_echo == 0);
            // time
            double time = mls.Time(block_idx, i_echo);
            if(v_attrib_info[0].required) Write<double>(it, time);

            // spherical coords
            if(v_attrib_info[1].required) Write<float>(it, mls.Range(block_idx, i_echo));
            if(v_attrib_info[2].required) Write<float>(it, mls.Theta(block_idx, i_echo));
            if(v_attrib_info[3].required) Write<float>(it, mls.Phi(block_idx, i_echo));
            if(display_echo)
            {
                cout << "["<< 100*block_idx/mls.NBlock() << "%] Time:" << time << "|sph:" <<
                        mls.Range(block_idx, i_echo) << " " << mls.Theta(block_idx, i_echo) << " " << mls.Phi(block_idx, i_echo) << endl;
            }

            // euclidian coords
            if(requires_xyz_sensor)
            {
                XPt3D Pt_sensor = mls.P(block_idx, i_echo);
                if(v_attrib_info[4].required) Write<float>(it, Pt_sensor.X);
                if(v_attrib_info[5].required) Write<float>(it, Pt_sensor.Y);
                if(v_attrib_info[6].required) Write<float>(it, Pt_sensor.Z + z_shift);
                if(display_echo) cout << "Pt_sensor: " << Pt_sensor << endl;
            }
            if(requires_xyz_origin)
            {
                XPt3D Pt_origin = mls.Cworld(block_idx, i_echo);
                if(v_attrib_info[7].required) Write<float>(it, Pt_origin.X - pivot_E);
                if(v_attrib_info[8].required) Write<float>(it, Pt_origin.Y - pivot_N);
                if(v_attrib_info[9].required) Write<float>(it, Pt_origin.Z + z_shift);
                if(display_echo) cout << "Pt_origin: " << Pt_origin << endl;
            }
            if(requires_xyz_ins)
            {
                XPt3D Pt_ins = mls.Oworld(block_idx, i_echo);
                if(v_attrib_info[10].required) Write<float>(it, Pt_ins.X - pivot_E);
                if(v_attrib_info[11].required) Write<float>(it, Pt_ins.Y - pivot_N);
                if(v_attrib_info[12].required) Write<float>(it, Pt_ins.Z + z_shift);
                if(display_echo) cout << "Pt_ins:" << Pt_ins << endl;
            }
            if(requires_xyz)
            {
                XPt3D Pt_ground = mls.Pworld(block_idx, i_echo);
                if(v_attrib_info[13].required) Write<float>(it, Pt_ground.X - pivot_E);
                if(v_attrib_info[14].required) Write<float>(it, Pt_ground.Y - pivot_N);
                if(v_attrib_info[15].required) Write<float>(it, Pt_ground.Z + z_shift);
                if(display_echo) cout << "Pt_ground: " << Pt_ground << endl;
                if(compare)
                {
                    XPt3D Pt_ground_frame = mls.Pworld_interpol_angles(block_idx, i_echo);
                    XPt3D err = Pt_ground - Pt_ground_frame;
                    double err_norm = sqrt(err.X*err.X + err.Y*err.Y + err.Z*err.Z);
                    if(err_norm > max_error) max_error = err_norm;
                    if(display_echo)
                    {
                        cout << "Pt_ground_frame: " << Pt_ground_frame << endl;
                        cout << "Error: " << err_norm << endl;
                    }
                }
            }
            // rest of the sbet
            if(requires_sbet)
            {
                SbetEvent sbet_event = mls.Sbet(block_idx, i_echo);
                if(v_attrib_info[16].required) Write<float>(it, sbet_event.m_xVelocity);
                if(v_attrib_info[17].required) Write<float>(it, sbet_event.m_yVelocity);
                if(v_attrib_info[18].required) Write<float>(it, sbet_event.m_zVelocity);
                if(v_attrib_info[19].required) Write<float>(it, sbet_event.m_roll);
                if(v_attrib_info[20].required) Write<float>(it, sbet_event.m_pitch);
                if(v_attrib_info[21].required) Write<float>(it, sbet_event.m_plateformHeading);
                if(v_attrib_info[22].required) Write<float>(it, sbet_event.m_wanderAngle);
                if(v_attrib_info[23].required) Write<float>(it, sbet_event.m_xAcceleration);
                if(v_attrib_info[24].required) Write<float>(it, sbet_event.m_yAcceleration);
                if(v_attrib_info[25].required) Write<float>(it, sbet_event.m_zAcceleration);
                if(v_attrib_info[26].required) Write<float>(it, sbet_event.m_xBodyAngularRate);
                if(v_attrib_info[27].required) Write<float>(it, sbet_event.m_yBodyAngularRate);
                if(v_attrib_info[28].required) Write<float>(it, sbet_event.m_zBodyAngularRate);
                if(display_echo) cout << "Sbet:" << sbet_event << endl;
            }
            // uncertainties
            if(requires_accuracy)
            {
                AccuracyEvent acc_event = mls.Accuracy(block_idx, i_echo);
                if(v_attrib_info[29].required) Write<float>(it, acc_event.m_northPositionRMSError);
                if(v_attrib_info[30].required) Write<float>(it, acc_event.m_eastPositionRMSError);
                if(v_attrib_info[31].required) Write<float>(it, acc_event.m_downPositionRMSError);
                if(v_attrib_info[32].required) Write<float>(it, acc_event.m_northVelocityRMSError);
                if(v_attrib_info[33].required) Write<float>(it, acc_event.m_eastVelocityRMSError);
                if(v_attrib_info[34].required) Write<float>(it, acc_event.m_downVelocityRMSError);
                if(v_attrib_info[35].required) Write<float>(it, acc_event.m_RollRMSError);
                if(v_attrib_info[36].required) Write<float>(it, acc_event.m_PitchRMSError);
                if(v_attrib_info[37].required) Write<float>(it, acc_event.m_headingRMSError);
                if(display_echo)
                {
                    cout << "Acc:" << acc_event << endl;
                }
            }
            // physical
            if(v_attrib_info[38].required) Write<float>(it, p_ampl->at(block_idx)[i_echo]);
            if(v_attrib_info[39].required) Write<float>(it, p_refl->at(block_idx)[i_echo]);
            if(v_attrib_info[40].required) Write<unsigned char>(it, p_dev->at(block_idx)[i_echo]);
            if(v_attrib_info[41].required) Write<unsigned char>(it, mls.NumEcho(block_idx, i_echo));
            if(v_attrib_info[42].required) Write<unsigned char>(it, mls.NbOfEcho(block_idx, i_echo));
            if(display_echo)
            {
                if(v_attrib_info[38].required) cout << "Ampl:" << p_ampl->at(block_idx)[i_echo];
                if(v_attrib_info[39].required) cout << "|Refl:" << p_refl->at(block_idx)[i_echo];
                if(v_attrib_info[40].required) cout << "|Dev:" << p_dev->at(block_idx)[i_echo];
                cout << "|" << mls.NumEcho(block_idx, i_echo) << "/" << mls.NbOfEcho(block_idx, i_echo) << endl;
            }
        }
        mls.Free(block_idx);
    }
    if(compare) cout << "max_error: " << max_error << endl;
    cout << "Writing " << mls.NTotalEcho() << " echos of size " << echo_bytesize << "=" << buffer_size/1000000 << "MB to " << out_name << endl;
    fileOut.write(buffer, buffer_size); // todo: split buffer if too big
    fileOut.close();
    delete buffer;
    cout << "Done in " << (double)(clock() - start) / (double) CLOCKS_PER_SEC << "s" << endl;
    return 0;
}

