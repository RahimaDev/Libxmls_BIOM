
#define WHAT "SelfRegistration: refines a trajectory to enhance registration on overlaps"

#include <ctime>
#include <iostream>
#include <limits>
#include "libXMls/XMls.h"
#include "libXBase/XPt2D.h"
#include "LiteGeom/LgPolygon2.hpp"
#include "LiteGeom/LgIntersection2.hpp"
#include "LiteGeom/LgDistance2.hpp"

using namespace std;

class Footprint:public Lg::Polygon2
{
public:
    int id;
    double time;
    Lg::Point2 center;
    vector<Footprint*> v_overlap;
};

//-----------------------------------------------------------------------------
int main(int argc, char **argv)
{
    cout << WHAT << endl;
    if(argc < 2)
    {
        cout << "Usage: " << argv[0] << "  sbet" << endl;
        //cout << "Usage: " << argv[0] << "  sbet laser_calib.xml ept_folder output_sbet [band_width=10]" << endl;
        cout << "sbet: path to the sbet (trajecto) file" << endl;
        //cout << "laser_calib.xml: calibration file for the laser" << endl;
        //cout << "ept_folder: folder containing the echo pulse tables (generated with EptExport)" << endl;
        //cout << "tri_threshold: threshold on max triangle edge size to add the triangle (default=0.5m)" << endl;
        //cout << "max_DP_error: maximum Douglas-Peucker error (default=0m=do not decimate)" << endl;
        //cout << "format: choose between: ply (default), off" << endl;
        //cout << "(start|end)_time: start and end time of the laser points to export (default=everything)" << endl;
        return 0;
    }

    int i_arg=1;
    // required
    string sbet(argv[i_arg++]);
    //string laser_calib(argv[i_arg++]);
    //string ept_folder(argv[i_arg++]);
    //string output_sbet(argv[i_arg++]);

    // optional
    float band_width=10.f, band_width2 = band_width*band_width;
    if(i_arg < argc) band_width = atof(argv[i_arg++]);
    int min_delta_block=10;

    clock_t start = clock();
    // constructor and infos accessible after construction
    cout << "Reading " << sbet << endl;
    XTrajecto traj(sbet);
    traj.Load();
    traj.PreComputeGeoref();
    cout << traj.Nevent() << " traj events in " << (double)(clock() - start) / (double) CLOCKS_PER_SEC << "s" << endl;
    //XMls mls(ept_folder, laser_calib, &traj);
    vector<Footprint> v_footprint;
    int step = 100;
    for(unsigned int i=0; i<traj.Nevent()-step; i+=step)
    {
        XArchiGeoref georef0=traj.GetGeoref(i), georef1=traj.GetGeoref(i+1);
        Footprint fp;
        XPt3D XP0=georef0.Translation(), XP1=georef1.Translation();
        XPt3D Xv0=georef0.Rotation().A, Xv1=georef1.Rotation().A;
        Lg::Point2 P0(XP0.X, XP0.Y), P1(XP1.X, XP1.Y), v0(Xv0.X, Xv0.Y), v1(Xv1.X, Xv1.Y);
        fp.push_back(P0-0.5*band_width*v0);
        fp.push_back(P0+0.5*band_width*v0);
        fp.push_back(P1+0.5*band_width*v1);
        fp.push_back(P1-0.5*band_width*v1);
        fp.center = 0.5*(P0+P1);
        fp.id = i;
        v_footprint.push_back(fp);
    }
    cout << v_footprint.size() << " footprints in " << (double)(clock() - start) / (double) CLOCKS_PER_SEC << "s" << endl;
    int i_inter=0;
    for(unsigned int i=0; i<v_footprint.size(); i++)
    {
        if(i%100 == 0) cout << '.' << flush;
        for(unsigned int j=i+min_delta_block; j<v_footprint.size(); j++)
        {
            //cout << v_footprint[i] << " inter " << v_footprint[j] << endl;
            //if(Lg::Intersects(v_footprint[i],v_footprint[j]))
            if(Lg::SquaredDistance(v_footprint[i].center, v_footprint[j].center) < band_width2)
            {
                v_footprint[i].v_overlap.push_back(&(v_footprint[j]));
                v_footprint[j].v_overlap.push_back(&(v_footprint[i]));
                i_inter++;
            }
        }
    }
    int min_overlap=v_footprint.size(), max_overlap=0, n_overlap=0, i_max=0;
    for(unsigned int i=0; i<v_footprint.size(); i++)
    {
        int overlap = v_footprint[i].v_overlap.size();
        n_overlap += overlap;
        if(overlap<min_overlap) min_overlap = overlap;
        if(overlap>max_overlap) {max_overlap = overlap; i_max=i;}
    }
    cout << endl << i_inter << " intersections found in " << (double)(clock() - start) / (double) CLOCKS_PER_SEC << "s" << endl;
    cout << "min " << min_overlap << " avg " << n_overlap/v_footprint.size() << " max " << max_overlap << endl;
    for(unsigned int i=0; i<v_footprint[i_max].v_overlap.size(); i++)
    {
        cout << v_footprint[i_max].v_overlap[i]->id << " ";
    }
    cout << endl;
    return 0;
}

