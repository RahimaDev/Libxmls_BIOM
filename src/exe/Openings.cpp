#define WHAT "Polygon: extracts planar polygons from a point cloud"

# define PI           3.14159265358979323846  /* pi */

#include <iostream>
#include <fstream>
#include <string>
#include <iterator>
#include <vector>
#include <list>
#include <sstream> // iss
#include <stdlib.h> // atof, atoi, rand(), abs()
#include <string> // c_str(), to_string()
#include <omp.h>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <math.h> // sqrt
#include <random>
#include <algorithm>
#include <stdio.h> //printf
//#include "LidarFormat/LidarDataContainer.h"
//#include "LidarFormat/LidarFile.h"
#include "libXMls/XMls.h"
#include "LiteGeom/LgPoint3.hpp"
#include "LiteGeom/LgPlane3.hpp"
#include "LiteGeom/LgDistance3.hpp"
#include "LiteGeom/LgPoint2.hpp"
#include "LiteGeom/LgTriangle3.hpp"
#include "LiteGeom/LgLine3.hpp"
#include "LiteGeom/LgIntersection3.hpp"

#include <pcl/common/geometry.h>



//*** PCL***************************
#include <pcl/point_cloud.h>
#include <pcl/PCLPointCloud2.h>

#include <pcl/point_types.h>
#include <pcl/io/ply_io.h>
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/registration/icp.h>
#include <pcl/io/ply_io.h>
#include <boost/make_shared.hpp>
#include <pcl/console/parse.h>
#include <pcl/common/transforms.h>
#include <pcl/visualization/pcl_visualizer.h>
//////////////////////////////////
#include <pcl/ModelCoefficients.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/model_types.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/common/centroid.h>
#include <pcl/common/transforms.h>
#include <pcl/common/common.h>
#include <pcl/visualization/cloud_viewer.h>

#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/correspondence.h>
#include <pcl/features/normal_3d_omp.h>
#include <pcl/features/shot_omp.h>
#include <pcl/features/board.h>
#include <pcl/filters/uniform_sampling.h>
#include <pcl/recognition/cg/hough_3d.h>
#include <pcl/recognition/cg/geometric_consistency.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/kdtree/impl/kdtree_flann.hpp>
#include <pcl/common/transforms.h>
#include <pcl/console/parse.h>
/////////////// CGAL /////////////////////////
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include "RectangleFitting.hpp"
#include "PlaneFitting.hpp"
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Boolean_set_operations_2.h>

#include <CGAL/algorithm.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Alpha_shape_2.h>
#include <pcl/segmentation/extract_clusters.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <pcl/filters/extract_indices.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;


typedef K::FT FT;
typedef K::Point_2  Point;
typedef K::Segment_2  Segment;
typedef K::Vector_2 Vector;

typedef K::Triangle_2 Triangle;
typedef CGAL::Polygon_2<K>                           Polygon;
typedef CGAL::Polygon_with_holes_2<K>                Polygon_with_holes;

typedef CGAL::Alpha_shape_vertex_base_2<K> Vb;
typedef CGAL::Alpha_shape_face_base_2<K>  Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb> Tds;
typedef CGAL::Delaunay_triangulation_2<K,Tds> Triangulation_2;

typedef CGAL::Alpha_shape_2<Triangulation_2>  Alpha_shape_2;
typedef Alpha_shape_2::Alpha_shape_edges_iterator Alpha_shape_edges_iterator;
typedef Alpha_shape_2::Alpha_shape_vertices_iterator Alpha_shape_vertices_iterator;
typedef Alpha_shape_2::Face_handle Face_handle;
typedef Alpha_shape_2::Edge Edge;
typedef Alpha_shape_2::Vertex_handle Vertex_handle;
//vector<Lg::Point3f> A_array, B_array; // TODO: remove awful globals

using namespace std;
///////////////////////////

void saveResultAsOBJ(std::string filename ,pcl::PointCloud<pcl::PointXYZ>::Ptr cloud)
{
    // view_mutex_.lock();
    // view_reserve_mutex_.lock();



    // get filename
    //std::string filename = "new.obj";

    std::ofstream file;
    file.open(filename.c_str());

    size_t lineID = 0;
    size_t pointID = 1;
    std::map<size_t,size_t> lines2points;

    for(size_t i=0; i<cloud->points.size(); i++)
    {
        if((i%2)==0)
        {Eigen::Vector3d P1 =  Eigen::Vector3d(cloud->points[i].x,cloud->points[i].y,cloud->points[i].z);
            Eigen::Vector3d P2 =  Eigen::Vector3d(cloud->points[i+1].x,cloud->points[i+1].y,cloud->points[i+1].z);

            file << "v " << P1.x() << " " << P1.y() << " " << P1.z() << std::endl;
            file << "v " << P2.x() << " " << P2.y() << " " << P2.z() << std::endl;

            lines2points[lineID] = pointID;
            ++lineID;
            pointID+=2;

        }
        else
            continue;
    }

    std::map<size_t,size_t>::const_iterator it = lines2points.begin();
    for(; it!=lines2points.end(); ++it)
    {
        file << "l " << it->second << " " << it->second+1 << std::endl;
    }

    file.close();

    //view_reserve_mutex_.unlock();
    //view_mutex_.unlock();
}




////////////////////

int factoriel_recursif(int val)
{
    if(val == 0)
        return 1;
    else
        return val*factoriel_recursif(val-1);
}

struct Param
{
    //string sbet, acc;
    XMls::Param xmls_param;
    string pc_ply, poly_ply;
    XSecond block_id;
    int lidar_width, sample_window_width, growing_window_width;
    double error_threshold, alpha_shape_radius;
    uint min_plane_size, n_sample;
    bool region_growing, Contrario;
    double pivot_E, pivot_N, pivot_H;

    Param():pc_ply("PC.ply"), poly_ply("poly.ply"), block_id(0),
        lidar_width(1021),
        sample_window_width(20),
        growing_window_width(8),
        error_threshold(0.02),
        alpha_shape_radius(0.5),
        min_plane_size(1800),
        n_sample(500),
        region_growing(false),
        Contrario(false),
        pivot_E(1050958.), pivot_N(6842000.), pivot_H(50.) {}
};

class Plane : public Lg::Plane3d
{
public:
    unsigned char rc, gc, bc;
    vector<XEchoIndex> v_inlier_idx;
    Plane(){}
    Plane(Lg::Point3d A, Lg::Point3d B, Lg::Point3d C):
        Lg::Plane3d(A, B, C), rc(rand () % 255), gc(rand () % 255), bc(rand () % 255), v_inlier_idx(0)
    {
        //unsigned int test = rand();
        //std::cout << "COUCOU " << test <<std::endl;
        this->Normalize();
        //std::cout << "au revOIR " << test << std::endl;
    }
};
struct Line{
    pcl::PointXYZ p1;
    pcl::PointXYZ p2;
    char label;



};
struct Opening
{
    Line L1;
    Line L2;
    Line L3;
    Line L4;
};
struct DATA
{
    Plane pl;
    std::vector<Opening> op;
};
void writeOutLines( string filename, std::vector<DATA> data )
{
    string fileEdgePoints = filename + "lines.txt";
    FILE *fp2 = fopen( fileEdgePoints.c_str(), "w");
    int h=-1;
    for(int p=0; p<data.size(); p++)
    {
        Plane pl=data[p].pl;
        std::vector<Opening>OP;
        OP=data[p].op;

        fprintf( fp2, "%c  %.6lf   %.6lf  %.6lf  %.6lf  %d \n   ", 'p', pl.Normal().x(), pl.Normal().y(),pl.Normal().z(),pl.d(),p );

        for(int m=0; m<OP.size(); m++)
        {
            Opening oo=OP[m];
            fprintf(fp2,"%c  %d %d \n",'O',m,p);
            fprintf( fp2, "%c %.6lf   %.6lf   %.6lf %.6lf    %.6lf %.6lf %c %d %d \n",'L',  oo.L1.p1.x, oo.L1.p1.y, oo.L1.p1.z, oo.L1.p2.x, oo.L1.p2.y, oo.L1.p2.z, 'V',m,p );
            fprintf( fp2, "%c  %.6lf   %.6lf   %.6lf %.6lf    %.6lf %.6lf %c %d %d \n",'L', oo.L2.p1.x, oo.L2.p1.y, oo.L2.p1.z, oo.L2.p2.x, oo.L2.p2.y, oo.L2.p2.z, 'H',m,p );
            fprintf( fp2, "%c  %.6lf   %.6lf   %.6lf %.6lf    %.6lf %.6lf %c %d %d \n",'L', oo.L3.p1.x, oo.L3.p1.y, oo.L3.p1.z, oo.L3.p2.x, oo.L3.p2.y, oo.L3.p2.z, 'H',m,p);
            fprintf( fp2, "%c  %.6lf   %.6lf   %.6lf %.6lf    %.6lf %.6lf %c %d %d\n",'L', oo.L4.p1.x, oo.L4.p1.y, oo.L4.p1.z, oo.L4.p2.x, oo.L4.p2.y, oo.L4.p2.z, 'V',m,p );
        }}
    fclose( fp2 );
}
class XYZ
{
public:
    XMls * mp_mls;
    XDoubleAttrib * mp_x, *mp_y, *mp_z;
    XShortAttrib * mp_id;
    XUCharAttrib * mp_proc, * mp_r, * mp_g, * mp_b;

    XYZ(XMls * p_mls, Param param):mp_mls(p_mls)
    {
        mp_x = mp_mls->AddEchoAttrib<XDoubleAttrib>("x");
        mp_y = mp_mls->AddEchoAttrib<XDoubleAttrib>("y");
        mp_z = mp_mls->AddEchoAttrib<XDoubleAttrib>("z");
        mp_id = mp_mls->AddEchoAttrib<XShortAttrib>("id");
        mp_proc = mp_mls->AddEchoAttrib<XUCharAttrib>("processed");
        mp_r = mp_mls->AddEchoAttrib<XUCharAttrib>("r");
        mp_g = mp_mls->AddEchoAttrib<XUCharAttrib>("g");
        mp_b = mp_mls->AddEchoAttrib<XUCharAttrib>("b");

        for(XEchoIndex echo_idx=0; echo_idx<mp_mls->NEcho(0); echo_idx++)
        {
            XPt3D Pw = mp_mls->Pworld(0, echo_idx);
            //XPt3D Cw = mls.Cworld(0, echo_idx);
            mp_x->at(0)[echo_idx] = Pw.X-param.pivot_E;
            mp_y->at(0)[echo_idx] = Pw.Y-param.pivot_N;
            mp_z->at(0)[echo_idx] = Pw.Z-param.pivot_H;
            mp_id->at(0)[echo_idx] = -1;
            mp_proc->at(0)[echo_idx] = 0;
            //float e0=Cw.X-param.pivot_E, n0=Cw.Y-param.pivot_N, h0=Cw.Z-param.pivot_H;
        }
        cout << mp_mls->NEcho(0) << " x,y,z coords precomputed" << endl;
    }
    inline unsigned int N(){return mp_mls->NEcho(0);}
    inline Lg::Point3d P(XEchoIndex echo_idx)
    {
        return Lg::Point3d(mp_x->at(0)[echo_idx],mp_y->at(0)[echo_idx],mp_z->at(0)[echo_idx]);
    }
    inline Plane getPlane(XEchoIndex i1, XEchoIndex i2, XEchoIndex i3)
    {
        return Plane(P(i1), P(i2), P(i3));
    }
    inline bool processed(XEchoIndex echo_idx)
    {
        return mp_proc->at(0)[echo_idx]>=0;
    }
    inline bool hasPlane(XEchoIndex echo_idx)
    {
        return mp_id->at(0)[echo_idx]>=0;
    }
    inline short & id(XEchoIndex echo_idx)
    {
        return mp_id->at(0)[echo_idx];
    }
    inline void setRGB(XEchoIndex echo_idx, unsigned char r, unsigned char g, unsigned char b)
    {
        mp_r->at(0)[echo_idx]=r;
        mp_g->at(0)[echo_idx]=g;
        mp_b->at(0)[echo_idx]=b;
    }
    inline unsigned char r(XEchoIndex echo_idx){return mp_r->at(0)[echo_idx];}
    inline unsigned char g(XEchoIndex echo_idx){return mp_g->at(0)[echo_idx];}
    inline unsigned char b(XEchoIndex echo_idx){return mp_b->at(0)[echo_idx];}
};

bool contains(vector<int> & v_idx, int idx)
{
    for(int & i:v_idx) if(i==idx) return true;
    return false;
}

vector<int> getDiffRandIdx(int nIdx, int maxIdx) {
    vector<int> v_idx;
    for(int i=0; i<nIdx; i++)
    {
        int randIdx=0;
        do randIdx = rand() % maxIdx; while(contains(v_idx,randIdx));
        v_idx.push_back(randIdx);
    }
    return v_idx;
}

int RandInterval(int min, int max)
{
    return rand() % (max-min+1) + min;
}

/// gather the indices of unprocessed echoes in the sensor neighborhood of seed_idx
vector<XEchoIndex> getSensorNeighborhood(XYZ & xyz, XEchoIndex seed_idx, int window_width)
{
    vector<XEchoIndex> v_ret;
    XPulseIndex seed_pulse_idx = xyz.mp_mls->IdxPulse(0,seed_idx);
    unsigned int n_pulse = xyz.mp_mls->NPulse(0);
    int ppl = xyz.mp_mls->PulsePerLine();
    for(int dl=-window_width; dl<=window_width; dl++)
    {
        for(int dc=-window_width; dc<=window_width; dc++)
        {
            if(dc != 0 || dl != 0)
            {
                XPulseIndex pulse_idx = seed_pulse_idx+dl*ppl+dc;
                if(pulse_idx>=0 && pulse_idx<n_pulse)
                {
                    XEchoIndex echo_idx = xyz.mp_mls->IdxFirstEcho(0, pulse_idx);
                    for(int i=0; i < xyz.mp_mls->NbOfEcho(0, pulse_idx); i++)
                    {
                        if(!xyz.hasPlane(echo_idx))
                            v_ret.push_back(echo_idx);
                        echo_idx++;
                    }
                }
            }
        }
    }
    return v_ret;
}


/// return a vector orthogonal to n by choosing the longest of n^x and n^y
/// returns (0,0,0) if n is (0,0,0)
Lg::Point3d best_m(Lg::Point3d n)
{
    if(abs(n.x())>abs(n.y()))
        return Lg::Point3d(-n.z(),0,n.x());
    return Lg::Point3d(0,n.z(),-n.y());
}

template <class OutputIterator>
void
alpha_edges( const Alpha_shape_2&  A,
             OutputIterator out)
{
    for(Alpha_shape_edges_iterator it =  A.alpha_shape_edges_begin();
        it != A.alpha_shape_edges_end();
        ++it){
        *out++ = A.segment(*it);
    }
}

std::list<Polygon_with_holes> alpha_shape( Plane const & pl, XYZ & xyz, double radius)
{

    /// choose a local frame for the plane
    Lg::Point3d n = pl.Normal();
    Lg::Point3d m = best_m(n);
    m.Normalize();
    Lg::Point3d k = n^m;
    k.Normalize();

    vector<K::Point_2> v_pt2(pl.v_inlier_idx.size());
    Lg::Point3d O = pl.Point();
    int i=0;
    for(auto & inlier_idx:pl.v_inlier_idx)
    {
        Lg::Point3d V = xyz.P(inlier_idx)-O;
        v_pt2[i++] = Point
                (V*m, V*k);
    }
    Alpha_shape_2 A(v_pt2.begin(), v_pt2.end(),
                    FT(0.07),
                    Alpha_shape_2::REGULARIZED);
    std::cout << "##### Join triangles #####" << std::endl;

    std::list<Polygon> triangles;
    for(typename Alpha_shape_2::Finite_faces_iterator fit = A.finite_faces_begin();
        fit != A.finite_faces_end();
        ++fit){

        if(A.classify(fit) == Alpha_shape_2::INTERIOR){
            Triangle triangle = A.triangle(fit);
            Polygon poly;
            poly.push_back(triangle.vertex(0));
            poly.push_back(triangle.vertex(1));
            poly.push_back(triangle.vertex(2));
            triangles.push_back(poly);
        }
    }

    std::list<Polygon_with_holes> res;
    CGAL::join(triangles.begin(), triangles.end(), std::back_inserter (res));

    for(typename std::list<Polygon_with_holes>::iterator it_ring = res.begin();
        it_ring != res.end(); it_ring++){
        Polygon outer = it_ring->outer_boundary();
        // std::cout << "=== Outer ===" << std::endl;
        for(typename Polygon::Vertex_const_iterator it_vertex = outer.vertices_begin();
            it_vertex != outer.vertices_end(); it_vertex++){
            //std::cout << it_vertex->x() << "\t" << it_vertex->y() << std::endl;
        }
        for(typename Polygon_with_holes::Hole_const_iterator it_hole = it_ring->holes_begin();
            it_hole != it_ring->holes_end(); it_hole++)
        {
            Polygon inner = *it_hole;
            //  std::cout << "=== Inner ===" << std::endl;
            for(typename Polygon::Vertex_const_iterator it_vertex = inner.vertices_begin();
                it_vertex != inner.vertices_end(); it_vertex++){
                // std::cout << it_vertex->x() << "\t" << it_vertex->y() << std::endl;
            }
        }
    }
    return res;
}


///////////////////////////
double angle (Lg::Point3f coefficients1, Lg::Point3f  coefficients2)
{

    double a=coefficients1.x();
    double b=coefficients1.y();
    double c=coefficients1.z();


    double A=coefficients2.x();
    double B=coefficients2.y();
    double C=coefficients2.z();






    //calculates the angle between the two planes

    //F=normal1.normal2
    double F=(a * A)+(b* B)+(c * C);

    //H=||normal1||
    double H=sqrt((a * a)+(b* b)+(c * c));

    // L=||normal2||
    double L=sqrt((A * A)+(B*B)+(C * C ));

    double alpha=std::acos(F/(H*L));
    double theta=(alpha*180)/PI;
    return theta;

}

bool inside(Polygon_with_holes P, Plane pl,Lg::Point3d pt)
{
    bool Inside=false;
    Lg::Point3d n = pl.Normal();
    Lg::Point3d m = best_m(n);
    m.Normalize();
    Lg::Point3d k = n^m;
    k.Normalize();


    Lg::Point3d O = pl.Point();
    Lg::Point3d V = pt-O;
    Point A=Point(V*m, V*k);
    Polygon outer =P.outer_boundary();
    bool In_outer = (outer.bounded_side(A) == CGAL::ON_BOUNDED_SIDE);
    if (In_outer ==true)
        Inside=true;
    for(typename Polygon_with_holes::Hole_const_iterator it_hole = P.holes_begin();
        it_hole != P.holes_end(); it_hole++)
    {
        Polygon inner = *it_hole;
        bool In_inner = (inner.bounded_side(A) == CGAL::ON_BOUNDED_SIDE);
        if(In_inner==true)
            Inside=true;}

    return Inside;
}
bool ray_cross_poly(Point s,XPt3D p,Plane pl, XPt3D Origine)
{Lg::Point3d n = pl.Normal();
    Lg::Point3d m = best_m(n);
    m.Normalize();
    Lg::Point3d k = n^m;
    k.Normalize();


    Lg::Point3d O = pl.Point();

    Lg::Point3d s1 = O + CGAL::to_double(s.x())*m + CGAL::to_double(s.y())*k;




    Lg::Point3d vec1, vec2;
    double val1=0;
    double val2=0;
    vec1.x()=p.X-s1.x();
    vec1.y()=p.Y-s1.y();
    vec1.z()=p.Z-s1.z();

    vec2.x()=Origine.X-s1.x();
    vec2.y()=Origine.Y-s1.y();
    vec2.z()=Origine.Z-s1.z();
    val1=(vec1.x()*pl.N().x())+(vec1.y()*pl.N().y())+(vec1.z()*pl.N().z());


    val2=(vec2.x()*pl.N().x())+(vec2.y()*pl.N().y())+(vec2.z()*pl.N().z());



    if(((val1>0.07)&&(val2<0))||((val1<-0.07)&&(val2>0)))

        return true;
    else
        return false;
}

bool valid_int_point(Polygon_with_holes P, Plane pl,Lg::Point3d pt,XPt3D p, XPt3D Origine)
{bool valid=false;
    if(inside(P, pl,pt)==true)
    {Polygon outer =P.outer_boundary();
        if(ray_cross_poly(outer[0], p, pl, Origine))
            valid=true;
        else
        {
            for(typename Polygon_with_holes::Hole_const_iterator it_hole = P.holes_begin();
                it_hole != P.holes_end(); it_hole++)
            {
                Polygon inner = *it_hole;
                if(ray_cross_poly(inner[0], p, pl, Origine))
                    valid=true;

            }}}
    return valid;
}
////////////////////////Contrario function //////////////
static float logcombi(int k, int n)
{
    if (k>=n || k<=0)
    {return(0.0);


    }
    if (n-k<k) k=n-k;
    double r = 0.0;
    for (int i = 1; i <= k; i++)
        r += log10((double)(n-i+1))-log10((double)i);

    return static_cast<float>(r);
}

/// tabulate logcombi(.,n)
static void makelogcombi_n(int n, std::vector<float> & l)
{
    l.resize(n+1);
    for (int k = 0; k <= n; k++)
        l[k] = logcombi(k,n);
}

/// tabulate logcombi(k,.)
static void makelogcombi_k(int k,int nmax, std::vector<float> & l)
{
    l.resize(nmax+1);
    for (int n = 0; n <= nmax; n++)
        l[n] = logcombi(k,n);
} 
double p_sigma(double sigma, double d, double V)
//{ return ((2*sigma*d)/V);
{return ((2*sigma*d*d)/V);

}
double NFA(int N,int k,double p, int n_test)
{ double nfa=0.0;

    //omp_set_num_threads(8);

    int j;
    //   std::cout<<"ddddd="<< omp_get_num_threads( )<<std::endl;

#pragma omp parallel for private(j)
    for( j=k; j<N; j++){

        /*{*/  double d=logcombi(j, N);
        double h=j*log(p)+(N-j)*log(1-p);
        double f=exp(h);
        //std::cout<<"F"<<f<< " "<<"h"<<h<<std::endl;
        nfa=nfa+(d*f);
    }
    //
    std::cout<<(n_test*nfa)<< " "<<"nfa"<<std::endl;
    return (n_test*nfa);
} 

////////////
bool linePlaneIntersection(Plane pl,XPt3D Cw, XPt3D Pw,Lg::Point3d &pt)
{
    Lg::Point3d V=Lg::Point3d(Pw.X-Cw.X,Pw.Y-Cw.Y,Pw.Z-Cw.Z);
    Lg::Point3d A=Lg::Point3d (Cw.X,Cw.Y,Cw.Z);
    V.Normalize();

    double vn = V*pl.Normal();
    if(vn == 0) return false; // parallel
    double x = pl.Dist() - A*pl.Normal();
    if(x<0) return false; // not in right direction
    pt = A+ x/vn * V;
    return true;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{   std::vector<Lg::Point3d>  Points1,Points2;
    std::vector<Lg::Point3d>Points,Points0;

    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud2(new pcl::PointCloud<pcl::PointXYZ>);

    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud1 (new pcl::PointCloud<pcl::PointXYZ>);
    pcl::PointCloud<pcl::PointXYZ>Cloud_C ;
    pcl::PointCloud<pcl::PointXYZ>::Ptr Cloud_w(new pcl::PointCloud<pcl::PointXYZ>);
    cout << WHAT << endl;
    Param param, param1;

    if(argc < 5)
    {
        cout << "Usage: " << argv[0] << "  sbet laser_calib.xml ept_folder block_id ";
        cout << "[pc_ply poly_ply sample_window_width growing_window_width error_threshold min_plane_size n_sample region_growing]" << endl;
        cout << "sbet: sbet file" << endl;
        cout << "laser_calib.xml: calibration file for the laser" << endl;
        cout << "ept_folder: folder containing the echo pulse tables (generated with EptExport)" << endl;
        cout << "block_id: id of the block to process (look at info.txt in ept_folder for valid range)" << endl;
        cout << "pc_ply: name of the output colored point cloud ply file (default="
             << param.pc_ply << ")" << endl;
        cout << "poly_ply: name of the output polygons ply file (default="
             << param.poly_ply << ")" << endl;
        cout << "sample_window_width: width of the RANSAC sampling window (default="
             << param.sample_window_width << ")" << endl;
        cout << "growing_window_width: width of the window used for region growing (default="
             << param.growing_window_width << ")" << endl;
        cout << "error_threshold: threshold separating plane in/outliers (default="
             << param.error_threshold << "m)" << endl;
        cout << "min_plane_size: minimum number of planes inliers (default="
             << param.min_plane_size << ")" << endl;
        cout << "n_sample: number of RANSAC iterations (default="
             << param.n_sample << ")" << endl;
        cout << "region_growing: perform a region growing instead of computing all point/plane distances (default="
             << param.region_growing << ")" << endl;
        return 0;
    }

    int i_arg=1;

    // required
    param.xmls_param.sbet = string(argv[i_arg++]);
    param.xmls_param.laser_calib = string(argv[i_arg++]);
    param.xmls_param.ept_folder = string(argv[i_arg++]);
    param.block_id = atof(argv[i_arg++]);
    param.Contrario=atoi(argv[i_arg++]);
    // optional
    if(i_arg < argc) param.pc_ply = string(argv[i_arg++]);
    if(i_arg < argc) param.poly_ply = string(argv[i_arg++]);
    if(i_arg < argc) param.sample_window_width = atof(argv[i_arg++]);
    if(i_arg < argc) param.growing_window_width = atof(argv[i_arg++]);
    if(i_arg < argc) param.error_threshold = atof(argv[i_arg++]);
    if(i_arg < argc) param.min_plane_size = atoi(argv[i_arg++]);
    if(i_arg < argc) param.n_sample = atoi(argv[i_arg++]);
    if(i_arg < argc) param.region_growing = atoi(argv[i_arg++]);

    int start_s = clock();
    srand(time(NULL));

    // constructor and infos accessible after construction
    XTrajecto traj(param.xmls_param.sbet);
    XMls mls(param.xmls_param.ept_folder, param.xmls_param.laser_calib, &traj);
    mls.Select(param.block_id);
    cout << mls.NBlock() << " block(s) selected" << endl;
    mls.Load(0); // load the first (and only) block
    XPt3D Pivot;
    Pivot.X=1050958.;
    Pivot.Y=6842000.;
    Pivot.Z=50.;
    ;
    XYZ xyz(&mls, param); // structure to precompute and access x,y,z coords in local frame


    for(XEchoIndex echo_idx=0; echo_idx<mls.NEcho(0); echo_idx++)
    {

        XPt3D Pw = mls.Pworld(0, echo_idx);
        float x= Pw.X-Pivot.X;
        float y= Pw.Y-Pivot.Y;
        float z=Pw.Z-Pivot.Z;

        pcl::PointXYZ pt(x,y,z);
        std::cout<<pt<<std::endl;
        cloud1->points.push_back(pt);
        
    }



    pcl::io::savePLYFile("STR.ply", *cloud1);

    double thre=0;
    /////////////////////////////////////////////// Contrario method for choosing the threshold
    if(param.Contrario)
    {
        pcl::PointXYZ min_pt,max_pt, centre;
        pcl::getMinMax3D(*cloud1, min_pt, max_pt);
        double volume=((max_pt.x - min_pt.x)*(max_pt.y - min_pt.y)*(max_pt.z - min_pt.z));
        double LD=sqrt(((max_pt.x - min_pt.x)*(max_pt.x - min_pt.x))+((max_pt.y - min_pt.y)*(max_pt.y - min_pt.y))+((max_pt.z - min_pt.z)*(max_pt.z - min_pt.z)));
        std::vector<double>vecth;
        std::vector<double>vecth_nfa;



        vecth.push_back(0.01);
        vecth.push_back(0.013);
        vecth.push_back(0.017);
        vecth.push_back(0.02);

        vecth.push_back(0.025);
        vecth.push_back(0.035);
        vecth.push_back(0.015);
        vecth.push_back(0.03);
        vecth.push_back(0.045);
        vecth.push_back(0.04);



        double minnfa=1.0;
        double maxnfa=-1;
        // vector<Plane> v_plane;
        vector<int> v_size;
        for(int l=0; l<vecth.size(); l++)

        {    vector<Plane> v_plane;

            do {
                cout << "PLANE " << v_plane.size() << endl;
                Plane best_plane;
                for (uint n = 0; n < param.n_sample; n++)
                {
                    if(n%100==0) cout << "Sample " << n << endl;
                    XEchoIndex seed_idx = 0;
                    vector<XEchoIndex> v_neigh;
                    do
                    {
                        do seed_idx = rand()%xyz.N(); while(xyz.hasPlane(seed_idx));
                        v_neigh = getSensorNeighborhood(xyz, seed_idx, param.sample_window_width);
                    } while(v_neigh.size()<5);
                    vector<int> v_idx = getDiffRandIdx(2, v_neigh.size());
                    Plane cur_plane = xyz.getPlane(seed_idx, v_neigh[v_idx[0]], v_neigh[v_idx[1]]);
                    vector<XEchoIndex> v_inlier_idx;
                    if(!param.region_growing) // plane distance only
                    {
                        for(XEchoIndex echo_idx=0; echo_idx<xyz.N(); echo_idx++)
                            if(!xyz.hasPlane(echo_idx) &&
                                    abs(Lg::SignedDistance(xyz.P(echo_idx), cur_plane)) < vecth[l]


                                    )
                                v_inlier_idx.push_back(echo_idx);
                    }


                    if (v_inlier_idx.size() > best_plane.v_inlier_idx.size())
                    {
                        best_plane = cur_plane;
                        best_plane.v_inlier_idx = v_inlier_idx;
                        cout << "Best region inlier: " << v_inlier_idx.size() << endl;
                    }
                }
                cout << "Best region total inlier: " << best_plane.v_inlier_idx.size() << endl;
                v_plane.push_back(best_plane);

                for (XEchoIndex & echo_idx:best_plane.v_inlier_idx)
                {
                    xyz.id(echo_idx) = v_plane.size();
                    xyz.setRGB(echo_idx, best_plane.rc, best_plane.gc, best_plane.bc);
                }
            }
            while (v_plane.back().v_inlier_idx.size() >param.min_plane_size); // boucle principale RANSAC
            v_size.push_back(v_plane.size());
            v_plane.clear();
        }
        for(int l=0; l<vecth.size(); l++)
        {
            double p=p_sigma(vecth[l],LD,volume);

            int n=xyz.N()/(param.min_plane_size);
            double nfa=NFA(n-v_size[l],v_size[l],p,500*v_size[l]);
            vecth_nfa.push_back(nfa);
        }

        double min=1000;
        int w;
        auto smallest = std::min_element(std::begin(vecth_nfa), std::end(vecth_nfa));

        int position=std::distance(std::begin(vecth_nfa), smallest);
        thre=vecth[position];
        std::cout<<"the threshold is"<<" "<<thre<<std::endl;
    }
    else
    {thre=param.error_threshold;
    }
    XYZ xyz1(&mls, param);
    vector<Plane> v_plane1;
    do {
        cout << "PLANE " << v_plane1.size() << endl;
        Plane best_plane;
        for (uint n = 0; n < param.n_sample; n++)
        {
            if(n%100==0) cout << "Sample " << n << endl;
            XEchoIndex seed_idx = 0;
            vector<XEchoIndex> v_neigh;
            do
            {
                do seed_idx = rand()%xyz1.N(); while(xyz1.hasPlane(seed_idx));
                v_neigh = getSensorNeighborhood(xyz1, seed_idx, param.sample_window_width);
            } while(v_neigh.size()<5);
            vector<int> v_idx = getDiffRandIdx(2, v_neigh.size());
            Plane cur_plane = xyz1.getPlane(seed_idx, v_neigh[v_idx[0]], v_neigh[v_idx[1]]);
            vector<XEchoIndex> v_inlier_idx;
            if(!param.region_growing) // plane distance only
            {
                for(XEchoIndex echo_idx=0; echo_idx<xyz1.N(); echo_idx++)
                    if(!xyz1.hasPlane(echo_idx) &&
                            abs(Lg::SignedDistance(xyz1.P(echo_idx), cur_plane)) <thre)
                        v_inlier_idx.push_back(echo_idx);
            }

            if (v_inlier_idx.size() > best_plane.v_inlier_idx.size())
            {
                best_plane = cur_plane;
                best_plane.v_inlier_idx = v_inlier_idx;
                cout << "Best region inlier: " << v_inlier_idx.size() << endl;
            }
        }
        cout << "Best region total inlier: " << best_plane.v_inlier_idx.size() << endl;
        v_plane1.push_back(best_plane);

        for (XEchoIndex & echo_idx:best_plane.v_inlier_idx)
        {
            xyz1.id(echo_idx) = v_plane1.size();
            xyz1.setRGB(echo_idx, best_plane.rc, best_plane.gc, best_plane.bc);
        }
    }
    while (v_plane1.back().v_inlier_idx.size() >param.min_plane_size); // boucle principale RANSAC


    vector<Plane> PLS, PLS1, pl1,pl2;
    Plane pl_xy=v_plane1[0];
    for(int i=0; i<v_plane1.size();i++){
        std::cout<<i<<std::endl;
        double alpha=angle(v_plane1[i].N(),pl_xy.N());
        std::cout<<alpha<<std::endl;
        if((alpha>=84)&&(alpha<=96))
        {PLS1.push_back(v_plane1[i]);}

    }


    XPt3D pt1=mls.Pworld(0, 0);;
    XPt3D pt2=mls.Pworld(0,mls.NEcho(0)-1);
    Lg::Point3d w;
    w.x()=pt1.X;
    w.y()=pt1.Y;
    w.z()=pt1.Z;
    Lg::Point3d z;
    z.x()=pt2.X;
    z.y()=pt2.Y;
    z.z()=pt2.Z;
    double length_traj=(z-w).Norm ();

    std::cout<<"length_traj="<<length_traj<<std::endl;
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr CC(new pcl::PointCloud<pcl::PointXYZRGB>);
    int J=-1;
    if(PLS1.size()>0)

    {std::vector<DATA>VEC;


        for(int i=0; i<PLS1.size(); i++)
        {
            DATA dat;
            std::vector<Lg::Point3d>  Points3;
            pcl::PointCloud<pcl::PointXYZ>::Ptr cl(new pcl::PointCloud<pcl::PointXYZ>);
            pcl::PointCloud<pcl::PointXYZ>::Ptr cll(new pcl::PointCloud<pcl::PointXYZ>);
            std::vector<Opening>VO;
            std::list<Polygon_with_holes> RES= alpha_shape( PLS1[i], xyz1,  0.05);
            int H=-1;
            for(typename std::list<Polygon_with_holes>::iterator it_ring = RES.begin();
                it_ring != RES.end(); it_ring++){
                Polygon outer = it_ring->outer_boundary();
                double Area=CGAL::to_double(outer.area());

                if((Area>=29))
                    //if((Area>=6))
                {





                    for(XEchoIndex echo_idx=0; echo_idx<mls.NEcho(0); echo_idx++)
                    { if(i!=J)
                        {J=i;
                            std::cout<<"le plane est="<<i<<std::endl;
                        }

                        XPt3D Cw = mls.Cworld(0, echo_idx);
                        XPt3D Pw = mls.Pworld(0, echo_idx);
                        Cw=Cw-Pivot;
                        Pw=Pw-Pivot;
                        Lg::Line3d ln(Cw.X,Cw.Y,Cw.Z,Pw.X,Pw.Y,Pw.Z);
                        Lg::Point3d PT;

                        if((linePlaneIntersection( PLS1[i], Cw,Pw, PT ))&&(valid_int_point(*it_ring, PLS1[i],PT,Pw ,Cw)))
                        {



                            Points3.push_back(PT);

                        }}}}
            cl->width    = Points3.size();
            cl->height   = 1;
            cl->is_dense = false;

            for(int f=0; f<Points3.size(); f++)
            {
                float x=Points3[f].x();//-Pivot.X;
                float y=Points3[f].y();//-Pivot.Y;
                float z=Points3[f].z();//-Pivot.Z;
                cl->points.push_back(pcl::PointXYZ(x,y,z));

            }


            if(cl->points.size()>0)
            {

                /////////////////////////////////////////////////////////////////////////////////
                pcl::search::KdTree<pcl::PointXYZ>::Ptr kdtree(new pcl::search::KdTree<pcl::PointXYZ>);
                kdtree->setInputCloud(cl);

                // Euclidean clustering object.
                pcl::EuclideanClusterExtraction<pcl::PointXYZ> clustering;
                // Set cluster tolerance to 2cm (small values may cause objects to be divided
                // in several clusters, whereas big values may join objects in a same cluster).
                //clustering.setClusterTolerance(0.5);
                clustering.setClusterTolerance(0.2);
                // Set the minimum and maximum number of points that a cluster can have.
                clustering.setMinClusterSize(100);

                clustering.setMaxClusterSize(250000);
                clustering.setSearchMethod(kdtree);
                clustering.setInputCloud(cl);
                std::vector<pcl::PointIndices> clusters;
                clustering.extract(clusters);

                // For every cluster...
                int currentClusterNum = 1;
                for (std::vector<pcl::PointIndices>::const_iterator g = clusters.begin(); g != clusters.end(); ++g)
                {
                    // ...add all its points to a new cloud...
                    pcl::PointCloud<pcl::PointXYZ>::Ptr cluster(new pcl::PointCloud<pcl::PointXYZ>);
                    pcl::PointCloud<pcl::PointXYZRGB>::Ptr CL(new pcl::PointCloud<pcl::PointXYZRGB>);
                    for (std::vector<int>::const_iterator point = g->indices.begin(); point != g->indices.end(); point++)
                        cluster->points.push_back(cl->points[*point]);
                    cluster->width = cluster->points.size();
                    cluster->height = 1;
                    cluster->is_dense = true;

                    // ...and save it to disk.
                    if (cluster->points.size() <= 0)
                        break;
                    std::cout << "Cluster " << currentClusterNum << " has " << cluster->points.size() << " points." << std::endl;
                    int32_t red = rand ()>> 16 & 0xFF;
                    int32_t green = rand () >> 8 & 0xFF;
                    int32_t blue = rand ()& 0xFF;

                    for(int k=0; k<cluster->points.size(); k++)
                    {pcl::PointXYZRGB pt;
                        pt.x=cluster->points[k].x;
                        pt.y=cluster->points[k].y;
                        pt.z=cluster->points[k].z;
                        pt.r=red;
                        pt.g=green;
                        pt.b=blue;
                        uint32_t rgb = (red << 16) | (green << 8) | (blue);
                        pt.rgb=rgb;
                        CL->points.push_back(pt);
                    }
                    *CC+=*CL;

                    Lg::Point3d N,O;
                    N=PLS1[i].Normal();


                    O=PLS1[i].Point();
                    Lg::Point3d m = best_m(N);
                    m.Normalize();
                    Lg::Point3d k = N^m;
                    k.Normalize();

                    std::vector<Eigen::Vector2d>vec;
                    std::cout<<"mama"<<std::endl;
                    for(int l=0; l<cluster->points.size(); l++)
                    {
                        Lg::Point3d pt=Lg::Point3d(cluster->points[l].x,cluster->points[l].y,cluster->points[l].z);






                        Lg::Point3d V = pt-O;
                        Lg::Point2d A=Lg::Point2d(V*m, V*k);
                        Eigen::Vector2d Poi=Eigen::Vector2d(A.x(),A.y());
                        vec.push_back(Poi);
                    }
                    auto RectangleFitting = std::make_shared<RectangleFittingProblem>();

                    RectangleFitting->setData(vec);

                    robest::RANSAC solver;
                    solver.solve(RectangleFitting);
                    Eigen::Vector2d v1,v2;
                    RectangleFitting->getResult(v1,v2);
                    Lg::Point2d p1=Lg::Point2d(v1.x(),v1.y());
                    //////////////////////////////////////////////
                    Lg::Point3d s1 = O + p1.x()*m + p1.y()*k;
                    pcl::PointXYZ P1=pcl::PointXYZ(s1.x(),s1.y(),s1.z());
                    Lg::Point2d p2=Lg::Point2d(v1.x(),v2.y());
                    Lg::Point3d s2 = O + p2.x()*m + p2.y()*k;
                    pcl::PointXYZ P2=pcl::PointXYZ(s2.x(),s2.y(),s2.z());

                    /////////////////

                    Lg::Point2d p4=Lg::Point2d(v2.x(),v2.y());
                    Lg::Point3d s4 = O + p4.x()*m + p4.y()*k;
                    pcl::PointXYZ P4=pcl::PointXYZ(s4.x(),s4.y(),s4.z());
                    /////////////////
                    Lg::Point2d p3=Lg::Point2d(v2.x(),v1.y());
                    Lg::Point3d s3 = O + p3.x()*m + p3.y()*k;
                    pcl::PointXYZ P3=pcl::PointXYZ(s3.x(),s3.y(),s3.z());
                    //////////////////////

                    if(( pcl::geometry::distance(P1,P3)>=0.6)&&(pcl::geometry::distance(P1,P2)>=0.6))
                    {
                        cloud2->points.push_back(P1);
                        cloud2->points.push_back(P2);
                        cloud2->points.push_back(P1);
                        cloud2->points.push_back(P3);
                        cloud2->points.push_back(P4);
                        cloud2->points.push_back(P3);
                        cloud2->points.push_back(P4);
                        cloud2->points.push_back(P2);

                        Line L1,L2,L3,L4;
                        Opening op;
                        L1.p1=P1;
                        L1.p2=P2;
                        L1.label='V';
                        op.L1=L1;
                        ////////////////////////////////

                        L4.p1=P3;
                        L4.p2=P4;
                        L4.label='V';
                        op.L4=L4;
                        //////////////////////////

                        L2.p1=P1;
                        L2.p2=P3;
                        L2.label='H';
                        op.L2=L2;
                        ////////////////////////////////

                        L3.p1=P2;
                        L3.p2=P4;
                        L3.label='H';
                        op.L3=L3;

                        VO.push_back(op);
                        *CC+=*CL;


                    }}}
            if(VO.size()>0)
            {
                dat.pl=PLS1[i];
                dat.op=VO;
                VEC.push_back(dat);
            }


        }





        pcl::io::savePLYFileASCII("CompC.ply",*CC);
        saveResultAsOBJ("Lines.obj",cloud2);
        writeOutLines("win",VEC);}
    return 0;
}



