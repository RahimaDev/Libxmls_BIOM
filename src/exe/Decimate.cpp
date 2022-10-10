
#define WHAT "Decimate a mesh in .OFF format using CGAL"

#include <iostream>
#include <fstream>
#include <ctime>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
// Stop-condition policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Polyhedron_3<Kernel> Surface_mesh;
namespace SMS = CGAL::Surface_mesh_simplification ;

using namespace std;

//-----------------------------------------------------------------------------
int main(int argc, char **argv)
{
    cout << WHAT << endl;
    if(argc < 3)
    {
        cout << "Usage: " << argv[0] << "  input.off [output.off=output.off] [stop_ratio=0.5 BoundaryWeight=0.5 VolumeWeight=0.5 ShapeWeight=0]" << endl;
        cout << "tri_threshold: threshold on max triangle edge size to add the triangle (default=0.5m)" << endl;
        return 0;
    }

    int i_arg=1;
    // required
    string input(argv[i_arg++]);
    string output = "output.off";
    if(i_arg < argc) output = string(argv[i_arg++]);

    // optional
    float stop_ratio=0.5;
    if(i_arg < argc) stop_ratio = atof(argv[i_arg++]);
    SMS::LindstromTurk_params LT_params;
    if(i_arg < argc) LT_params.BoundaryWeight = atof(argv[i_arg++]);
    if(i_arg < argc) LT_params.VolumeWeight = atof(argv[i_arg++]);
    if(i_arg < argc) LT_params.ShapeWeight = atof(argv[i_arg++]);

    clock_t start = clock();
    Surface_mesh surface_mesh;

    std::ifstream is(input.c_str());
    std::ofstream os(output.c_str());
    string line;
    getline(is, line);
    if(line != "OFF") cout << "ERROR: Not an OFF file, starts with " << line << endl;
    os << "OFF" << endl;
    getline(is, line);
    while(line[0] == '#') {os << line << endl; getline(is, line);}
    is.seekg(0);
    is >> surface_mesh ;
    // This is a stop predicate (defines when the algorithm terminates).
    SMS::Count_ratio_stop_predicate<Surface_mesh> stop(stop_ratio);

    // This the actual call to the simplification algorithm.
    // The surface mesh and stop conditions are mandatory arguments.
    // The index maps are needed because the vertices and edges
    // of this surface mesh lack an "id()" field.
    int r = SMS::edge_collapse
              (surface_mesh
              ,stop
              ,CGAL::parameters::vertex_index_map(get(CGAL::vertex_external_index,surface_mesh)) // for CGAL <4.6, remove ::parameters
               .halfedge_index_map  (get(CGAL::halfedge_external_index  ,surface_mesh))
               .get_cost (SMS::LindstromTurk_cost<Surface_mesh>(LT_params))
               .get_placement(SMS::LindstromTurk_placement<Surface_mesh>(LT_params))
              );

    std::cout << "\nFinished...\n" << r << " edges removed.\n"
              << (surface_mesh.size_of_halfedges()/2) << " final edges.\n" ;

    os << '#'; // we copied the header to keep the comments
    os << surface_mesh ;
    cout << "Done in " << (double)(clock() - start) / (double) CLOCKS_PER_SEC << "s" << endl;
    return 0;
}


