#pragma once

/// common includes and typedefs for KinWa

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/algorithm.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/squared_distance_2.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <list>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned int, K> Vb;
typedef CGAL::Triangulation_data_structure_2<Vb> Tds;
typedef K::FT FT;
typedef K::Point_2  Point;
typedef K::Segment_2  Segment;
typedef CGAL::Delaunay_triangulation_2<K,Tds> DT2;

typedef CGAL::Alpha_shape_vertex_base_2<K> AVb;
typedef CGAL::Alpha_shape_face_base_2<K>  AFb;
typedef CGAL::Triangulation_data_structure_2<AVb, AFb> ATds;
typedef CGAL::Delaunay_triangulation_2<K,ATds> ADT2;
typedef CGAL::Alpha_shape_2<ADT2>  Alpha_shape_2;
typedef Alpha_shape_2::Alpha_shape_edges_iterator Alpha_shape_edges_iterator;

struct param_t
{
    std::string ept_dir, output_ply;
    int block_id;
};

// manifold weighting
inline float w(float n)
{
    if(n<1.) return 0.6f-0.2f*n; //isolated points
    if(n<2.) return 0.8f-0.4f*n; // antennas
    return .5f*(n-2.f); // manifold->non-manifold
}
inline float dw(float n)
{
    if(n<1.) return -0.2f;
    if(n<2.) return -0.4f;
    return 0.1f;
}
