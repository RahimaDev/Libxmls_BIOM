
#define WHAT "Create a 2D Delaunay triangulation of a slice of Lidar data in local coords in ept"

#include "libXMls/XEchoPulseTables.h"
//#include "libXMls/write_geojson.hpp"
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
//typedef K::Iso_rectangle_2 Iso_rectangle_2;
//typedef K::Segment_2 Segment_2;
//typedef K::Ray_2 Ray_2;
//typedef K::Line_2 Line_2;

using namespace std;

// manifold weighting
float w(float n)
{
    if(n<1.) return 0.6f-0.2f*n; //isolated points
    if(n<2.) return 0.8f-0.4f*n; // antennas
    return .5f*(n-2.f); // manifold->non-manifold
}
float dw(float n)
{
    if(n<1.) return -0.2f;
    if(n<2.) return -0.4f;
    return 0.1f;
}

struct param_t
{
    string ept_dir, output_ply;
    int block_id;
};

template <typename T> void Write(char * & it, T data)
{
    *reinterpret_cast<T*>(it) = data;
    it += sizeof(T);
}

template <> void Write<Point>(char * & it, Point P)
{
    Write<float>(it, P.x());
    Write<float>(it, P.y());
    Write<float>(it, 0.f);
}

class AbstractTriSet
{
public:
    AbstractTriSet(){}
    virtual unsigned int Ntri()=0;
    virtual char * Buffer()=0;
    inline unsigned int VertexBufferSize() {return 9*sizeof(float)*Ntri();}
};

class DTTriSet:public AbstractTriSet
{
public:
    DT2 * mp_dt;
    DTTriSet(DT2 * p_dt):mp_dt(p_dt){}
    virtual unsigned int Ntri() {return mp_dt->number_of_faces();}
    virtual char * Buffer()
    {
        char * buffer = new char[VertexBufferSize()], * it = buffer;
        DT2::Finite_faces_iterator fit = mp_dt->finite_faces_begin();
        for(;fit!=mp_dt->finite_faces_end();fit++)
        {
            for(int i=0; i<3; i++) Write<Point>(it, fit->vertex(i)->point());
        }
        return buffer;
    }
};

class ASTriSet:public AbstractTriSet
{
public:
    Alpha_shape_2 * mp_as;
    unsigned int m_n_tri;
    ASTriSet(Alpha_shape_2 * p_as):mp_as(p_as),m_n_tri(0)
    {
        // we want to export triangles and singular edges
        Alpha_shape_2::Finite_faces_iterator fit = mp_as->finite_faces_begin();
        for(;fit!=mp_as->finite_faces_end();fit++)
            if(mp_as->classify(fit) == Alpha_shape_2::INTERIOR) m_n_tri++;
        Alpha_shape_2::Finite_edges_iterator eit = mp_as->finite_edges_begin();
        for(;eit!=mp_as->finite_edges_end();eit++)
            if(mp_as->classify(*eit) == Alpha_shape_2::SINGULAR) m_n_tri++;
    }
    virtual unsigned int Ntri() {return m_n_tri;}
    virtual char * Buffer()
    {
        char * buffer = new char[VertexBufferSize()], * it = buffer;
        Alpha_shape_2::Finite_faces_iterator fit = mp_as->finite_faces_begin();
        for(;fit!=mp_as->finite_faces_end();fit++)
            if(mp_as->classify(fit) == Alpha_shape_2::INTERIOR)
                for(int i=0; i<3; i++) Write<Point>(it, fit->vertex(i)->point());
        Alpha_shape_2::Finite_edges_iterator eit = mp_as->finite_edges_begin();
        for(;eit!=mp_as->finite_edges_end();eit++)
            if(mp_as->classify(*eit) == Alpha_shape_2::SINGULAR)
            {
                Segment seg = mp_as->segment(*eit);
                Point S = seg.source(), T = seg.target();
                for(int i=0; i<2; i++) Write<Point>(it, S);
                Write<Point>(it, T);
            }
        return buffer;
    }
};

void WritePly(AbstractTriSet * p_triset, param_t param)
{
    ofstream fileOut(param.output_ply);
    if(!fileOut.good())
    {
        cout << "Cannot open " + param.output_ply + " for writing\n";
        return;
    }

    // write text header
    unsigned int n_face = p_triset->Ntri(), n_vertex = 3*n_face;
    fileOut << "ply\nformat binary_little_endian 1.0" << endl;
    fileOut << "comment Generated from " << param.ept_dir << " block " << param.block_id << endl; // TODO
    fileOut << "element vertex " << n_vertex << endl;
    fileOut << "property float x" << endl;
    fileOut << "property float y" << endl;
    fileOut << "property float z" << endl;
    fileOut << "element face " << n_face << endl;
    fileOut << "property list uchar int vertex_indices" << endl;
    fileOut << "end_header" << endl;
    // vertex list
    char * buffer = p_triset->Buffer();
    unsigned int vertex_buffer_size = p_triset->VertexBufferSize();
    cout << "Writing " << n_vertex << " vertices of size " <<
            vertex_buffer_size << "=" << 1.e-6*vertex_buffer_size << "MB" << endl;
    fileOut.write(buffer, vertex_buffer_size);
    delete buffer;
    // triangle list
    unsigned int tri_bytesize = sizeof(unsigned char)+3*sizeof(int);
    unsigned long tri_buffer_size = tri_bytesize*n_face;
    buffer = new char[tri_buffer_size];
    char * it = buffer;
    for(unsigned int tri_idx=0; tri_idx<n_face; tri_idx++)
    {
        Write<unsigned char>(it, 3);
        for(int i=0; i<3; i++) Write<int>(it, 3*tri_idx+i);
    }
    cout << "Writing " << n_face << " triangles of size " <<
            tri_buffer_size << "=" << 1.e-6*tri_buffer_size << "MB" << endl;
    cout << "Total " << 1.e-6*(p_triset->VertexBufferSize()+tri_buffer_size) << "MB" << endl;
    fileOut.write(buffer, tri_buffer_size);
    fileOut.close();
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

//// Priority Edges ////
typedef Alpha_shape_2::Vertex AsVertex;
typedef Alpha_shape_2::Vertex_handle AsVertexHandle;
typedef Alpha_shape_2::Finite_vertices_iterator AsVertexIt;
typedef Alpha_shape_2::Edge AsEdge;
typedef Alpha_shape_2::Finite_edges_iterator AsEdgeIt;

// fonctor for edges comparison
struct AsEdgeCmp {
    bool operator()(const AsEdge& lhs, const AsEdge& rhs) const {
        // iterators do not support comparison but pointers do
        if(lhs.second == rhs.second) return lhs.first < rhs.first;
        return lhs.second < rhs.second;
    }
};

// fonctor for vertices comparison
struct AsVertexCmp {
    bool operator()(const AsVertex& lhs, const AsVertex& rhs) const {
        if(lhs.point().x() == rhs.point().x())
            return lhs.point().y() < rhs.point().y();
        return lhs.point().x() < rhs.point().x();
    }
};

// class to deal with edges in a priority queue with dynamic priority
class PriorityEdge;
// corresponding vertices
class PriorityVertex
{
public:
    float x,y;
    vector<PriorityEdge*> mvp_edge;
    PriorityVertex(float x_, float y_):x(x_),y(y_){}
    int valence() const	{return mvp_edge.size();}
    int onValence() const
    {
        int n=0;
        for(auto & p_edge:mvp_edge) if(p_edge->m_on) n++;
        return n;
    }
    void update()
    {
        for(auto & p_edge:mvp_edge) p_edge->updatePriority();
    }
};
template <> void Write<PriorityVertex>(char * & it, PriorityVertex P)
{
    Write<float>(it, P.x);
    Write<float>(it, P.y);
    Write<float>(it, 0.f);
}

typedef std::multimap<float, PriorityEdge*> pedge_queue;

class PriorityEdge
{
    pedge_queue * mp_queue;
    pedge_queue::iterator m_queue_it;
    float m_length;

public:
    //int m_idx;
    PriorityVertex * mp_source, * mp_target;
    bool m_on; // is the edge kept or not
    PriorityEdge(pedge_queue * p_queue,
                 PriorityVertex * p_source,
                 PriorityVertex * p_target, on=false)
    mp_queue(p_queue), m_queue_it(pedge_queue::iterator()),
    mp_source(p_source), mp_target(p_target), m_on(on)
    {
        float dx=mp_target->x-mp_source->x;
        float dy=mp_target->y-mp_source->y;
        m_length=sqrt(dx*dx+dy*dy);
        p_source->mvp_edge.push_back(this);
        p_target->mvp_edge.push_back(this);
    }
    float priority();
    void add(float priority);
    void add(){add(priority())};
    void remove();
    inline bool isValid() const {return m_queue_it!=pedge_queue::iterator();}
    float length() const {return m_length;}
    void updatePriority();
    void switchState();
};

void PriorityEdge::priority()
{
    int n_source = mp_source->onValence();
    int n_target = mp_target->onValence();
    if(m_on) return -length()+w(n_source-1)-w(n_source)+w(n_target-1)-w(n_target);
    return length()+w(n_source+1)-w(n_source)+w(n_target+1)-w(n_target));
}

void PriorityEdge::add(float priority)
{
    m_queue_it = mp_queue->insert(std::make_pair(priority, this));
}

void PriorityEdge::remove()
{
    if(isValid()) mp_queue->erase(m_queue_it);
    else cout << "Removing a removed pedge" << endl;
    m_queue_it = pedge_queue::iterator();
}

void PriorityEdge::updatePriority()
{
    remove();
    add();
}

void PriorityEdge::switchState()
{
    m_on = !m_on;
    // TODO: slight opti = prevent this from being updated twice
    mp_source->update();
    mp_target->update();
}

class PedgeTriSet:public AbstractTriSet
{
public:
    vector<PriorityEdge> mv_pedge;
    unsigned int m_n_tri;
    PedgeTriSet(vector<PriorityEdge> & v_pedge):mv_pedge(v_pedge),m_n_tri(0)
    {
        for(auto & pedge_it:mv_pedge) if(pedge_it.m_on) m_n_tri++;
    }
    virtual unsigned int Ntri() {return m_n_tri;}
    virtual char * Buffer()
    {
        char * buffer = new char[VertexBufferSize()], * it = buffer;
        for(auto & pedge_it:mv_pedge) if(pedge_it.m_on)
        {
            for(int i=0; i<2; i++) Write<Point>(it, *pedge_it.mp_source);
            Write<Point>(it, *pedge_it.mp_target);
        }
        return buffer;
    }
};

//-----------------------------------------------------------------------------
int main(int argc, char **argv)
{
    cout << WHAT << endl;
    if(argc < 3)
    {
        cout << "Usage: " << argv[0] << "  ept_dir output.ply [block_id=-1]" << endl;
        cout << "block_id: identifier of the block on which to compute the DT (default=-1=first block)" << endl;
        return 0;
    }

    int i_arg=1;
    param_t param;
    // required
    param.ept_dir = string(argv[i_arg++]);
    param.output_ply = string(argv[i_arg++]);

    // optional
    param.block_id=-1;
    if(i_arg < argc) param.block_id = atoi(argv[i_arg++]);

    clock_t start = clock();
    XEchoPulseTable ept(argv[1]);
    cout << ept.NPulseAttrib() << "/" << ept.NEchoAttrib() << " pulse/echo attributes found" << endl;

    // Select all blocks = create blocks structure but do not load them, only required attributes are accessible
    ept.Select(param.block_id);
    ept.Load(0); // load the first (and only) block
    cout << ept.NBlock() << " blocks selected" << endl;
    cout << ept.NTotalPulse() << "/" << ept.NTotalEcho() << " pulses/echoes" << endl;

    // this is the simplest way to make an existing attribute accessible
    XFloatAttrib * p_refl_attrib = ept.GetEchoAttrib<XFloatAttrib>("reflectance");
    cout << "p_refl_attrib->m_ept_path=" << p_refl_attrib->EptPath() << endl;

    int n_pulse = ept.PulsePerLine();
    if(false)
    {
        vector< pair<Point,unsigned> > v_pt;
        for(XPulseIndex i_pulse=0; i_pulse<n_pulse; i_pulse++)
        {
            XPt3D P = ept.P(0, ept.IdxLastEcho(0,i_pulse));
            v_pt[i_pulse] = make_pair( Point(P.X,P.Y), i_pulse );
        }
        DT2 dt(v_pt.begin(), v_pt.end());
        DTTriSet * p_triset = new DTTriSet(&dt);
        WritePly(p_triset, param);
    }
    // build PriorityEdge/Vertex structures
    vector<PriorityVertex> v_pvertex;
    vector<PriorityEdge> v_pedge;
    pedge_queue pqueue;
    {
        vector<Point> v_pt(n_pulse);
        for(XPulseIndex i_pulse=0; i_pulse<n_pulse; i_pulse++)
        {
            XPt3D P = ept.P(0, ept.IdxLastEcho(0,i_pulse));
            v_pt[i_pulse] = Point(P.X,P.Y);
        }
        //    vector<Point> v_pt;
        //    v_pt.push_back(Point(0., 0.1));
        //    v_pt.push_back(Point(0.,-0.1));
        //    v_pt.push_back(Point(-0.2,0.));
        //    v_pt.push_back(Point(-0.11,0.));
        //    v_pt.push_back(Point( 0.11,0.));
        //    v_pt.push_back(Point( 0.2,0.));
        Alpha_shape_2 A(v_pt.begin(), v_pt.end(), FT(100),
                        Alpha_shape_2::GENERAL);
        cout << "Optimal alpha: " << *A.find_optimal_alpha(1)<<endl;
        map<AsVertexHandle, PriorityVertex*> vertex_map;
        AsVertexIt vit = A.finite_vertices_begin();
        for(;vit!=A.finite_vertices_end();vit++) if(A.classify(*vit) != Alpha_shape_2::EXTERIOR)
        {
            v_pvertex.push_back(PriorityVertex(vit->point().x(),vit->point().y()))
                    vertex_map.insert(make_pair( eit, &v_pvertex.back() ));
        }
        AsEdgeIt eit = A.finite_edges_begin();
        for(;eit!=A.finite_edges_end();eit++) if(A.classify(*eit) != Alpha_shape_2::EXTERIOR)
        {
            PriorityVertex * p_source=NULL, p_target=NULL;
            int i = eit->second; // idx of the vertex opposite to the edge in the triangle
            auto source_it = vertex_map.find(eit->first->vertex((i+1)%3));
            if (source_it != vertex_map.end()) p_source=source_it->second;
            auto target_it = vertex_map.find(eit->first->vertex((i+2)%3));
            if (target_it != vertex_map.end()) p_target=target_it->second;
            v_pedge.push_back(PriorityEdge(pqueue, p_source, p_target));
        }
        // add all edges to the priority queue.
        // needs to be done after structure is built because priority computation might depend on it
        for(auto & pedge_it:v_pedge) pedge_it.add();
    } // destroy temporaries (Alpha_shape_2 and v_pt)

    double energy = v_pvertex.size()*w(0); // all edges off init=only isolated points
    // greedy gradient descent
    int iter = 0;
    while(pqueue.begin()->first < 0. & iter++ < 2*v_pedge.size())
    {
        PriorityEdge & best_pedge = *pqueue.begin()->second;
        energy += pqueue.begin()->first;
        cout << "it " << iter << (best_pedge.m_on?" removing ":" adding ") << best_pedge.m_idx << " length " <<
                best_pedge.length() << " pri " << best_pedge.priority() <<
                " energy " << energy << endl;
        best_pedge.switchState();
    }
    //ASTriSet triset(&A);
    PedgeTriSet triset(v_pedge);
    WritePly(&triset, param);

    //vector<bool> v_b(dt.number_of_edges()); //should we keep edge b
    cout << "Done in " << (double)(clock() - start) / (double) CLOCKS_PER_SEC << "s" << endl;
    return 0;
}


