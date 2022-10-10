
#define WHAT "Create a 2D Delaunay triangulation of a slice of Lidar data in local coords in ept"

#include "PlyExport.h"
#include "libXMls/XEchoPulseTables.h"
//#include "libXMls/write_geojson.hpp"

using namespace std;

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

double Ew(vector<PriorityVertex*> & v_pvertex)
{
	double ret=0.;
	for(auto & pvertex_it:v_pvertex) ret += w(pvertex_it->onValence());
	return ret;
}

double Elen(vector<PriorityEdge*> v_pedge)
{
	double ret=0.;
	for(auto & pedge_it:v_pedge) if(pedge_it->m_on) ret += pedge_it->length();
	return ret;
}
    
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

    XPulseIndex n_pulse = ept.PulsePerLine();
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
    vector<PriorityVertex*> v_pvertex;
    vector<PriorityEdge*> v_pedge;
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
        for(;vit!=A.finite_vertices_end();vit++) //if(A.classify(*vit) != Alpha_shape_2::EXTERIOR)
        {
            v_pvertex.push_back(new PriorityVertex(vit->point().x(),vit->point().y()));
            vertex_map.insert(make_pair( vit, v_pvertex.back() ));
        }
        AsEdgeIt eit = A.finite_edges_begin();
        for(;eit!=A.finite_edges_end();eit++) if(A.classify(*eit) != Alpha_shape_2::EXTERIOR)
        {
            PriorityVertex * p_source=NULL, * p_target=NULL;
            int i = eit->second; // idx of the vertex opposite to the edge in the triangle
            auto source_it = vertex_map.find(eit->first->vertex((i+1)%3));
            if (source_it != vertex_map.end()) p_source=source_it->second;
            auto target_it = vertex_map.find(eit->first->vertex((i+2)%3));
            if (target_it != vertex_map.end()) p_target=target_it->second;
            v_pedge.push_back(new PriorityEdge(&pqueue, p_source, p_target));
        }
        // add all edges to the priority queue.
        // needs to be done after structure is built because priority computation might depend on it
        for(auto & pedge_it:v_pedge) pedge_it->add();
    } // destroy temporaries (Alpha_shape_2 and v_pt)
	cout << v_pvertex.size() << " vertices " << v_pedge.size() << " edges" << endl;
	
    // greedy gradient descent adding
    double E = v_pvertex.size()*w(0); // all edges off init=only isolated points
    unsigned int iter = 0;
    while(pqueue.begin()->first < 0. && iter++ < 2*v_pedge.size())
    {
        PriorityEdge & best_pedge = *pqueue.begin()->second;
        E += pqueue.begin()->first;
        cout << "it " << iter << (best_pedge.m_on?"-":"+") <<
                " len=" << best_pedge.length() << " pri=" <<
                pqueue.begin()->first << " E=" << E << endl;
        best_pedge.switchState();
    }
    // check energy
    double Elen_gp = Elen(v_pedge), Ew_gp = Ew(v_pvertex), E_gp=Elen_gp+Ew_gp;
    cout << "Elen_gp=" << Elen_gp << " Ew_gp=" << Ew_gp << " E=" << E_gp << " D=" << E_gp-E << endl;
    
    // greedy gradient descent removing
    for(auto & pedge_it:v_pedge) if(!pedge_it->m_on) pedge_it->switchState(); // all edges ON
    E = Elen(v_pedge) + Ew(v_pvertex);
    iter = 0;
    while(pqueue.begin()->first < 0. && iter++ < 2*v_pedge.size())
    {
        PriorityEdge & best_pedge = *pqueue.begin()->second;
        E += pqueue.begin()->first;
        cout << "it " << iter << (best_pedge.m_on?"-":"+") <<
                " len=" << best_pedge.length() << " pri=" <<
                pqueue.begin()->first << " E=" << E << endl;
        best_pedge.switchState();
    }
    // check energy
    double Elen_gm = Elen(v_pedge), Ew_gm = Ew(v_pvertex), E_gm=Elen_gm+Ew_gm;
    cout << "Elen_gp=" << Elen_gm << " Ew_gp=" << Ew_gm << " E=" << E_gm << " D=" << E_gm-E << endl;

    //ASTriSet triset(&A);
    PedgeTriSet triset(v_pedge);
    WritePly(&triset, param);

    //vector<bool> v_b(dt.number_of_edges()); //should we keep edge b
    cout << "Done in " << (double)(clock() - start) / (double) CLOCKS_PER_SEC << "s" << endl;
    return 0;
}


