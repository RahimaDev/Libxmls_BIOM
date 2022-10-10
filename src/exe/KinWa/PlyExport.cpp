
/// Export a .ply file from a generic class giving access to individual triangles

#include "PlyExport.h"
using namespace std;

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

char * DTTriSet::Buffer()
{
    char * buffer = new char[VertexBufferSize()], * it = buffer;
    DT2::Finite_faces_iterator fit = mp_dt->finite_faces_begin();
    for(;fit!=mp_dt->finite_faces_end();fit++)
    {
        for(int i=0; i<3; i++) Write<Point>(it, fit->vertex(i)->point());
    }
    return buffer;
}

ASTriSet::ASTriSet(Alpha_shape_2 * p_as):
    mp_as(p_as),m_n_tri(0)
{
    // we want to export triangles and singular edges
    Alpha_shape_2::Finite_faces_iterator fit = mp_as->finite_faces_begin();
    for(;fit!=mp_as->finite_faces_end();fit++)
        if(mp_as->classify(fit) == Alpha_shape_2::INTERIOR) m_n_tri++;
    Alpha_shape_2::Finite_edges_iterator eit = mp_as->finite_edges_begin();
    for(;eit!=mp_as->finite_edges_end();eit++)
        if(mp_as->classify(*eit) == Alpha_shape_2::SINGULAR) m_n_tri++;
}

char * ASTriSet::Buffer()
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

template <> void Write<PriorityVertex>(char * & it, PriorityVertex P)
{
    Write<float>(it, P.x);
    Write<float>(it, P.y);
    Write<float>(it, 0.f);
}


PedgeTriSet::PedgeTriSet(vector<PriorityEdge*> & vp_pedge):
    mvp_pedge(vp_pedge),m_n_tri(0)
{
    for(auto & pedge_it:mvp_pedge) if(pedge_it->m_on) m_n_tri++;
}

char * PedgeTriSet::Buffer()
{
    char * buffer = new char[VertexBufferSize()], * it = buffer;
    for(auto & pedge_it:mvp_pedge) if(pedge_it->m_on)
    {
        for(int i=0; i<2; i++) Write<PriorityVertex>(it, *pedge_it->mp_source);
        Write<PriorityVertex>(it, *pedge_it->mp_target);
    }
    return buffer;
}

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

// Do I need to keep that ?
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
