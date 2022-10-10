#pragma once

/// Export a .ply file from a generic class giving access to individual triangles

#include "common.h"
#include "PriorityEdge.h"

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
    virtual char * Buffer();
};

class ASTriSet:public AbstractTriSet
{
public:
    Alpha_shape_2 * mp_as;
    unsigned int m_n_tri;
    ASTriSet(Alpha_shape_2 * p_as);
    virtual unsigned int Ntri() {return m_n_tri;}
    virtual char * Buffer();
};

class PedgeTriSet:public AbstractTriSet
{
public:
    std::vector<PriorityEdge*> mvp_pedge;
    unsigned int m_n_tri;
    PedgeTriSet(std::vector<PriorityEdge*> & vp_pedge);
    virtual unsigned int Ntri() {return m_n_tri;}
    virtual char * Buffer();
};

void WritePly(AbstractTriSet * p_triset, param_t param);
