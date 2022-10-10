#pragma once

/// vertices of priority edges

#include "common.h"

// forward declaration
class PriorityEdge;

// corresponding vertices
class PriorityVertex
{
public:
    float x,y;
    std::vector<PriorityEdge*> mvp_edge;
    PriorityVertex(float x_, float y_):x(x_),y(y_),mvp_edge(0){}
    int valence() const	{return mvp_edge.size();}
    int onValence() const;
    void update();
};

std::ostream& operator<<(std::ostream& os, const PriorityVertex& v);
