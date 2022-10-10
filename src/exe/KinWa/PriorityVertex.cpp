
#include "PriorityVertex.h"
#include "PriorityEdge.h"

int PriorityVertex::onValence() const
{
    int n=0;
    for(auto & p_edge:mvp_edge) if(p_edge->m_on) n++;
    return n;
}

void PriorityVertex::update()
{
    for(auto & p_edge:mvp_edge) p_edge->updatePriority();
}

std::ostream& operator<<(std::ostream& os, const PriorityVertex& v)
{
    os << v.x << ',' << v.y << '-' << v.onValence() << '/' << v.valence();
    return os;
}
