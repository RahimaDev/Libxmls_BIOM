
#include "PriorityEdge.h"

PriorityEdge::PriorityEdge(pedge_queue * p_queue,
                           PriorityVertex * p_source,
                           PriorityVertex * p_target,
                           bool on):
    mp_queue(p_queue), m_queue_it(pedge_queue::iterator()),
    mp_source(p_source), mp_target(p_target), m_on(on)
{
    float dx=mp_target->x-mp_source->x;
    float dy=mp_target->y-mp_source->y;
    m_length=sqrt(dx*dx+dy*dy);
    p_source->mvp_edge.push_back(this);
    p_target->mvp_edge.push_back(this);
}

float PriorityEdge::priority() const
{
    int n_source = mp_source->onValence();
    int n_target = mp_target->onValence();
    if(m_on) return -length()+w(n_source-1)-w(n_source)+w(n_target-1)-w(n_target);
    return length()+w(n_source+1)-w(n_source)+w(n_target+1)-w(n_target);
}

void PriorityEdge::add(float priority)
{
    m_queue_it = mp_queue->insert(std::make_pair(priority, this));
}

void PriorityEdge::remove()
{
    if(isValid()) mp_queue->erase(m_queue_it);
    else std::cout << "Removing a removed pedge" << std::endl;
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

std::ostream& operator<<(std::ostream& os, const PriorityEdge& e)
{
    os << "len=" << e.length() << ",pri=" << e.priority() << '|' << (e.m_on?"on":"off");
    return os;
}
