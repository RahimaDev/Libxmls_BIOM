#pragma once

/// edges with priorities = energy decrease associated to switching its state

#include "common.h"
#include "PriorityVertex.h"

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
                 PriorityVertex * p_target,
                 bool on=false);
    float priority() const;
    void add(float priority);
    void add(){add(priority());}
    void remove();
    inline bool isValid() const {return m_queue_it!=pedge_queue::iterator();}
    float length() const {return m_length;}
    void updatePriority();
    void switchState();
};

std::ostream& operator<<(std::ostream& os, const PriorityEdge& e);
