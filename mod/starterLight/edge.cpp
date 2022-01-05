#include "edge.h"

bool Edge::operator==(const Edge &e) const {
    if ( (this->p1 == e.p1 && this->p2 == e.p2) || ((this->p1 == e.p2 && this->p2 == e.p1)) ) return true;
    return false;
}

