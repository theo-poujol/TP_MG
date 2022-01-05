#ifndef EDGE_H
#define EDGE_H

#include "point.h"
class Edge
{
public:
    Edge(Point &p1, Point &p2) : p1(&p1), p2(&p2) {}
    Point * p1, *p2;

    bool operator == (const Edge &e) const;

};

#endif // EDGE_H
