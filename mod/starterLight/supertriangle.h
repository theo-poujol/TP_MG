#ifndef SUPERTRIANGLE_H
#define SUPERTRIANGLE_H

#include "triangle.h"

class SuperTriangle : public Triangle {
public :
    SuperTriangle(Point A, Point B, Point C) : Triangle(A*2.0f, B*2.0f, C*2.0f) {}
};

#endif // SUPERTRIANGLE_H
