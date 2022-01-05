#ifndef TRIANGLE_H
#define TRIANGLE_H

#include <iostream>
#include "math.h"
#include "point.h"
#include "edge.h"

class Triangle
{
public:
    Triangle(Point &A, Point &B, Point &C) : A(A), B(B), C(C) {}
    Point A;
    Point B;
    Point C;

    Edge * AB = new Edge(A, B);
    Edge * BC = new Edge(B, C);
    Edge * CA = new Edge(C, A);

    float a_norm = A.norm();
    float b_norm = B.norm();
    float c_norm = C.norm();

    float circum_x = (a_norm * (C.y - B.y) + b_norm * (A.y - C.y) + c_norm * (B.y - A.y)) / (A.x * (C.y - B.y) + B.x * (A.y - C.y) + C.x * (B.y - A.y));
    float circum_y = (a_norm * (C.x - B.x) + b_norm * (A.x - C.x) + c_norm * (B.x - A.x)) / (A.y * (C.x - B.x) + B.y * (A.x - C.x) + C.y * (B.x - A.x));

    // Origine du cercle circonscrit
    Point * circum = new Point(circum_x/2, circum_y/2);

    // Rayon du cercle circonscrit
    float circum_radius = A.dist(*circum);

    bool isInside(Point &p);

    bool operator == (const Triangle &t) const;

    bool operator < (const Triangle &t) const;
};

#endif // TRIANGLE_H
