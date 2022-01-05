#ifndef POINT_H
#define POINT_H

#include "mainwindow.h"

class Point
{
public:
    Point(float x, float y) : x(x), y(y) {}
    float x;
    float y;
    VertexHandle * vh;

     bool operator == (const Point &other) const;
     Point& operator * (float v) noexcept;
     bool operator < (const Point &other) const;
     bool operator > (const Point &other) const;
     float norm() const;
     float dist(const Point &other) const;
     void setVh (MyMesh::VertexHandle &vh);

};

#endif // POINT_H
