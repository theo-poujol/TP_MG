#include "point.h"

#include <iostream>
#include "math.h"

bool Point::operator == (const Point &other) const {
    if (this->x == other.x && this->y == other.y) return true;
    return false;
}

Point& Point::operator * (float v) noexcept {
    float new_x = this->x * v;
    float new_y = this->y * v;
//        float new_z = this->z * v;
    Point * p = new Point(new_x, new_y);
    return *p;
}


bool Point::operator < (const Point &other) const {

    float r1 = sqrt( pow(this->x, 2) + pow(this->y, 2));
    float r2 = sqrt( pow(other.x, 2) + pow(other.y, 2));

    return (r1 < r2);
}


bool Point::operator > (const Point &other) const {

    float r1 = sqrt( pow(this->x, 2) + pow(this->y, 2));
    float r2 = sqrt( pow(other.x, 2) + pow(other.y, 2));

    return (r1 > r2);
}


float Point::norm() const {return x*x + y*y;}

float Point::dist(const Point &other) const {
    float dx = this->x - other.x;
    float dy = this->y - other.y;
    return pow(dx, 2) + pow(dy, 2);
}

void Point::setVh(MyMesh::VertexHandle &vh) {
    this->vh = &vh;
}
