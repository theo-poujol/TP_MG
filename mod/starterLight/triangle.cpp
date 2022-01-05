#include "triangle.h"

// Fonction qui détermine si le point p est à l'intérieur du cercle circonscrit du triangle
// en calculant la distance entre le centre du cercle et ce point
// Si cette distance est <= au rayon du cercle alors p est dans le cercle
bool Triangle::isInside(Point &p) {
    float dist = p.dist(*this->circum);

    return abs(dist) <= circum_radius;
}

bool Triangle::operator == (const Triangle &t) const {
    return (this->A == t.A && this->B == t.B && this->C == t.C);
}

bool Triangle::operator < (const Triangle &t) const {
    return (this->A < t.A && this->B < t.B && this->C < t.C);
}

