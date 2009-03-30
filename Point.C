
#include "Point.h"

#include <cmath>
#include <iostream>
using namespace std;

Point::Point(unsigned _x, unsigned _y) : 
    x(_x),
    y(_y)
{ 
    ; //nothing to do.
}


/** Distance between this point and another. */
double Point::distance(const Point *other) const {
    double dx = (double)other->x - (double)x;
    double dy = (double)other->y - (double)y;
    double temp = dx * dx + dy * dy;
    return sqrt(temp);
}
