#ifndef POINT_H
#define POINT_H


/**
 * Simple 1 dimensional point class for testing medoids algorithms.
 */
class Point {
public:
    const unsigned x;  //x position
    const unsigned y;  //y position
    
    /** 
     * Constructs a new point with position (x,y).
     */
    Point(unsigned _x, unsigned _y);

    /** Distance between this point and another. */
    double distance(const Point *other) const;
    
};


#endif //POINT_H
