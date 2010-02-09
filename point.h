#ifndef CLUSTER_TEST_POINT_H
#define CLUSTER_TEST_POINT_H

#include "libra-config.h"

#ifdef LIBRA_HAVE_MPI
#include <mpi.h>
#endif // LIBRA_HAVE_MPI 

#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include "partition.h"

namespace cluster {
  
  /// Simple 2 dimensional point class for testing medoids algorithms.
  struct point {
  public:
    int x, y;
    
    /// New point with position (x,y)
    point(int x, int y);

    /// New point at (0,0)
    point();

    /// Copy constructor
    point(const point& other);

    // Distance between this point and another. 
    double distance(const point& other) const {
      int dx = other.x - x;
      int dy = other.y - y;
      return ::sqrt((dx*dx) + (dy*dy));
    }
  
    point& operator+=(const point& other) {
      x += other.x;  y += other.y;
      return *this;
    }

    point operator+(const point& other) const {
      point result = *this;
      result += other;
      return result;
    }
  
    point& operator*=(const point& other) {
      x *= other.x;  y *= other.y;
      return *this;
    }

    point operator*(const point& other) const {
      point result = *this;
      result *= other;
      return result;
    }
  
    point& operator=(const point& other) { 
      x = other.x;  y = other.y;
      return *this;
    }

    bool operator==(const point& other) { 
      return x == other.x && y == other.y;
    }

    bool operator!=(const point& other) { 
      return !(*this == other);
    }

    void operator/=(int divisor) {
      x = static_cast<int>(round(x / (double)divisor));
      y = static_cast<int>(round(y / (double)divisor));
    }

    point operator/(int divisor) const {
      point p(x,y);
      p /= divisor;
      return p;
    }

#ifdef LIBRA_HAVE_MPI
    /// Returns the size of a packed point
    int packed_size(MPI_Comm comm) const;
  
    /// Packs a point into an MPI packed buffer
    void pack(void *buf, int bufsize, int *position, MPI_Comm comm) const;

    /// Unpacks a point from an MPI packed buffer
    static point unpack(void *buf, int bufsize, int *position, MPI_Comm comm);
#endif // LIBRA_HAVE_MPI
  };

  std::ostream& operator<<(std::ostream& out, const point& p);

  ///
  /// Parses a string containing points in parentheses, like this:
  ///  "(1, 1)  (2, 2) (3, 3)"
  /// Appends parsed points to points vector.
  ///
  void parse_points(const std::string& str, std::vector<point>& points);

  ///
  /// Draws a set of points in ascii with console colors.  Colors are assigned based on
  /// the partition provided.  Indices in points vector should correspond to ids in the 
  /// partition.
  /// 
  void draw(std::string label, std::vector<point>& points, const cluster::partition& parts, 
            std::ostream& out = std::cout);


  /// Distance bt/w two points
  struct point_distance {
    double operator()(const point& left, const point& right) {
      return left.distance(right);
    }
  };

} // namespace cluster

#endif // CLUSTER_TEST_POINT_H

