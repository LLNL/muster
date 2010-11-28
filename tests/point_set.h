#ifndef POINT_SET_H
#define POINT_SET_H

#ifdef HAVE_CONFIG_H
#include "muster-config.h"
#endif // HAVE_CONFIG_H

#include <vector>
#include <string>
#include <iostream>
#include "point.h"

namespace cluster {

  class point_set {
  public:
    point_set();
    ~point_set();
    
    void add_point(point p);
    void add_point(double x, double y) { add_point(point(x,y)); }

    ///
    /// Range normalization
    ///
    void normalize();

    ///
    /// Parses a string containing points in parentheses, like this:
    ///  "(1, 1)  (2, 2) (3, 3)"
    /// Appends parsed points to points vector.
    ///
    void parse_points(const std::string& str);

    ///
    /// Parses a line containing a point with dimensions separated by a commas:
    ///  "x, y"
    /// Appends parsed point to points vector.
    ///
    void parse_point_csv(const std::string& str);
    
    ///
    /// Loads a csv file.
    ///
    void load_csv_file(std::istream& input);


    ///
    /// Write a csv file out, optionally including information about
    /// a clustering of the points in the file.
    ///
    void write_csv_file(std::ostream& out, partition *clustering=NULL);

    /// Get a reference to the actual points vector here.
    std::vector<point>& points() { return points_; }

    /// Allow indexing.
    point& operator[](size_t i) { return points_[i]; }

    /// Number of points in this point set.
    size_t size()  { return points_.size(); }

    double min_x() { return min_x_; }
    double max_x() { return max_x_; }
    double min_y() { return min_y_; }
    double max_y() { return max_y_; }
    
  private:
    std::vector<point> points_;

    double min_x_, max_x_;
    double min_y_, max_y_;

#ifdef MUSTER_HAVE_MPI
    friend void scatter(point_set& points, int root, MPI_Comm comm);
#endif
  }; // point_set
  

#ifdef MUSTER_HAVE_MPI
  void scatter(point_set& points, int root, MPI_Comm comm);
#endif

} // namespace cluster

#endif // POINT_SET_H
