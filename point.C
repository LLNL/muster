#include "point.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#ifdef HAVE_MPI
#include "mpi_utils.h"
#endif // HAVE_MPI
#include "color.h"

#include <iomanip>
#include <cstdlib>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/algorithm/string.hpp>

using boost::numeric::ublas::matrix;
using namespace boost;
using namespace std;

namespace cluster {

  point::point(double _x, double _y) : x(_x), y(_y) { }

  point::point() : x(0), y(0) { }

  point::point(const point& other) : x(other.x), y(other.y) { }


#ifdef HAVE_MPI
  /// Returns the size of a packed point
  int point::packed_size(MPI_Comm comm) const {
    return 2 * mpi_packed_size(1, MPI_DOUBLE, comm);
  }
  
  /// Packs a point into an MPI packed buffer
  void point::pack(void *buf, int bufsize, int *position, MPI_Comm comm) const {
    PMPI_Pack(const_cast<double*>(&x), 1, MPI_DOUBLE, buf, bufsize, position, comm);
    PMPI_Pack(const_cast<double*>(&y), 1, MPI_DOUBLE, buf, bufsize, position, comm);
  }

  /// Unpacks a point from an MPI packed buffer
  point point::unpack(void *buf, int bufsize, int *position, MPI_Comm comm) {
    point p;
    PMPI_Unpack(buf, bufsize, position, &p.x, 1, MPI_DOUBLE,  comm);
    PMPI_Unpack(buf, bufsize, position, &p.y, 1, MPI_DOUBLE,  comm);
    return p;
  }
#endif // HAVE_MPI


  ostream& operator<<(ostream& out, const point& p) {
    return out << "(" << setw(2) << p.x << "," << setw(2) << p.y << ")";
  }


  void parse_points(const string& str, vector<point>& points) {
    string trimmed = trim_copy_if(str, is_any_of("() "));

    vector<string> parts;
    split(parts, trimmed, is_any_of("(,"));

    vector<double> values(parts.size());
    for (size_t i=0; i < parts.size(); i++) {
      parts[i]  = trim_copy_if(parts[i], is_any_of("() ,"));
      values[i] = strtod(parts[i].c_str(), NULL);
    }

    for (size_t i=1; i < values.size(); i += 2) {
      points.push_back(point(values[i-1], values[i]));
    }
  }


  void draw(string label, vector<point>& points, const partition& parts, ostream& out) {
    static const char *colors[] = {
      Blue, Green, Cyan, Red, Purple, Brown, Light_Gray, 
      Light_Blue, Light_Green, Light_Cyan, Light_Red, Light_Purple, 
      White, Dark_Gray, Yellow
    };
    static const size_t num_colors = sizeof(colors) / sizeof(const char *);

    double max_x=0;
    double max_y=0;

    for (size_t i=0; i < points.size(); i++) {
      max_x = max(points[i].x, max_x);
      max_y = max(points[i].y, max_y);
    }

    matrix<int> pmat((size_t)max_x+1, (size_t)max_y+1);
    for (size_t i=0; i < pmat.size1(); i++) {
      for (size_t j=0; j < pmat.size2(); j++) {
        pmat(i,j) = -1;
      }
    }

    for (size_t i=0; i < points.size(); i++) {
      size_t x = static_cast<size_t>(points[i].x);
      size_t y = static_cast<size_t>(points[i].y);

      // make sure not to overwrite medoids.
      if (!(pmat(x,y) != -1 && parts.is_medoid(pmat(x,y)))) {
        pmat((size_t)points[i].x, (size_t)points[i].y) = i;
      }
    }
  
    const double min_hbar = 80;
    out << label << " ";
    for (size_t x = 0; x <= max(max_x, min_hbar)-label.size()-1; x++) out << "-";
    out << endl;

    for (int y = (int)max_y; y >= 0; y--) {
      for (int x = 0; x <= max_x; x++) {
        int oid = pmat(x,y);
        if (oid >= 0) {
          int cid = parts.cluster_ids[oid];
          out << colors[cid % num_colors];
          out << (parts.is_medoid(oid) ? "o" : "+");
          out << None;
        } else {
          out << " ";
        }
      }
      out << endl;
    }

    for (int x = 0; x <= max(max_x, min_hbar); x++) out << "-";
    out << endl;
  }


} // namespace cluster
