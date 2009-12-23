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
using boost::numeric::ublas::matrix;
using namespace std;

#include "string_utils.h"
using namespace stringutils;

namespace cluster {

  point::point(int _x, int _y) : x(_x), y(_y) { }

  point::point(const point& other) : x(other.x), y(other.y) { }


#ifdef HAVE_MPI
  /// Returns the size of a packed point
  int point::packed_size(MPI_Comm comm) const {
    return 2 * mpi_packed_size(1, MPI_INT, comm);
  }
  
  /// Packs a point into an MPI packed buffer
  void point::pack(void *buf, int bufsize, int *position, MPI_Comm comm) const {
    PMPI_Pack(const_cast<int*>(&x), 1, MPI_INT, buf, bufsize, position, comm);
    PMPI_Pack(const_cast<int*>(&y), 1, MPI_INT, buf, bufsize, position, comm);
  }

  /// Unpacks a point from an MPI packed buffer
  point point::unpack(void *buf, int bufsize, int *position, MPI_Comm comm) {
    point p;
    PMPI_Unpack(buf, bufsize, position, &p.x, 1, MPI_INT,  comm);
    PMPI_Unpack(buf, bufsize, position, &p.y, 1, MPI_INT,  comm);
    return p;
  }
#endif // HAVE_MPI


  ostream& operator<<(ostream& out, const point& p) {
    return out << "(" << setw(2) << p.x << "," << setw(2) << p.y << ")";
  }


  void parse_points(const string& str, vector<point>& points) {
    string trimmed = trim(str, "() ");

    vector<string> parts;
    split(trimmed, "(,", parts);

    vector<int> values(parts.size());
    for (size_t i=0; i < parts.size(); i++) {
      parts[i]  = trim(parts[i], "() ,");
      values[i] = strtol(parts[i].c_str(), NULL, 0);
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

    int max_x=0;
    int max_y=0;
    for (size_t i=0; i < points.size(); i++) {
      max_x = max(points[i].x, max_x);
      max_y = max(points[i].y, max_y);
    }
  
    matrix<int> pmat(max_x+1, max_y+1);
    for (size_t i=0; i < pmat.size1(); i++) {
      for (size_t j=0; j < pmat.size2(); j++) {
        pmat(i,j) = -1;
      }
    }

    for (size_t i=0; i < points.size(); i++) {
      pmat(points[i].x, points[i].y) = i;
    }
  
    const int min_hbar = 80;
    out << label << " ";
    for (size_t x = 0; x <= max(max_x, min_hbar)-label.size()-1; x++) out << "-";
    out << endl;

    for (int y = max_y; y >= 0; y--) {
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
