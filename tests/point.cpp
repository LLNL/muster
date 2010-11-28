//////////////////////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2010, Lawrence Livermore National Security, LLC.  
// Produced at the Lawrence Livermore National Laboratory  
// LLNL-CODE-433662
// All rights reserved.  
//
// This file is part of Muster. For details, see http://github.com/tgamblin/muster. 
// Please also read the LICENSE file for further information.
//
// Redistribution and use in source and binary forms, with or without modification, are
// permitted provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright notice, this list of
//    conditions and the disclaimer below.
//  * Redistributions in binary form must reproduce the above copyright notice, this list of
//    conditions and the disclaimer (as noted below) in the documentation and/or other materials
//    provided with the distribution.
//  * Neither the name of the LLNS/LLNL nor the names of its contributors may be used to endorse
//    or promote products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
// OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
// LAWRENCE LIVERMORE NATIONAL SECURITY, LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
// ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//////////////////////////////////////////////////////////////////////////////////////////////////

///
/// @file point.cpp
/// @author Todd Gamblin tgamblin@llnl.gov
///
#include "point.h"
#include "muster-config.h"

#ifdef MUSTER_HAVE_MPI
#include "mpi_utils.h"
#endif // MUSTER_HAVE_MPI
#include "color.h"

#include <limits>
#include <iomanip>
#include <cstdlib>
#include <boost/numeric/ublas/matrix.hpp>

using boost::numeric::ublas::matrix;
using namespace std;

namespace cluster {
  
  point::point(double _x, double _y) : x(_x), y(_y) { }

  point::point() : x(0), y(0) { }

  point::point(const point& other) : x(other.x), y(other.y) { }


#ifdef MUSTER_HAVE_MPI
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
  
  MPI_Datatype point::mpi_datatype() {
    static MPI_Datatype type = MPI_DATATYPE_NULL;
    if (type == MPI_DATATYPE_NULL) {
      MPI_Type_contiguous(2, MPI_DOUBLE, &type);
      MPI_Type_commit(&type); 
    }
    return type;
  }

#endif // MUSTER_HAVE_MPI


  ostream& operator<<(ostream& out, const point& p) {
    return out << "(" << setw(2) << p.x << "," << setw(2) << p.y << ")";
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
