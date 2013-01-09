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
/// @file point_set.cpp
/// @author Todd Gamblin tgamblin@llnl.gov
///
#include "point_set.h"

#include <limits>
#include <algorithm>

#include <boost/algorithm/string.hpp>

#ifdef MUSTER_HAVE_MPI
#include "mpi_utils.h"
#include "mpi_bindings.h"
#endif // MUSTER_HAVE_MPI

using namespace boost;
using namespace std;


namespace cluster {

  point_set::point_set()
    : min_x_(numeric_limits<double>::min()),
      max_x_(numeric_limits<double>::max()),
      min_y_(numeric_limits<double>::min()),
      max_y_(numeric_limits<double>::max())
  { }

  point_set::~point_set() { }

  void point_set::add_point(point p) {
    min_x_ = min(p.x, min_x_);
    max_x_ = max(p.x, max_x_);
    min_y_ = min(p.y, min_y_);
    max_y_ = max(p.y, max_y_);

  }

  void point_set::normalize() {
    for (size_t i=0; i < points_.size(); i++) {
      points_[i].x = (points_[i].x - min_x_)/(max_x_ - min_x_);
      points_[i].y = (points_[i].y - min_y_)/(max_y_ - min_y_);
    }
  }

  void point_set::parse_points(const string& str) {
    string trimmed = trim_copy_if(str, is_any_of("() "));

    vector<string> parts;
    split(parts, trimmed, is_any_of("(,"));

    vector<double> values(parts.size());
    for (size_t i=0; i < parts.size(); i++) {
      parts[i]  = trim_copy_if(parts[i], is_any_of("() ,"));
      values[i] = strtod(parts[i].c_str(), NULL);
    }

    for (size_t i=1; i < values.size(); i += 2) {
      points_.push_back(point(values[i-1], values[i]));
    }
  }

  void point_set::parse_point_csv(const string& str) {
    double X, Y;
    char  *err;
    string trimmed = trim_copy_if(str, is_any_of(","));

    vector<string> parts;
    split(parts, str, is_any_of(","));

    if (parts.size() < 2)
      return;

    X = strtod(parts[0].c_str(), &err);
    if (*err) return;

    Y = strtod(parts[1].c_str(), &err);
    if (*err) return;

    points_.push_back(point(X,Y));
  }


  void point_set::load_csv_file(istream& in) {
    string line;
    while (in.good()) {
      getline(in, line);
      parse_point_csv(line);
    }
  }


  void point_set::write_csv_file(ostream& out, partition *clustering) {
    out.setf(std::ios::fixed);
    for (size_t i = 0; i < points_.size(); i++) {
      out << points_[i].x << ", " << points_[i].y;
      if (clustering) {
        out << ", " << clustering->cluster_ids[i];
      }
      out << endl;
    }
  }

#ifdef MUSTER_HAVE_MPI
  void scatter(point_set& points, int root, MPI_Comm comm) {
    int rank, size;
    CMPI_Comm_rank(comm, &rank);
    CMPI_Comm_size(comm, &size);

    size_t points_per_proc = points.size() / size;
    if (rank == root && (points_per_proc * size != points.size())) {
      cerr << "Size must evenly divide points size." << endl;
      exit(1);
    }

    CMPI_Scatter(&points[0], points_per_proc, point::mpi_datatype(),
                 &points[0], points_per_proc, point::mpi_datatype(),
                 root, comm);
  }
#endif // MUSTER_HAVE_MPI
} // namespace cluster
