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
/// @file point_set.h
/// @author Todd Gamblin tgamblin@llnl.gov
///
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
