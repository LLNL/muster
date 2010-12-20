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
/// @file point.h
/// @author Todd Gamblin tgamblin@llnl.gov
///
#ifndef MUSTER_TEST_POINT_H
#define MUSTER_TEST_POINT_H

#include "muster-config.h"

#ifdef MUSTER_HAVE_MPI
#include <mpi.h>
#endif // MUSTER_HAVE_MPI 

#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include "partition.h"

namespace cluster {
  
  /// Simple 2 dimensional point class for testing medoids algorithms.
  struct point {
  public:
    double x, y;
    
    /// New point with position (x,y)
    point(double x, double y);

    /// New point at (0,0)
    point();

    /// Copy constructor
    point(const point& other);

    // Distance between this point and another. 
    double distance(const point& other) const {
      double dx = other.x - x;
      double dy = other.y - y;
      return ::sqrt(dx*dx + dy*dy);
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

    void operator/=(double divisor) {
      x /= divisor;
      y /= divisor;
    }

    point operator/(double divisor) const {
      point p(x,y);
      p /= divisor;
      return p;
    }

#ifdef MUSTER_HAVE_MPI
    /// Returns the size of a packed point
    int packed_size(MPI_Comm comm) const;
  
    /// Packs a point into an MPI packed buffer
    void pack(void *buf, int bufsize, int *position, MPI_Comm comm) const;

    /// Unpacks a point from an MPI packed buffer
    static point unpack(void *buf, int bufsize, int *position, MPI_Comm comm);
    
    /// Get an MPI datatype for a point.
    static MPI_Datatype mpi_datatype();
#endif // MUSTER_HAVE_MPI
  };

  std::ostream& operator<<(std::ostream& out, const point& p);
    
  ///
  /// Draws a set of points in ascii with console colors.  Colors are assigned based on
  /// the partition provided.  Indices in points vector should correspond to ids in the 
  /// partition.
  /// 
  void draw(std::string label, std::vector<point>& points, const cluster::partition& parts, 
            std::ostream& out = std::cout);


  /// Distance bt/w two points
  struct point_distance {
    double operator()(const point& left, const point& right) const {
      return left.distance(right);
    }
  };

} // namespace cluster

#endif // MUSTER_TEST_POINT_H

