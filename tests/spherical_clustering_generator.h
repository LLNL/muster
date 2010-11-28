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
/// @file spherical_clustering_generator.cpp
/// @author Todd Gamblin tgamblin@llnl.gov
///
#ifndef SPHERICAL_CLUSTERING_GENERATOR_H
#define SPHERICAL_CLUSTERING_GENERATOR_H

#include <vector>
#include "point.h"
#include "gaussian.h"
#include "random.h"
#include <climits>


namespace cluster {
  
  class spherical_clustering_generator {
    std::vector<gaussian_generator_2d> generators;
    double default_stddev;
    double xscale;
    double yscale;
    size_t cur_generator;

  public:
    spherical_clustering_generator() 
      : default_stddev(1.0), xscale(1.0), yscale(1.0), cur_generator(0) { }

    void add_cluster(point center, double stddev = -1) {
      if (stddev < 0) stddev = default_stddev;
      generators.push_back(gaussian_generator_2d(center.x, center.y, stddev, xscale, yscale));
    }
    
    /// Get the next point from the next cluster in round-robin order.
    point next_point() {
      if (!generators.size()) {
        return point(0,0);
      }      
      point next = generators[cur_generator].next_point();
      cur_generator = (cur_generator + 1) % generators.size();
      return next;
    }

    /// Get the next point for a specific cluster
    point next_point(size_t i) {
      return generators[i].next_point();
    }
  
    void set_default_stddev(double dsd) {
      default_stddev = dsd;
    }
    
    void set_scale(double scale) {
      xscale = yscale = scale;
    }
    
    void set_xscale(double scale) {
      xscale = scale;
    }

    void set_yscale(double scale) {
      yscale = scale;
    }

    size_t size() const {
      return generators.size();
    }
    
  };

} // namespace cluster


#endif // SPHERICAL_CLUSTERING_GENERATOR_H
