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
/// @file gaussian.h
/// @author Todd Gamblin tgamblin@llnl.gov
/// 
#ifndef GAUSSIAN_GENERATOR_H
#define GAUSSIAN_GENERATOR_H

#include <cfloat>
#include <sys/time.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

#include "point.h"
#include "random.h"

namespace cluster {
  ///
  /// Gaussian cluster with cutoff radius.
  ///
  class gaussian_generator_2d {
    typedef boost::mt19937 rand_type;

  public:    
    gaussian_generator_2d(double _x, double _y, double stddev, double xscale=1.0, double yscale=1.0) 
      : x(_x * xscale), y(_y * yscale), xs(xscale), ys(yscale), normal(0, stddev), radius(rand, normal) { 
      rand.seed(get_time_seed());
    }

    void seed(rand_type::result_type s) {
      rand.seed(s);
    }

    point next_point() {
      double r = radius();
      double theta = rand() / (double)rand.max() * 2 * M_PI;
      double px = x + r * cos(theta);
      double py = y + r * sin(theta);
      
      return point(px * xs, py * ys);
    }

    point center() {
      return point(x,y);
    }

  private:
    double x, y;
    double xs, ys;
    rand_type rand;
    
    typedef boost::normal_distribution<double> dist_type;
    dist_type normal;
    boost::variate_generator<rand_type, dist_type> radius;
  };


} // namespace

#endif // GAUSSIAN_GENERATOR_H
