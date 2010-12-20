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
/// @file random.h
/// @author Todd Gamblin tgamblin@llnl.gov
/// @brief Helper functions for taking random samples and seeding 
///        RNGs from the system clock.
///
#ifndef MUSTER_RANDOM_H
#define MUSTER_RANDOM_H

#include <sys/time.h>

namespace cluster {

  ///
  /// This is Knuth's algorithm R for taking a sample of indices from
  /// 0 to numElements.  We sample size indices from this (the superset)
  /// and put them in the subset's mapping.
  ///
  /// @param numElements    total number of elements to select from
  /// @param sample_size    number of elements to select
  /// @param out            destination for selected elements 
  /// @param random         model of STL Random Number Generator.
  ///                       must be callable as random(N), returning a random number in [0,N).
  ///
  template <class OutputIterator, class Random>
  void random_subset(size_t numElements, size_t sample_size, OutputIterator out, Random& random) {
    size_t first = 0;
    size_t remaining = numElements;
    size_t m = sample_size;
  
    while (m > 0) {
      if ((size_t)random(remaining) < m) {
        *out = first;
        ++out;
        --m;
      }
      --remaining;
      ++first;
    }
  }


  ///
  /// Returns a reasonably distributed seed for random number generators.
  /// Based on the product of the seconds and usec in gettimeofday().
  ///
  inline long get_time_seed() {
    struct timeval time;
    gettimeofday(&time, 0);
    return time.tv_sec * time.tv_usec;
  }

} // namespace cluster

#endif // MUSTER_RANDOM_H
