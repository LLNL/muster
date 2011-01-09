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
#include <tr1/unordered_map>

namespace cluster {

  ///
  /// This is Knuth's algorithm R for taking a sample of numElements numbers.
  /// We sample <code>sample_size</code> indices from [0..numElements) and write them to
  /// <code>out</code>.
  ///
  /// Note that this algorithm scales linaerly with the number of objects sampled from
  /// (that is, numElements).
  ///
  /// @sa fast_sample()
  ///
  /// @param numElements    total number of elements to select from
  /// @param sample_size    number of elements to select
  /// @param out            destination for selected elements, must model output iterator
  /// @param random         model of STL Random Number Generator.
  ///                       must be callable as random(N), returning a random number in [0,N).
  ///
  template <class OutputIterator, class Random>
  void algorithm_r(size_t numElements, size_t sample_size, OutputIterator out, Random& random) {
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
  /// This is a fast algorithm for random sampling that scales with the number of elements 
  /// sampled (sample_size).  Its performance does not depend on the number of elements sampled
  /// from (numElements), so it will perform better than algorithm R in most cases.
  ///
  /// Specifically, this algorithm will perform better when sample_size is small compared to 
  /// numElements. Algorithm R will perform better when sample_size is a large percentage of
  /// the elements to be sampled.  Use this algoritm when you need O(sample_size) performance
  /// and sample_size is constant in the limit with respect to num_elements.
  ///
  /// @param numElements    total number of elements to select from
  /// @param sample_size    number of elements to select
  /// @param out            destination for selected elements, must model output iterator.
  /// @param random         model of STL Random Number Generator.
  ///                       must be callable as random(N), returning a random number in [0,N).
  ///
  template <class OutputIterator, class Random>
  void fast_sample(size_t numElements, size_t sample_size, OutputIterator out, Random& random) {
    // Map keeps track of numbers whose identities have been swapped.
    std::tr1::unordered_map<size_t, size_t> swaps;
    
    for (size_t s = 0; s < sample_size; s++) {
      size_t remaining = numElements - s;  // number of elements remaining to be picked
      size_t pick = random(remaining);     // random pick from remaining elts

      *out = swaps.count(pick) ? swaps[pick] : pick;
      ++out;

      // Make as-yet unpicked number eligible to be picked next time
      size_t last = remaining - 1;
      swaps[pick] = swaps.count(last) ? swaps[last] : last;
    }
  }
  

  ///
  /// Returns a seed for random number generators based on the product
  /// of sec and usec from gettimeofday().
  ///
  inline long get_time_seed() {
    struct timeval time;
    gettimeofday(&time, 0);
    return time.tv_sec * time.tv_usec;
  }

} // namespace cluster

#endif // MUSTER_RANDOM_H
