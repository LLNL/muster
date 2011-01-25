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
/// @file random_test.cpp
/// @author Todd Gamblin tgamblin@llnl.gov
/// 
#include <stdint.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <boost/random.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "random.h"
#include "matrix_utils.h"
#include "timing.h"

using namespace cluster;
using namespace std;
using boost::numeric::ublas::matrix;

const size_t reps = 5;


int main(int argc, char **argv) {
  uint32_t seed = get_time_seed();
  typedef boost::mt19937 random_t;
  random_t random;
  boost::random_number_generator<random_t> rng(random);

  // ranges of powers of two to test.
  size_t startN = 8, endN = 22;
  size_t starts = 5;

  // Put labels in the 0th column and 0th row.
  matrix<timing_t> algorithm_r_timings(endN - startN+1, endN - starts);
  fill(algorithm_r_timings.begin1(), algorithm_r_timings.end1(), 0);
  vector<size_t> results;
  for (size_t N=startN; N < endN; N++) {
    algorithm_r_timings(N-startN+1, 0) = 1<<N;
    for (size_t s=starts; s < N; s++) {
      algorithm_r_timings(0, s-starts+1) = 1<<s;
    }
  }
  matrix<timing_t> fast_timings(algorithm_r_timings);
  matrix<double>   speedup(algorithm_r_timings);

  // Run sampling algorithms for different values of s and N
  for (size_t N=startN; N < endN; N++) {
    for (size_t s=starts; s < N; s++) {
      timing_t start;
      size_t Nx = N-startN+1, sx = s-starts+1;

      // reserve space for the numbers so we don't include memory allocation time
      results.reserve(1 << s);

      // run algorithm r
      random.seed(seed);
      start = get_time_ns();
      for (size_t i=0; i < reps; i++) {
        results.clear();
        algorithm_r(1 << N, 1 << s, back_inserter(results), rng);
      }
      algorithm_r_timings(Nx, sx) = get_time_ns() - start;

      // run fast sampling algorithm
      random.seed(seed);
      start = get_time_ns();
      for (size_t i=0; i < reps; i++) {
        results.clear();
        fast_sample(1 << N, 1 << s, back_inserter(results), rng);
      }
      fast_timings(Nx, sx) = get_time_ns() - start;
      
      // record how much faster the fast sampling algorithm was.
      speedup(Nx, sx) =  algorithm_r_timings(Nx, sx) / (double)fast_timings(Nx, sx);
    }
  }

  // matrix of algorithm r timings
  cout << endl;
  output(algorithm_r_timings);

  // matrix of fast sample timings
  cout << endl;
  output(fast_timings);

  // matrix of speedup values
  cout << endl;
  output(speedup);
}
