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
/// @file multi_gather_test.cpp
/// @author Todd Gamblin tgamblin@llnl.gov
/// 
#include <mpi.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <iterator>
#include <iomanip>
#include <sstream>

#include <boost/random.hpp>

#include "multi_gather.h"
#include "point.h"
#include "random.h"
#include "Timer.h"

using namespace std;
using namespace cluster;


template <class OutputIterator>
void generate_points_for_rank(int rank, OutputIterator out) {
  for (int i=0; i < (rank % 10 + 1); i++) {
    int p = rank * 10 + i;
    *out++ = point(p,p);
  }
}


int main(int argc, char **argv) {
  Timer timer;
  MPI_Init(&argc, &argv);
  
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  bool verbose = false;
  bool timing  = false;

  for (int i=1; i < argc; i++) {
    if (string(argv[i]) == "-v") {
      verbose = true;
    }

    if (string(argv[i]) == "-t") {
      timing = true;
    }
  }

  typedef boost::mt19937 random_t;
  random_t random(1);  // make sure all ranks generate the same numbers.
  boost::random_number_generator<random_t> rng(random);

  vector<point> points;   // local points to send
  vector<point> dest;     // destination vector for gathered points
  vector<int> sources;    // ranks we received from, so we can check the local points vector

  generate_points_for_rank(rank, back_inserter(points));

  timer.record("init");

  // now fire off <size> gathers, each with ~size/2 elements.
  multi_gather<point> gather(MPI_COMM_WORLD);
  for (int root=0; root < size; root++) {
    vector<int> cur_sources;
    algorithm_r(size, (int)ceil(sqrt((double)size)), back_inserter(cur_sources), rng);
    gather.start(points.begin(), points.end(), cur_sources.begin(), cur_sources.end(), dest, root);

    if (rank == root) {
      // record sources so we can check later.
      cur_sources.swap(sources);
    }
  }

  timer.record("start_gathers");
  
  if (verbose) {
    cerr << rank << " gathering from [";
    for (size_t i=0; i < sources.size(); i++) cerr << setw(3) << sources[i] << " ";
    cerr << "]" << endl;
    timer.record("verbose");
  }

  if (verbose) cerr << rank << " calling multi_gather::finish()." << endl;

  gather.finish();
  timer.record("finish_gathers");

  if (verbose) cerr << rank << " finished gather." << endl;
  
  int passed = 1;
  size_t index = 0;
  for (size_t i=0; passed && i < sources.size(); i++) {
    vector<point> expected;
    generate_points_for_rank(sources[i], back_inserter(expected));
    
    for (size_t j=0; j < expected.size(); j++) {
      if (index >= dest.size() || dest[index] != expected[j]) {
        passed = 0;
      }
      index++;
    }
  }
  timer.record("check");

  if (verbose) {
    ostringstream msg;

    msg << rank << " Expected: ";
    for (size_t i=0; i < sources.size(); i++) {
      generate_points_for_rank(sources[i], ostream_iterator<point>(msg, " "));
    }
    msg << endl;

    msg << rank << " Found:   ";
    for (size_t i=0; i < dest.size(); i++) {
      msg << " " << dest[i];
    }
    msg << endl;
    cerr << msg.str();
    timer.record("verbose");
  }

  int num_passed = 0;
  MPI_Allreduce(&passed, &num_passed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Finalize();

  bool all_passed = (num_passed == size);
  if (rank == 0) {
    cerr << (all_passed ? "PASSED" : "FAILED") << endl;
  }

  if (timing && rank == 0) timer.dump(cerr);

  return all_passed ? 0 : 1;
}
