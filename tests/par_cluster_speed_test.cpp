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
/// @file par_cluster_speed_test.cpp
/// @author Todd Gamblin tgamblin@llnl.gov
/// 
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <sys/time.h>

#include <boost/random.hpp>

#include "timing.h"
#include "point.h"
#include "bic.h"
#include "par_kmedoids.h"

using namespace cluster;
using namespace std;


void usage() {
  cerr << "Usage: par-cluster-speed-test [-htvi] [-n num-points] [-c clusters] [-s scale]" << endl;
  cerr << "  Compare parallel clustering with sequential clustering." << endl;
  cerr << "Options:" << endl;
  cerr << "  -h         Show this message." << endl;
  cerr << "  -x         Use BIC-scored versions of PAM and CAPEK." << endl;
  cerr << "  -t         Save details timing info in a file." << endl;
  cerr << "  -n         Number of points per process." << endl;
  cerr << "               Default is 1." << endl;
  cerr << "  -i         Initial sample size in CAPEK (before 2*k is added)." << endl;
  cerr << "               Default is 40." << endl;
  cerr << "  -r         Number of trials per k in CAPEK." << endl;
  cerr << "               Default is 5." << endl;
  cerr << "  -k         Max number of clusters to search for." << endl;
  cerr << "               Default is 10." << endl;
  exit(1);
}

size_t objects_per_process = 1;
size_t num_clusters = 10;
size_t init_size = 40;
size_t max_reps = 5;
bool use_bic = false;
bool timing = false;

/// Uses getopt to read in arguments.
void get_args(int *argc, char ***argv, int rank) {
  int c;
  char *err;

  while ((c = getopt(*argc, *argv, "htxn:i:r:k:")) != -1) {
    switch (c) {
    case 'h':
      if (rank == 0) usage();
      exit(1);
      break;
    case 'x':
      use_bic = true;
      break;
    case 't':
      timing = true;
      break;
    case 'n':
      objects_per_process = strtol(optarg, &err, 0);
      if (*err) usage();
      break;
    case 'i':
      init_size = strtol(optarg, &err, 0);
      if (*err) usage();
      break;
    case 'r':
      max_reps = strtol(optarg, &err, 0);
      if (*err) usage();
      break;
    case 'k':
      num_clusters = strtol(optarg, &err, 0);
      if (*err) usage();
      break;
    default:
      if (rank == 0) usage();
      exit(1);
      break;
    }
  }

  // adjust params
  *argc -= optind;
  *argv += optind;
}




int main(int argc, char **argv) { 
  MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  get_args(&argc, &argv, rank);

  // generator to make points.  Seed differently on each rank.
  typedef boost::mt19937 random_t;
  random_t random(get_time_seed() + rank);  
  boost::random_number_generator<random_t> rng(random);

  // vector of local points
  vector<point> points;
  
  // generate randomly distributed, zero-centered points.
  // @todo use random gaussian points?
  for (size_t i=0; i < objects_per_process; i++) {
    int x = rng(5000+1) - 2500;
    int y = rng(5000+1) - 2500;
    points.push_back(point(x,y));
  }

  par_kmedoids parkm;
  parkm.set_init_size(init_size);
  parkm.set_max_reps(max_reps);

  // trials of whole algorithm, to account for any variability
  const size_t trials = 10;

  long long start = get_time_ns();
  for (size_t i=0; i < trials; i++) {
    parkm.capek(points, point_distance(), num_clusters);
  }

  double total = get_time_ns() - start;
  double avg = total / trials;
  
  if (rank == 0) {
    if (timing) {
      ostringstream timing_filename;
      timing_filename << "times-" << size;
      ofstream timing_file(timing_filename.str().c_str());
      parkm.get_timer().write(timing_file);
    }

    cout << "PROCS:   " << size << " " << endl;
    cout << "TOTAL:   " << total / 1e9 << endl;
    cout << "AVERAGE: " << avg   / 1e9 << endl;
  }
}
