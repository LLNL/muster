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
/// @file par_cluster_test.cpp
/// @author Todd Gamblin tgamblin@llnl.gov
/// 
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

#include "point.h"
#include "bic.h"
#include "par_kmedoids.h"

using namespace cluster;
using namespace std;

void usage() {
  cerr << "Usage: par-cluster-test [-hxv] [-n num-points] [-k clusters] [-i initial-size] [-r reps] [-s seed]" << endl;
  cerr << "  Compare parallel clustering with sequential clustering." << endl;
  cerr << "Options:" << endl;
  cerr << "  -h         Show this message." << endl;
  cerr << "  -x         Use BIC-scored versions of clustering algorithms." << endl;
  cerr << "  -v         Verbose mode.  Draws actual clusterings and outputs timings." << endl;
  cerr << "  -n         Number of points per process." << endl;
  cerr << "               Default is 1." << endl;
  cerr << "  -k         Max number of clusters to search for." << endl;
  cerr << "               Default is number of processes * points per process." << endl;
  cerr << "  -i         Initial sample size in CAPEK (before 2*k is added)." << endl;
  cerr << "               Default is 40." << endl;
  cerr << "  -r         Number of repeated trials per k in CAPEK." << endl;
  cerr << "               Default is 5." << endl;
  cerr << "  -s         Seed value for random number generator at start of test." << endl;
  cerr << "               Default will be generated based on gettimeofday, see documentation." << endl;
  exit(1);
}

size_t objects_per_process = 1;
size_t num_clusters = 0;
size_t init_size = 40;
size_t max_reps = 5;
bool use_bic = false;
bool verbose = false;
bool use_seed = false;
size_t seed = 0;

/// Uses getopt to read in arguments.
void get_args(int *argc, char ***argv, int rank) {
  int c;
  char *err;

  while ((c = getopt(*argc, *argv, "hxn:i:r:k:vs:")) != -1) {
    switch (c) {
    case 'h':
      if (rank == 0) usage();
      exit(1);
      break;
    case 'x':
      use_bic = true;
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
    case 'v':
      verbose = true;
      break;
    case 's':
      use_seed = true;
      seed = strtol(optarg, &err, 0);
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


vector<point> points;
kmedoids km;
par_kmedoids parkm(MPI_COMM_WORLD);

void print_cluster_info(const cluster::partition& gathered, const dissimilarity_matrix& distance) {
  cout << "seq k: " << km.medoid_ids.size()
       << ", par k: " << gathered.medoid_ids.size()
       << ", Mirkin distance: " << setprecision(3) << mirkin_distance(km, gathered)
       << endl;
  
  if (verbose) {
    ostringstream pam_msg;
    pam_msg << "PAM"
            << ", " << km.medoid_ids.size() << " clusters"
            << ", Avg. dissimilarity: " << km.average_dissimilarity()
            << ", BIC: " << bic(km, matrix_distance(distance), 2);
    
    ostringstream parkm_msg;
    parkm_msg << "CAPEK" 
              << ", " << gathered.medoid_ids.size() << " clusters"
              << ", Avg. dissimilarity: " << parkm.average_dissimilarity()
              << ", BIC: " << parkm.bic_score();
    
    draw(pam_msg.str(), points, km);
    draw(parkm_msg.str(), points, gathered);
    cout << endl;
    parkm.get_timer().write();
    cout << endl;
  }
}


/// Test using points instead of QGrams to make sure clustering 
/// algorithms work.
int main(int argc, char **argv) { 
  MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  get_args(&argc, &argv, rank);

  //put 5-pt crosses inthe vector, offset by increasing distances
  point ref(1,1);
  point stencil[] = {
    point( 0, 0),
    point( 0, 1), 
    point( 0,-1), 
    point(-1, 0), 
    point( 1, 0)
  };
  size_t stencil_size = sizeof(stencil) / sizeof(point);


  size_t num_objects = size * objects_per_process;
  if (num_clusters == 0) {
    num_clusters = num_objects / stencil_size + 5;
  }

  num_clusters = min(num_clusters, num_objects);

  vector<point> my_points;
  for (size_t i=0; points.size() < num_objects; i++) {
    for (size_t s=0; s < stencil_size && points.size() < num_objects; s++) {
      point p = ref + stencil[s];
      if (rank == (int)(points.size() / objects_per_process)) {
        my_points.push_back(p);
      }
      points.push_back(p);
    }
    ref += point(i+4, 0);
  }

  dissimilarity_matrix distance;
  build_dissimilarity_matrix(points, point_distance(), distance);

  parkm.set_init_size(init_size);
  parkm.set_max_reps(max_reps);
  if (use_seed) {
      parkm.set_seed(seed);
  }

  for (size_t k=1; k <= num_clusters; k++) {
    vector<point> medoids;
    cluster::partition gathered;


    km.pam(distance, k);
    parkm.capek(my_points, point_distance(), k, &medoids);
    parkm.gather(gathered);

    if (rank == 0) {
      cout << "Max k = " << k << endl;
      cout << "  No BIC:    ";
      print_cluster_info(gathered, distance);
    }

    km.xpam(distance, k, 2);
    parkm.xcapek(my_points, point_distance(), k, 2, &medoids);
    parkm.gather(gathered);

    if (rank == 0) {
      cout << "  Using BIC: ";
      print_cluster_info(gathered, distance);
      cout << endl;
    }
  }
}

