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
/// @file pont_cluster_test.cpp
/// @author Todd Gamblin tgamblin@llnl.gov
///
#include <vector>
#include <iostream>
#include <cstdlib>
using namespace std;

#include "point.h"
#include "dissimilarity.h"
#include "kmedoids.h"
#include "color.h"
#include "point.h"
using namespace cluster;


void usage() {
  cerr << "Usage: point-cluster-test [-hv] [-c clusters] [-k max_k] [-i initial-size] [-r reps]" << endl;
  cerr << "  Compare sequential clustering techniques with regular pattern of 5-point clusters." << endl;
  cerr << "Options:" << endl;
  cerr << "  -h         Show this message." << endl;
  cerr << "  -v         Verbose mode.  Draws actual clusterings and outputs timings." << endl;
  cerr << "  -c         Number of clusters to generate." << endl;
  cerr << "               Default is 5." << endl;
  cerr << "  -k         Max number of clusters to search for." << endl;
  cerr << "               Default is same as number of clusters generated." << endl;
  cerr << "  -i         Initial sample size in clara (before 2*k is added)." << endl;
  cerr << "               Default is 40." << endl;
  cerr << "  -r         Number of repeated trials per k in clara." << endl;
  cerr << "               Default is 5." << endl;
  exit(1);
}

size_t num_clusters = 5;
int max_k = -1;
size_t init_size = 40;
size_t max_reps = 5;
bool verbose = false;

/// Uses getopt to read in arguments.
void get_args(int *argc, char ***argv) {
  int c;
  char *err;

  while ((c = getopt(*argc, *argv, "hvc:k:i:r:")) != -1) {
    switch (c) {
    case 'h':
      usage();
      break;
    case 'v':
      verbose = true;
      break;
    case 'c':
      num_clusters = strtol(optarg, &err, 0);
      if (*err) usage();
      break;
    case 'k':
      max_k = strtol(optarg, &err, 0);
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
    default:
      usage();
      break;
    }
  }

  // adjust params
  *argc -= optind;
  *argv += optind;
}


int main(int argc, char **argv) {   
  get_args(&argc, &argv);

  cerr << "max k = " << max_k << endl;

  if (max_k < 0) {
    max_k = num_clusters;
  }

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

  vector<point> points;  //vector of test points
  for (size_t i=0; i < num_clusters; i++) {
    for (size_t s=0; s < stencil_size; s++) {
      point p = ref + stencil[s];
      points.push_back(p);
    }
    ref += point(i+4, 0);
  }

  cout << "Testing with " << points.size() << " points "
       << "with " << num_clusters << " clusters "
       << "for k = 1 to " << max_k 
       << endl;

  dissimilarity_matrix dmatrix;
  build_dissimilarity_matrix(points, point_distance(), dmatrix);

  kmedoids km, clara, xkm, xclara;
  for (size_t k = 1; k <= (size_t)max_k; k++) {
    km.pam(dmatrix, k);
    clara.clara(points, point_distance(), k);
    xkm.xpam(dmatrix, k, 2);
    xclara.xclara(points, point_distance(), k, 2);

    cout << "k: " << k << endl;
    cout << "  mirkin(pam, clara):   " << mirkin_distance(km, clara) << endl;
    cout << "  mirkin(xpam, pam):    " << mirkin_distance(xkm, km) << endl;
    cout << "  mirkin(xpam, clara):  " << mirkin_distance(xkm, clara) << endl;
    cout << "  mirkin(xpam, xclara): " << mirkin_distance(xkm, xclara) << endl;

    if (verbose) {
      draw("PAM", points, km);
      draw("CLARA", points, clara);
      draw("XPAM", points, xkm);
      draw("XCLARA", points, xclara);
    }

    cout << endl;
  }
}
