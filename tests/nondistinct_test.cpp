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
/// @file nondistinct_test.cpp
/// @author Todd Gamblin tgamblin@llnl.gov
/// 
#include <iostream>
#include <iomanip>

#include "kmedoids.h"
#include "matrix_utils.h"
#include "point.h"
#include "Timer.h"

using namespace cluster;
using namespace std;


int main(int argc, char **argv) {
  size_t max_k = 10;
  if (argc > 1) {
    max_k = atoi(argv[1]);
  }

  // make a set of points with some duplicate elements
  std::vector<point> points;
  for (size_t i=0; i < 16; i++) {
    size_t type = i % 3;
    switch (type) {
    case 0:
      points.push_back(point(1,1));
      break;
    case 1:
      points.push_back(point(3,3));
      break;
    case 2:
      points.push_back(point(4,4));
      break;
    }
  }

  dissimilarity_matrix distance;
  build_dissimilarity_matrix(points, point_distance(), distance);

  kmedoids cluster;
  for (size_t k=1; k <= 10; k++) {
    cluster.pam(distance, k);

    vector<object_id> medoids(cluster.medoid_ids);
    sort(medoids.begin(), medoids.end());
    size_t num_unique = unique(medoids.begin(), medoids.end()) - medoids.begin();
    if (num_unique != cluster.medoid_ids.size()) {
      cerr << "Error: medoids are not distinct:" << endl;
      cerr << "[";
      copy(cluster.medoid_ids.begin(), cluster.medoid_ids.end(), 
           ostream_iterator<object_id>(cerr, " "));
      cerr << "]" << endl;

      cout << "FAILED" << endl;
      exit(1);
    }

    // ensure that each cluster has at least one element
    cluster_list clist;
    cluster.to_cluster_list(clist);
    for (size_t i=0; i < clist.size(); i++) {
      if (!clist[i].size()) {
        cerr << "Error: clustering contains empty clusters." << endl;
        cerr << cluster << endl;
        cout << "FAILED" << endl;
        exit(1);
      }
    }
  }

  cout << "PASSED" << endl;
  exit(0);
}
