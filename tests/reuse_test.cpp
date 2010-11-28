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
/// @file reuse_test.cpp
/// @author Todd Gamblin tgamblin@llnl.gov
///
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

#include "kmedoids.h"
#include "point.h"
#include "bic.h"

using namespace cluster;
using namespace std;


/// Test using points instead of QGrams to make sure clustering 
/// algorithms work.
int main(int argc, char **argv) { 
  //put 5-pt crosses inthe vector, offset by increasing distances
  vector<point> points;
  point ref(1,1);
  point stencil[] = {
    point( 0, 0),
    point( 0, 1), 
    point( 0,-1), 
    point(-1, 0), 
    point( 1, 0)
  };
  size_t stencil_size = sizeof(stencil) / sizeof(point);

  size_t num_objects = 64;
  if (argc > 1) {
    num_objects = strtol(argv[1], NULL, 0);
  }

  size_t max_clusters = num_objects / stencil_size + 5; // go 5 over to test BIC

  for (size_t i=0; points.size() < num_objects; i++) {
    for (size_t s=0; s < stencil_size && points.size() < num_objects; s++) {
      point p = ref + stencil[s];
      points.push_back(p);
    }
    ref += point(i+4, 0);
  }

  cerr << "num_objects  = " << num_objects << endl;
  cerr << "max_clusters = " << max_clusters << endl;

  cerr << points.size() << " points: ";
  copy(points.begin(), points.end(), ostream_iterator<point>(cerr, " "));
  cerr << endl;

  dissimilarity_matrix distance;
  build_dissimilarity_matrix(points, point_distance(), distance);

  cerr << "matrix size is " << distance.size1() << "x" << distance.size2() << endl;

  kmedoids km;
  for (size_t k=1; k <= max_clusters; k++) {
    km.pam(distance, k);

    ostringstream pam_msg;
    pam_msg << "PAM"
            << ", " << km.medoid_ids.size() << " clusters"
            << ", Avg. dissimilarity: " << km.average_dissimilarity()
            << ", BIC: " << bic(km, matrix_distance(distance), 2);
    
    draw(pam_msg.str(), points, km);
    cout << endl;
  }
}
