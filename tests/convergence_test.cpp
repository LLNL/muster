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
/// @file convergence_test.cpp
/// @author Todd Gamblin tgamblin@llnl.gov
/// 
#include <iostream>
#include <iomanip>

#include "kmedoids.h"
#include "matrix_utils.h"
#include "point.h"
#include "point_set.h"
#include "Timer.h"

using namespace cluster;
using namespace std;

struct test {
  const char *points;
  size_t k;
};

int main(int argc, char **argv) {
  const test tests[] = {
    {"(1,  1) (1,  2) (0,  1) (2,  1) (19,  1) (19,  2) (19,  0) (18,  1) (38,  0)"
     "(39, 1) (58, 1) (58, 0) (79, 2) (101, 1) (100, 1) (124, 2) (124, 0)",
     3
     },
    {"( 1, 1) ( 1, 2) ( 0, 1) ( 2, 1) ( 5, 1) ( 5, 2) ( 5, 0) ( 4, 1) ( 6, 1) (10, 1)"
     "(10, 2) (10, 0) ( 9, 1) (16, 1) (16, 2) (16, 0) (15, 1) (17, 1) (23, 2) (23, 0)"
     "(31, 1) (31, 0) (30, 1) (32, 1) (40, 1) (40, 2) (39, 1) (41, 1) (50, 1) (50, 0)"
     "(49, 1) (51, 1) (61, 2) (61, 0) (62, 1) (72, 1) (74, 1) (86, 1) (86, 2) (86, 0)"
     "(87, 1) (100, 1) (100, 2) (100, 0) (99, 1) (115, 1) (115, 2) (115, 0)",
     4},
    {"( 1, 1) ( 1, 2) ( 1, 0) ( 0, 1) ( 2, 1) ( 5, 1) ( 5, 2) ( 5, 0) ( 6, 1) (10, 1) (10, 2) (10, 0)"
     "( 9, 1) (16, 1) (16, 2) (16, 0) (15, 1) (17, 1) (23, 2) (22, 1) (24, 1) (31, 1) (31, 2) (40, 0)"
     "(39, 1) (50, 1) (50, 2) (51, 1) (61, 1) (61, 2) (61, 0) (60, 1) (73, 1) (73, 2) (73, 0) (72, 1)"
     " (86, 1) (100, 1) (100, 2) (100, 0) (101, 1) (115, 2) (115, 0) (114, 1)",
     6},
    {"( 1, 1) ( 1, 2) ( 1, 0) ( 0, 1) ( 2, 1) ( 5, 1) ( 5, 2) ( 5, 0) ( 4, 1) ( 6, 1) (10, 1) (10, 2)"
     "(10, 0) ( 9, 1) (11, 1) (16, 1) (16, 2) (16, 0) (15, 1) (17, 1) (23, 1) (23, 2) (23, 0) (22, 1)"
     "(24, 1) (31, 1) (31, 2) (31, 0) (30, 1) (32, 1) (40, 1) (40, 2) (40, 0) (39, 1) (41, 1) (50, 1)"
     "(50, 2) (50, 0) (49, 1) (51, 1) (61, 1) (61, 2) (61, 0) (60, 1) (62, 1) (73, 1) (73, 2) (73, 0)"
     "(72, 1) (74, 1) (86, 1) (86, 2) (86, 0) (85, 1) (87, 1) (100, 1) (100, 2) (100, 0) (99, 1) "
     "(101, 1) (115, 1) (115, 2) (115, 0) (114, 1)", 
     6}
  };
  const size_t num_tests = sizeof(tests) / sizeof(test);

  for (size_t i=0; i < num_tests; i++) {
    point_set points;
    points.parse_points(tests[i].points);
    const size_t k = tests[i].k;

    kmedoids cluster;
    dissimilarity_matrix distance;
    build_dissimilarity_matrix(points.points(), point_distance(), distance);
    
    cerr << "Starting Test " << i << " with k=" << k <<endl;
    cluster.pam(distance, k);
    cerr << "Finished."    << endl;

    ostringstream label;
    label << "Trial " << i;
    draw(label.str(), points.points(), cluster);
  }
}
