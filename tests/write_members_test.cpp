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
/// @file write_members_test.cpp
/// @author Todd Gamblin tgamblin@llnl.gov
///
#include <iostream>
#include <iterator>
#include <sstream>
#include <algorithm>
#include <string.h>
#include "partition.h"

using namespace cluster;
using namespace std;

static bool passed = true;

template <class T, class Iterator2>
void validate(T* expected_begin, T* expected_end, Iterator2 actual_begin) {
  if (!equal(expected_begin, expected_end, actual_begin)) {
    passed = false;
    cerr << "expected: ";
    copy(expected_begin, expected_end, ostream_iterator<T>(cerr, " "));
    cerr << endl;
    cerr << "found:    ";
    int len = expected_end - expected_begin;
    copy(actual_begin, actual_begin + len, ostream_iterator<T>(cerr, " "));
    cerr << endl;
  }
}


int main(int argc, char **argv) {
  cluster::partition p;

  // indices         0  1  2  3  4  5  6  7  8  9 10 11 12 13
  medoid_id ids[] = {0, 0, 0, 1, 0, 1, 1, 1, 2, 1, 2, 2, 2, 0};
  const size_t num_ids = sizeof(ids) / sizeof(medoid_id);
  
  object_id medoids[] = {0, 3, 8};
  const size_t num_meds = sizeof(medoids) / sizeof(object_id);
  for (size_t i=0; i < num_meds; i++) {
    p.medoid_ids.push_back(medoids[i]); 
  }
  for (size_t i=0; i < num_ids; i++) {
    p.cluster_ids.push_back(ids[i]);
  }
  

  vector<medoid_id> out_ids;
  p.write_members(0, back_inserter(out_ids));
  medoid_id expected0[] = {0, 1, 2, 4, 13};
  validate(expected0, expected0 + 5, out_ids.begin());
  
  out_ids.clear();
  p.write_members(1, back_inserter(out_ids));
  medoid_id expected1[] = {3, 5, 6, 7, 9};
  validate(expected1, expected1 + 5, out_ids.begin());
  
  out_ids.clear();
  p.write_members(2, back_inserter(out_ids));
  medoid_id expected2[] = {8, 10, 11, 12};
  validate(expected2, expected2 + 4, out_ids.begin());


  ostringstream out;
  out << p.members(0);
  const char *expected = "0-2 4 13";
  validate(expected, expected + strlen(expected), out.str().c_str());
  
  out.str("");
  out << p.members(1);
  expected = "3 5-7 9";
  validate(expected, expected + strlen(expected), out.str().c_str());
  
  out.str("");
  out << p.members(2);
  expected = "8 10-12";
  validate(expected, expected + strlen(expected), out.str().c_str());
  

  if (passed) {
    cerr << "PASSED" << endl;
    exit(0);
  } else {
    cerr << "FAILED" << endl;
    exit(1);
  }
}
