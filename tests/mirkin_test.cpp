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
/// @file mirkin_test.cpp
/// @author Todd Gamblin tgamblin@llnl.gov
/// 
#include <algorithm>
#include <iterator>
#include <iostream>
#include <set>
#include <cstdlib>
#include <sstream>

#include "partition.h"

using namespace std;
using namespace cluster;

void makeClusterings(cluster_list& c1, cluster_list& c2) {
  const char *clusters1[] = {
    "0 2 4 6",
    "1 3",
    "5 7"
  };

  const char *clusters2[] = {
    "0 2 4 6",
    "1 3 7",
    "5"
  };
  
  c1.resize(sizeof(clusters1)/sizeof(char*));
  c2.resize(sizeof(clusters2)/sizeof(char*));

  for (size_t i=0; i < c1.size(); i++) {
    istringstream s1(clusters1[i]);
    istringstream s2(clusters2[i]);
    copy(istream_iterator<size_t>(s1), istream_iterator<size_t>(), inserter(c1[i], c1[i].begin()));
    copy(istream_iterator<size_t>(s2), istream_iterator<size_t>(), inserter(c2[i], c2[i].begin()));
  }
}


int main(int argc, char **argv) {
  cluster_list c1, c2;
  makeClusterings(c1, c2);

  cout << c1 << endl;
  cout << c2 << endl;
  cout << mirkin_distance(c1, c2) << endl;

  expand(c1);
  expand(c2);

  cout << c1 << endl;
  cout << c2 << endl;
  cout << mirkin_distance(c1, c2) << endl;

  expand(c1);
  expand(c2);

  cout << c1 << endl;
  cout << c2 << endl;
  cout << mirkin_distance(c1, c2) << endl;
}
