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
/// @file par_kmedoids.cpp
/// @author Todd Gamblin tgamblin@llnl.gov
///
#include "par_kmedoids.h"

#include <cstdlib>
#include <stdint.h>
#include <sys/time.h>
using namespace std;

#include "random.h"

namespace cluster {

  par_kmedoids::par_kmedoids(MPI_Comm comm) 
    : par_partition(comm),
      seed_set(false),
      total_dissimilarity(numeric_limits<double>::infinity()),
      best_bic_score(0),
      init_size(40),
      max_reps(5),
      epsilon(1e-15)
  { }

  void par_kmedoids::set_seed(uint32_t s) {
      random.seed(s);
      seed_set = true;
  }

  void par_kmedoids::set_epsilon(double e) {
    epsilon = e;
  }


  double par_kmedoids::average_dissimilarity() {
    int size;
    CMPI_Comm_size(comm, &size);
    return total_dissimilarity / size;
  }

  double par_kmedoids::bic_score() {
    return best_bic_score;
  }

  void par_kmedoids::seed_random_uniform(MPI_Comm comm) {
    int rank;
    CMPI_Comm_rank(comm, &rank);

    // same seed on all processes.
    uint32_t seed = get_time_seed();
    CMPI_Bcast(&seed, 1, MPI_INT, 0, comm);
    random.seed(seed);
    seed_set = true;
  }

} // namespace cluster
