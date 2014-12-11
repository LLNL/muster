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
/// @file kmedoids.cpp
/// @author Todd Gamblin tgamblin@llnl.gov
///
#include "kmedoids.h"
using namespace cluster;

#include <vector>
#include <sstream>

#include <algorithm>
#include <numeric>
#include <cassert>
#include <cstdlib>
#include <sys/time.h>
using namespace std;

#include "random.h"
#include "bic.h"
#include "matrix_utils.h"

namespace cluster {
  
  kmedoids::kmedoids(size_t num_objects) 
    : partition(num_objects), 
      random(get_time_seed()),
      rng(random),
      total_dissimilarity(std::numeric_limits<double>::infinity()),
      sort_medoids(true),
      epsilon(1e-15),
      init_size(40),
      max_reps(5),
      xcallback(NULL)
  { }


  kmedoids::~kmedoids() {  }

  double kmedoids::average_dissimilarity() const {
    return total_dissimilarity / cluster_ids.size();
  }

  void kmedoids::set_seed(unsigned long s) {
    random.seed(s);
  }

  void kmedoids::set_sort_medoids(bool sort) {
    sort_medoids = sort;
  }

  void kmedoids::set_epsilon(double e) {
    epsilon = e;
  }
  
  void kmedoids::set_xcallback(void (*xpc)(const partition& part, double bic)) {
    xcallback = xpc;
  }

  void kmedoids::init_medoids(size_t k, const dissimilarity_matrix& distance) {
    medoid_ids.clear();
    // find first oject: object minimum dissimilarity to others
    object_id first_medoid = 0;
    double min_dissim = DBL_MAX;
    for (size_t i=0; i < distance.size1(); i++) {
      double total = 0.0;
      for (size_t j=0; j < distance.size2(); j++) {
        total += distance(i,j);
      }
      if (total < min_dissim) {
        min_dissim   = total;
        first_medoid = i;
      }
    }
    
    // add first object to medoids and compute medoid ids.
    medoid_ids.push_back(first_medoid);
    assign_objects_to_clusters(matrix_distance(distance));

    // now select next k-1 objects according to KR's BUILD algorithm
    for (size_t cur_k = 1; cur_k < k; cur_k++) {
      object_id best_obj = 0;
      double max_gain = 0.0;
      for (size_t i=0; i < distance.size1(); i++) {
        if (is_medoid(i)) continue;

        double gain = 0.0;
        for (size_t j=0; j < distance.size1(); j++) {
          double Dj = distance(j, medoid_ids[cluster_ids[j]]);  // distance from j to its medoid
          gain += max(Dj - distance(i,j), 0.0);                 // gain from selecting i  
        }

        if (gain >= max_gain) {   // set the next medoid to the object that 
          max_gain = gain;        // maximizes the gain function.
          best_obj = i;
        }
      }
      
      medoid_ids.push_back(best_obj);
      assign_objects_to_clusters(matrix_distance(distance));
    }
  }


  double kmedoids::cost(medoid_id i, object_id h, const dissimilarity_matrix& distance) const {
    double total = 0;
    for (object_id j = 0; j < cluster_ids.size(); j++) {
      object_id mi  = medoid_ids[i];                // object id of medoid i
      double    dhj = distance(h, j);               // distance between object h and object j
      
      object_id mj1 = medoid_ids[cluster_ids[j]];   // object id of j's nearest medoid
      double    dj1 = distance(mj1, j);             // distance to j's nearest medoid

      // check if distance bt/w medoid i and j is same as j's current nearest medoid.
      if (distance(mi, j) == dj1) {
        double dj2 = DBL_MAX;
        if (medoid_ids.size() > 1) {   // look at 2nd nearest if there's more than one medoid.
          object_id mj2 = medoid_ids[sec_nearest[j]];  // object id of j's 2nd-nearest medoid
          dj2 = distance(mj2, j);                      // distance to j's 2nd-nearest medoid
        }
        total += min(dj2, dhj) - dj1;

      } else if (dhj < dj1) {
        total += dhj - dj1;
      }
    }
    return total;
  }


  void kmedoids::pam(const dissimilarity_matrix& distance, size_t k, const object_id *initial_medoids) {
    if (k > distance.size1()) {
      throw std::logic_error("Attempt to run PAM with more clusters than data.");
    }

    if (distance.size1() != distance.size2()) {
      throw std::logic_error("Error: distance matrix is not square!");
    }
    
    // first get this the right size.
    cluster_ids.resize(distance.size1());

    // size cluster_ids appropriately and randomly pick initial medoids
    if (initial_medoids) {
      medoid_ids.clear();
      copy(initial_medoids, initial_medoids + k, back_inserter(medoid_ids));
    } else {
      init_medoids(k, distance);
    }

    // set tolerance equal to epsilon times mean magnitude of distances.
    // Note that distances *should* all be non-negative.
    double tolerance = epsilon * sum(distance) / (distance.size1() * distance.size2());

    while (true) {
      // initial cluster setup
      total_dissimilarity = assign_objects_to_clusters(matrix_distance(distance));

      //vars to keep track of minimum
      double minTotalCost = DBL_MAX;
      medoid_id minMedoid = 0;
      object_id minObject = 0;

      //iterate over each medoid
      for (medoid_id i=0; i < k; i++) {
        //iterate over all non-medoid objects
        for (object_id h = 0; h < cluster_ids.size(); h++) {
          if (is_medoid(h)) continue;

          //see if the total cost of swapping i & h was less than min
          double curCost = cost(i, h, distance);
          if (curCost < minTotalCost) {
            minTotalCost = curCost;
            minMedoid = i;
            minObject = h;
          }
        }
      }

      // bail if we can't gain anything more (we've converged)
      if (minTotalCost >= -tolerance) break;

      // install the new medoid if we found a beneficial swap
      medoid_ids[minMedoid] = minObject;
      cluster_ids[minObject] = minMedoid;
    }
    
    if (sort_medoids) sort();
  }


  double kmedoids::xpam(const dissimilarity_matrix& distance, size_t max_k, size_t dimensionality) {
    double best_bic = -DBL_MAX;   // note that DBL_MIN isn't what you think it is.

    for (size_t k = 1; k <= max_k; k++) {
      kmedoids subcall;
      subcall.pam(distance, k);
      double cur_bic = bic(subcall, matrix_distance(distance), dimensionality);

      if (xcallback) xcallback(subcall, cur_bic);

      if (cur_bic > best_bic) {
        best_bic = cur_bic;
        swap(subcall);
      }
    }
    return best_bic;
  }


} // namespace cluster  
