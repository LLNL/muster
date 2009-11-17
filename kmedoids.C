#include "kmedoids.h"
using namespace cluster;

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <sys/time.h>
using namespace std;

#include "counter.h"


namespace cluster {

  
  kmedoids::kmedoids(size_t num_objects) : partition(num_objects) {
    struct timeval time;
    gettimeofday(&time, 0);
    random.seed(time.tv_sec * time.tv_usec);
  }


  kmedoids::~kmedoids() {
    //nothing necessary.
  }


  void kmedoids::init_medoids(size_t k) {
    medoids.clear();
    random_subset(cluster_ids.size(), k, back_inserter(medoids), random);
    assert(medoids.size() == k);
  }


  /// Finds the cost of swapping object oh with medoid j.
  /// PRE: object[oi] is the medoid with index cluster_id[oi] in medoids.
  double kmedoids::cost(medoid_id mi, object_id oh, object_id oj, 
                        const dissimilarity_matrix& distance) const {
    //oi is the object id of medoid i, mi is the medoid id.
    object_id oi = medoids[mi];
  
    double cost;
    if (cluster_ids[oj] == mi) {
      medoid_id mj2 = nearest_medoid(oj, matrix_distance(distance), mi).first;
      object_id oj2 = medoids[mj2];
      if (distance(oj, oj2) < distance(oj, oh)) {
        cost = distance(oj, oj2) - distance(oj, oi);
      } else {
        cost = distance(oj, oh) - distance(oj, oi);            
      }
    
    } else {
      medoid_id mj2 = cluster_ids[oj];
      object_id oj2 = medoids[mj2];
      if (distance(oj, oj2) < distance(oj, oh)) {
        cost = 0.0;
      } else {
        cost = distance(oj, oh) - distance(oj, oj2);
      }
    }
  
    return cost;
  }


  /// Total cost of swapping object h with medoid i.
  /// Sums costs of this exchagne for all objects j.
  double kmedoids::total_cost(medoid_id i, object_id h, const dissimilarity_matrix& distance) const {
    double sum =0;
    for (object_id j = 0; j < cluster_ids.size(); j++) {
      if (is_medoid(j) || j == h) continue;   //skip medoids and self
      sum += cost(i, h, j, distance);    //add cost of swapping object h with object j
    }
    return sum;
  }


  void kmedoids::pam(dissimilarity_matrix distance, size_t k) {
    if (k > distance.size1()) {
      throw std::logic_error("Attempt to instantiate kmedoids with more clusters than data.");
    }

    if (distance.size1() != distance.size2()) {
      throw std::logic_error("Error: distance matrix is now square!");
    }
    
    // first get this the right size.
    cluster_ids.resize(distance.size1());

    // size cluster_ids appropriately and randomly pick initial medoids
    init_medoids(k);

    // initial cluster setup
    average_dissimilarity = assign_objects_to_clusters(matrix_distance(distance));
    if (k == 1) return;  // bail here if we only need one cluster.

    size_t counter = 0;
    while (true) {
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
          double curCost = total_cost(i, h, distance);
          if (curCost < minTotalCost) {
            minTotalCost = curCost;
            minMedoid = i;
            minObject = h;
          }
        }
      }

      // bail if we can't gain anything more (we've converged)
      if (minTotalCost >= 0 || counter++ > 10000) break;  // TODO: hack!  investigate convergence.

      //replace a medoid if it gains us something
      medoids[minMedoid] = minObject;
      cluster_ids[minObject] = minMedoid;

      //put objects in cluster w/nearest medoid
      average_dissimilarity = assign_objects_to_clusters(matrix_distance(distance));
    }
  }


} // namespace cluster  
