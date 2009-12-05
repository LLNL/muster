#include "kmedoids.h"
using namespace cluster;

#include <mpi.h>
#include <vector>
#include <fstream>

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <sys/time.h>
using namespace std;

#include "counter.h"
#include "matrix_utils.h"

// relative error for cost function.  Avoids roundoff error.
static const double epsilon = 1e-15;

namespace cluster {
  
  kmedoids::kmedoids(size_t num_objects) : partition(num_objects) {
    struct timeval time;
    gettimeofday(&time, 0);
    random.seed(time.tv_sec * time.tv_usec);
  }


  kmedoids::~kmedoids() {  }


  void kmedoids::init_medoids(size_t k) {
    medoids.clear();
    random_subset(cluster_ids.size(), k, back_inserter(medoids), random);
    assert(medoids.size() == k);
  }


  double kmedoids::cost(medoid_id i, object_id h, const dissimilarity_matrix& distance) const {
    static ofstream branches("branches");

    double total = 0;
    for (object_id j = 0; j < cluster_ids.size(); j++) {
      if (is_medoid(j) || j == h) continue;      //skip medoids and self

      object_id mi  = medoids[i];                // object id of medoid i
      double    dhj = distance(h, j);            // distance between object h and object j
      
      object_id mj1 = medoids[cluster_ids[j]];   // object id of j's nearest medoid
      double    dj1 = distance(mj1, j);          // distance to j's nearest medoid

      // check if distance bt/w medoid i and j is same as j's current nearest medoid.
      if (distance(mi, j) == dj1) {
        double dj2 = DBL_MAX;
        if (medoids.size() > 1) { // look at 2nd nearest if there's more than one medoid.
          object_id mj2 = medoids[sec_nearest[j]];  // object id of j's 2nd-nearest medoid
          dj2 = distance(mj2, j);                   // distance to j's 2nd-nearest medoid
        }
        total += min(dj2, dhj) - dj1;

      } else if (dhj < dj1) {
        total += dhj - dj1;
      }
    }
    return total;
  }


  void kmedoids::pam(dissimilarity_matrix distance, size_t k) {
    if (k > distance.size1()) {
      throw std::logic_error("Attempt to instantiate kmedoids with more clusters than data.");
    }

    if (distance.size1() != distance.size2()) {
      throw std::logic_error("Error: distance matrix is not square!");
    }
    
    // first get this the right size.
    cluster_ids.resize(distance.size1());

    // size cluster_ids appropriately and randomly pick initial medoids
    init_medoids(k);

    // set tolerance equal to epsilon times mean magnitude of distances.
    // Note that distances *should* all be non-negative.
    double tolerance = epsilon * sum(distance) / (distance.size1() * distance.size2());

    while (true) {
      // initial cluster setup
      average_dissimilarity = assign_objects_to_clusters(matrix_distance(distance));

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
      medoids[minMedoid] = minObject;
      cluster_ids[minObject] = minMedoid;
    }
  }


} // namespace cluster  
