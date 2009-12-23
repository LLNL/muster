#include "kmedoids.h"
using namespace cluster;

#include <mpi.h>
#include <vector>
#include <sstream>

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <sys/time.h>
using namespace std;

#include "counter.h"
#include "matrix_utils.h"

#include "point.h"

// relative error for cost function.  Avoids roundoff error.
static const double epsilon = 1e-15;

namespace cluster {
  
  kmedoids::kmedoids(size_t num_objects) 
    : partition(num_objects), 
      total_dissimilarity(std::numeric_limits<double>::infinity()),
      sort_medoids(true)
  {
    struct timeval time;
    gettimeofday(&time, 0);
    random.seed(time.tv_sec * time.tv_usec);
  }


  kmedoids::~kmedoids() {  }

  double kmedoids::average_dissimilarity() {
    return total_dissimilarity / cluster_ids.size();
  }

  void kmedoids::set_sort_medoids(bool sort) {
    sort_medoids = sort;
  }


  void kmedoids::init_medoids(size_t k) {
    medoid_ids.clear();
    random_subset(cluster_ids.size(), k, back_inserter(medoid_ids), random);
    assert(medoid_ids.size() == k);

    ostringstream msg;
    msg << "initial medoids: [";
    copy(medoid_ids.begin(), medoid_ids.end(), ostream_iterator<object_id>(msg, " "));
    msg << "]" << endl;
    cerr << msg.str();

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


  void kmedoids::pam(dissimilarity_matrix distance, size_t k, const object_id *initial_medoids) {
    if (k > distance.size1()) {
      throw std::logic_error("Attempt to instantiate kmedoids with more clusters than data.");
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
      init_medoids(k);
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


} // namespace cluster  
