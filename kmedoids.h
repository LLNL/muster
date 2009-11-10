
#ifndef K_MEDOIDS_H
#define K_MEDOIDS_H

#include <vector>
#include <set>
#include <iostream>
#include <stdexcept>
#include <cfloat>

#include "random.h"
#include "dissimilarity.h"
#include "MersenneTwister.h"
#include "partition.h"

namespace cluster {

  /// 
  /// Implementation of classical clustering methods PAM and CLARA, from 
  /// "Finding Groups in Data", by Kaufman and Rousseeuw.  
  /// 
  class kmedoids : public partition {
  public:
    ///
    /// Constructor.  Can optionally specify number of objects to be clustered.
    /// and this will start out with all of them in one cluster.
    /// 
    /// The random number generator associated with this kmedoids object is seeded
    /// with the time in microseconds since the epoch.
    /// 
    kmedoids(size_t num_objects = 0);
    ~kmedoids();

    /// Classic K-Medoids clustering, using the Partitioning-Around-Medoids (PAM)
    /// algorithm as described in Kaufman and Rousseeuw. 
    /// Parameters:
    ///   distance     dissimilarity matrix for all objects to cluster
    ///   k            number of clusters to produce
    void pam(dissimilarity_matrix distance, size_t k);
    
    ///
    /// CLARA clustering algorithm, as per Kaufman and Rousseuw and
    /// R. Ng and J. Han, "Efficient and Effective Clustering Methods 
    /// for Spatial Data Mining."
    /// 
    /// Template parameters (inferred from args):
    ///   T              Type of objects to be clustered.
    ///   D              Dissimilarity metric type.  D should be callable 
    ///                  on (T, T) and should return a double.
    /// 
    /// Parameters:
    ///   objects        Objects to cluster
    ///   dmetric        Distance metric to build dissimilarity matrices with
    ///   k              Number of clusters to partition
    ///   sample_size    defaults to 40+2*k, per Kaufman and Rousseeuw's recommendation
    ///   iterations     Number of times to run PAM with sampled dataset
    ///
    template <class T, class D>
    void clara(const std::vector<T> objects, D dmetric,
               size_t k, size_t init_size = 40, size_t iterations=5) {

      size_t sample_size = init_size + 2*k;
    
      // Just run plain KMedoids once if sampling won't gain us anything
      if (objects.size() <= sample_size) {
        dissimilarity_matrix mat;
        build_dissimilarity_matrix(objects, dmetric, mat);
        pam(mat, k);
        return;
      }

      // get everything the right size before starting.
      medoids.resize(k);
      cluster_ids.resize(objects.size());

      // medoids and clusters for best partition so far.
      partition best_partition;

      //run KMedoids on a sampled subset ITERATIONS times
      double best_dissim = DBL_MAX;
      for (size_t i = 0; i < iterations; i++) {
        // Take a random sample of objects, store sample in a vector
        std::vector<size_t> sample_to_full;
        random_subset(objects.size(), sample_size, back_inserter(sample_to_full), random);

        // Build a distance matrix for PAM
        dissimilarity_matrix distance;
        build_dissimilarity_matrix(objects, sample_to_full, dmetric, distance);

        // Actually run PAM on the subset
        kmedoids subcall;
        subcall.pam(distance, k);

        // copy medoids from the subcall to local data, being sure to translate indices
        for (size_t i=0; i < medoids.size(); i++) {
          medoids[i] = sample_to_full[subcall.medoids[i]];
        }

        // sync up the cluster_ids matrix with the new medoids by assigning
        // each object to its closest medoid.  Remember the quality of the clustering.
        average_dissimilarity = assign_objects_to_clusters(lazy_distance(objects, dmetric));

        // keep the best clustering found so far around
        if (average_dissimilarity < best_dissim) {
          this->swap(best_partition);
          best_dissim = average_dissimilarity;
        } 
      }
      
      this->swap(best_partition);
      average_dissimilarity = best_dissim;
    }    

    protected:
    MTRand random;             /// Random number generator for this algorithm
    double average_dissimilarity;   /// Avg dissimilarity for last clustering run.

    /// Assigns medoids randomly from the input objects.
    void init_medoids(size_t k);

    /// Finds the cost of swapping object oh with medoid j.
    /// PRE: object[oi] is the medoid with index cluster_id[oi] in medoids.
    double cost(medoid_id mi, object_id oh, object_id oj, const dissimilarity_matrix& distance) const;

    /// Total cost of swapping object h with medoid i.
    /// Sums costs of this exchagne for all objects j.
    double total_cost(medoid_id i, object_id h, const dissimilarity_matrix& distance) const;

    /// Assign each object to the cluster with the closest medoid.
    /// Returns:
    ///   Average dissimilarity of objects w/their medoids.
    /// 
    /// D should be a callable object that computes distances.  
    ///   In PAM, it should use the distance matrix.  
    ///   In CLARA, it should compute lazily.
    template <class D>
    double assign_objects_to_clusters(D distance) {
      // first point the medoid objects at the right clusters
      for (medoid_id mi = 0; mi < medoids.size(); mi++) {
        cluster_ids[medoids[mi]] = mi;
      }

      // now go through and assign each object to nearest medoid, keeping track of total 
      // dissimilarity, too.
      double sum = 0;
      for (object_id i=0; i < cluster_ids.size(); i++) {
        if (!is_medoid(i)) {
          std::pair<medoid_id, double> nearest = nearest_medoid(i, distance);
          cluster_ids[i] = nearest.first;
          sum += nearest.second;
        }
      }

      return sum / cluster_ids.size();
    }

    /// Finds the nearest medoid to the object oi.
    /// Returns:
    ///   - Index of the nearest medoid to oi
    ///   - distance from oi to this medoid.
    ///   
    /// D should be a callable object that computes distances.  
    ///   In PAM, it should use the distance matrix.  
    ///   In CLARA, it should compute lazily.
    template <class D>
    std::pair<medoid_id, double> nearest_medoid(
      object_id oi, D distance, medoid_id exclude=std::numeric_limits<size_t>::max()
    ) const {
      medoid_id nearest = 0;
      double min_distance = DBL_MAX;
      
      for (size_t i=0; i < medoids.size(); i++) {
        if (i == exclude) continue;
        double d = distance(oi, medoids[i]);
        if (d < min_distance) {
          min_distance = d;
          nearest = i;
        }
      }
      return std::make_pair(nearest, min_distance);
    }


    /// Functor for computing distance lazily from an object array and
    /// a distance metric.  Use this for CLARA, where we don't want to
    /// precompute the entire distance matrix.
    template <class T, class D>
    struct lazy_distance_functor {
      const std::vector<T>& objects;
      D dissimilarity;

      lazy_distance_functor(const std::vector<T>& objs, D d)
        : objects(objs), dissimilarity(d) { }

      double operator()(size_t i, size_t j) {
        return dissimilarity(objects[i], objects[j]);
      }
    };

    /// Type-inferred syntactic sugar for constructing lazy_distance_functors.
    template <class T, class D>
    lazy_distance_functor<T,D> lazy_distance(const std::vector<T>& objs, D dist) {
      return lazy_distance_functor<T,D>(objs, dist);
    }

    /// Simple adapter for passing a matrix by reference to the template
    /// functions above.  Avoids wholesale copy of distance matrix.
    struct matrix_distance {
      const dissimilarity_matrix& mat;
      matrix_distance(const dissimilarity_matrix& m) : mat(m) { }
      double operator()(size_t i, size_t j) { return mat(i,j); }
    };
  };


} // namespace cluster

#endif //K_MEDOIDS_H
