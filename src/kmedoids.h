#ifndef K_MEDOIDS_H
#define K_MEDOIDS_H

#include <vector>
#include <set>
#include <iostream>
#include <stdexcept>
#include <cfloat>

#include <boost/random.hpp>

#include "random.h"
#include "dissimilarity.h"
#include "partition.h"
#include "bic.h"

namespace cluster {

  /// 
  /// Implementations of the classic clustering algorithms PAM and CLARA, from 
  /// <i>Finding Groups in Data</i>, by Kaufman and Rousseeuw.
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

    /// Get the average dissimilarity of objects w/their medoids for the last run.
    double average_dissimilarity() const;
    
    /// Set whether medoids will be sorted by object id after clustering is complete.
    /// Defaults to true.
    void set_sort_medoids(bool sort_medoids);

    /// Set tolerance for convergence.  This is relative error, not absolute error.  It will be
    /// multiplied by the mean distance before it is used to test convergence.
    /// Defaults to 1e-15; may need to be higher if there exist clusterings with very similar quality.
    void set_epsilon(double epsilon);

    /// 
    /// Classic K-Medoids clustering, using the Partitioning-Around-Medoids (PAM)
    /// algorithm as described in Kaufman and Rousseeuw. 
    /// Parameters:
    ///   distance         dissimilarity matrix for all objects to cluster
    ///   k                number of clusters to produce
    ///   initial_medoids  Optionally supply k initial object ids to be used as initial medoids.
    /// 
    void pam(const dissimilarity_matrix& distance, size_t k, const object_id *initial_medoids = NULL);

    ///
    /// Classic K-Medoids clustering, using the Partitioning-Around-Medoids (PAM)
    /// algorithm as described in Kaufman and Rousseeuw. Runs PAM from 1 to max_k and selects
    /// the best k using the bayesian information criterion.  Sets this partition to the best
    /// partition found using PAM from 1 to k.
    /// 
    /// Based on X-Means, see Pelleg & Moore, 2000.
    /// 
    /// Parameters:
    ///   distance         dissimilarity matrix for all objects to cluster
    ///   max_k            Upper limit on number of clusters to find.
    ///   dimensionality   Number of dimensions in clustered data, for BIC.
    ///
    /// Return value:
    ///   This routine returns the best BIC value found (the bic value of the final partitioning).
    ///
    double xpam(const dissimilarity_matrix& distance, size_t max_k, size_t dimensionality);

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
    void clara(const std::vector<T>& objects, D dmetric, size_t k) {
      size_t sample_size = init_size + 2*k;
    
      // Just run plain KMedoids once if sampling won't gain us anything
      if (objects.size() <= sample_size) {
        dissimilarity_matrix mat;
        build_dissimilarity_matrix(objects, dmetric, mat);
        pam(mat, k);
        return;
      }

      // get everything the right size before starting.
      medoid_ids.resize(k);
      cluster_ids.resize(objects.size());

      // medoids and clusters for best partition so far.
      partition best_partition;

      //run KMedoids on a sampled subset max_reps times
      total_dissimilarity = DBL_MAX;
      for (size_t i = 0; i < max_reps; i++) {
        // Take a random sample of objects, store sample in a vector
        std::vector<size_t> sample_to_full;
        random_subset(objects.size(), sample_size, back_inserter(sample_to_full), rng);

        // Build a distance matrix for PAM
        dissimilarity_matrix distance;
        build_dissimilarity_matrix(objects, sample_to_full, dmetric, distance);

        // Actually run PAM on the subset
        kmedoids subcall;
        subcall.set_sort_medoids(false); // skip sort for subcall since it's not needed
        subcall.pam(distance, k);  

        // copy medoids from the subcall to local data, being sure to translate indices
        for (size_t i=0; i < medoid_ids.size(); i++) {
          medoid_ids[i] = sample_to_full[subcall.medoid_ids[i]];
        }

        // sync up the cluster_ids matrix with the new medoids by assigning
        // each object to its closest medoid.  Remember the quality of the clustering.
        double dissimilarity = assign_objects_to_clusters(lazy_distance(objects, dmetric));
        
        // keep the best clustering found so far around
        if (dissimilarity < total_dissimilarity) {
          swap(best_partition);
          total_dissimilarity = dissimilarity;
        } 
      }
      
      if (sort_medoids) sort();   // just do one final ordering of ids.
    }    


    ///
    /// Takes existing clustering and reassigns medoids by taking closest medoid to mean
    /// of each cluster.  This is O(n) and can give better representatives for CLARA clusterings.
    /// This is needed to apply gaussian-model BIC criteria to clusterings.
    /// 
    /// This function requires that T support algebraic operations.
    /// Specifically, T must support enough to construct a mean:
    ///    - addition        T + T = T
    ///    - scalar division T / c = T
    ///
    template <class T, class D>
    void center_medoids(const std::vector<T>& objects, D distance) {
      // First sum up elements in all clusters to get the mean for each
      std::vector<T>      means(medoid_ids.size());
      std::vector<size_t> counts(medoid_ids.size());
      for (size_t i=0; i < cluster_ids.size(); i++) {
        medoid_id m = cluster_ids[i];
        means[m] = counts[m] ? (means[m] + objects[i]) : objects[i];
        counts[m]++;
      }
      
      // Now, turn cumulative sums into means and calculate distance from medoids to means
      std::vector<double> shortest(means.size());
      for (size_t m=0; m < means.size(); m++) {
        means[m] = means[m] / counts[m];
        shortest[m] = distance(means[m], objects[medoid_ids[m]]);
      }
      
      // Find closest medoids to each mean element, preferring existing medoids if there are ties.
      for (size_t i=0; i < cluster_ids.size(); i++) {
        medoid_id m = cluster_ids[i];
        double d = distance(objects[i], means[m]);
        if (d < shortest[m]) {
          medoid_ids[m] = i;
          shortest[m]   = d;
        }
      }
    }


    ///
    /// TODO: figure out better BIC criterion for this.
    ///
    template <class T, class D>
    double xclara(const std::vector<T>& objects, D dmetric, size_t max_k, size_t dimensionality) {
      double best_bic = -DBL_MAX;   // note that DBL_MIN isn't what you think it is.

      for (size_t k = 1; k <= max_k; k++) {
        kmedoids subcall;
        subcall.clara(objects, dmetric, k);
        center_medoids(objects, dmetric);
        double cur_bic = bic(subcall, lazy_distance(objects, dmetric), dimensionality);

        if (xcallback) xcallback(subcall, cur_bic);
        if (cur_bic > best_bic) {
          best_bic = cur_bic;
          swap(subcall);
        }
      }
      return best_bic;
    }


    void set_init_size(size_t sz) { init_size = sz; }
    void set_max_reps(size_t r) { max_reps = r; }


    /// Set callback function for XPAM and XCLARA.  default is none.
    void set_xcallback(void (*)(const partition& part, double bic));

    protected:
    typedef boost::mt19937 random_type;                /// Type for RNG used in this algorithm
    random_type random;                                /// Randomness source for this algorithm
    
    /// Adaptor for STL algorithms.
    typedef boost::random_number_generator<random_type, unsigned long> rng_type;
    rng_type rng;

    std::vector<medoid_id> sec_nearest;      /// Index of second closest medoids.  Used by PAM.
    double total_dissimilarity;              /// Total dissimilarity bt/w objects and their medoid
    bool sort_medoids;                       /// Whether medoids should be canonically sorted by object id.
    double epsilon;                          /// Normalized sensitivity for convergence
    size_t init_size;                        /// initial sample size (before 2*k)
    size_t max_reps;                         /// initial sample size (before 2*k)


    /// Callback for each iteration of xpam.  is called with the current clustering and its BIC score.
    void (*xcallback)(const partition& part, double bic);

    /// KR BUILD algorithm for assigning initial medoids to a partition.
    void init_medoids(size_t k, const dissimilarity_matrix& distance);

    /// Total cost of swapping object h with medoid i.
    /// Sums costs of this exchagne for all objects j.
    double cost(medoid_id i, object_id h, const dissimilarity_matrix& distance) const;


    /// Assign each object to the cluster with the closest medoid.
    /// Returns:
    ///   Total dissimilarity of objects w/their medoids.
    /// 
    /// DM should be a callable object that computes distances between indices, as a distance 
    /// matrix would.  Algorithms are free to use real distance matrices (as in PAM) or to compute
    /// lazily (as in CLARA medoid assignment).
    template <class DM>
    double assign_objects_to_clusters(DM distance) {
      if (sec_nearest.size() != cluster_ids.size()) {
        sec_nearest.resize(cluster_ids.size());
      }
      
      // go through and assign each object to nearest medoid, keeping track of total dissimilarity.
      double total_dissimilarity = 0;
      for (object_id i=0; i < cluster_ids.size(); i++) {
        double    d1, d2;  // smallest, second smallest distance to medoid, respectively
        medoid_id m1, m2;  // index of medoids with distances d1, d2 from object i, respectively

        d1 = d2 = DBL_MAX;
        m1 = m2 = medoid_ids.size();
        for (medoid_id m=0; m < medoid_ids.size(); m++) {
          double d = distance(i, medoid_ids[m]);
          if (d < d1 || medoid_ids[m] == i) {  // prefer the medoid in case of ties.
            d2 = d1;  m2 = m1;
            d1 = d;   m1 = m;
          } else if (d < d2) {
            d2 = d;   m2 = m;
          }
        }

        cluster_ids[i] = m1;
        sec_nearest[i] = m2;
        total_dissimilarity += d1;
      }

      return total_dissimilarity;
    }
  };


} // namespace cluster

#endif //K_MEDOIDS_H
