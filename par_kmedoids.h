#ifndef PAR_KMEDOIDS_H
#define PAR_KMEDOIDS_H

#include <mpi.h>
#include <vector>
#include <set>

#include "kmedoids.h"
#include "multi_gather.h"

namespace cluster {

  class par_kmedoids : public kmedoids {
  public:
    ///
    /// Constructs a parallel kmedoids object and seeds its random number generator.
    /// This is a collective operation, and needs to be called by all processes.
    ///
    /// par_kmedoids assumes that there is one object to be clustered per process.
    ///
    par_kmedoids() : kmedoids() { }

    
    ///
    /// Parallel version of the CLARA clustering algorithm.  Assumes that objects
    /// to be clustered are fully distributed across parallel process, one object
    /// per process.  
    ///
    /// Template parameters (inferred from args):
    ///   T              Type of objects to be clustered.
    ///                  T must also support the following operations:
    ///                    - int packed_size(MPI_Comm comm) const
    ///                    - void pack(void *buf, int bufsize, int *position, MPI_Comm comm) const
    ///                    - static T unpack(void *buf, int bufsize, int *position, MPI_Comm comm)
    ///   
    ///   D              Dissimilarity metric type.  D should be callable 
    ///                  on (T, T) and should return a double.
    /// 
    /// Parameters:
    ///   objects        Objects to cluster
    ///   dmetric        Distance metric to build dissimilarity matrices with
    ///   max_k          Max number of clusters to find.
    ///   sample_size    defaults to 40+2*k, per Kaufman and Rousseeuw's recommendation
    ///   iterations     Number of times to run PAM with sampled dataset
    ///
    /// The parallel version of CLARA will run [2..max_k] * iterations instances of PAM on sets
    /// of sample_size objects distributed over all processes in the system.
    ///
    template <class T, class D>
    void par_clara(const T& object, D dmetric, size_t max_k, MPI_Comm comm,
                   size_t init_size = 40, size_t iterations=5) {
      
      int size, rank;
      MPI_Comm_size(comm, &size);
      MPI_Comm_rank(comm, &rank);

      seed_random(comm); // seed RN generator uniformly across ranks.

      // Aggregate samples to worker processes to do PAM clustering.
      std::vector<T>  my_objects;
      multi_gather<T> gather(comm);

      int root = 0;
      for (size_t k = 1; k < max_k; k++) {
        for (size_t trial = 0; trial < 5; trial++) {
          // size of the sample to cluster.
          size_t sample_size = init_size + 2 * k;
          
          // Generate a set of indices for members of this k-medoids trial
          std::set<int> samples;
          random_subset(size, sample_size, inserter(samples, samples.begin()), random);

          // if this rank is either the root or if it's a sampled process, start a gather.
          if (samples.count(rank) || rank == root) {
            gather.start(object, samples.begin(), samples.end(), my_objects, root);
            root++;
          }
        }
      }
      
      
      


      /*
      
      
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
      */
    }    

    /// 
    /// Seeds random number generators across all processes with the same number,
    /// taken from the time in microseconds since the epoch on the process 0.
    /// 
    void seed_random(MPI_Comm comm);
    
  protected:




  };

} // namespace cluster

#endif // PAR_KMEDOIDS_H
