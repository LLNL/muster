#ifndef PAR_KMEDOIDS_H
#define PAR_KMEDOIDS_H

#include "kmedoids.h"
#include <mpi.h>

namespace cluster {

  class par_kmedoids : public kmedoids {
    
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
    void par_clara(T object, D dmetric, size_t max_k, MPI_Comm comm,
                   size_t sample_size = 0, size_t iterations=5) {

      if (!sample_size) sample_size = 40+2*k;
      
      int size, rank;
      MPI_Comm_size(comm, &size);
      MPI_Comm_rank(comm, &rank);

      // Generate all samples using the same sequence of random numbers across all processes.
      vector<int> samples;
      random_subset(size, sample_size * iterations, back_inserter(sample), random);

      // Aggregate samples to worker processes to do PAM clustering.
      // Attempt to distribute workers evenly through ranks by using mod sets.
      // Iterate multiple times through this if there are not enough workers to 
      // do everything in parallel.
      std::vector<T> my_objects;
      
      
      
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


    ///
    ///   Gathers T to rank <root> from all ranks in the source_ranks vector.
    ///   Results will be stored in the dest vector on the root, but the dest
    ///   vector is untouched on other processes.
    ///   
    ///   Communication is done asynchronously, and requests made are appended to 
    ///   the reqs vector.
    ///   
    ///   asynchronously.  Requests (for later completion with MPI_Waitall() or equivalent)
    ///   are pl
    ///   
    ///   T must support the following operations:
    ///     - int packed_size(MPI_Comm comm) const
    ///     - void pack(void *buf, int bufsize, int *position, MPI_Comm comm) const
    ///     - static T unpack(void *buf, int bufsize, int *position, MPI_Comm comm)
    ///
    template <class Iterator, class T>
    void gather(const T& my_object, Iterator start_rank, Iterator end_rank, std::vector<T> dest, 
                   MPI_Comm comm, std::vector<MPI_Request> reqs, int root=0) 
    {
      int size, rank;
      MPI_Comm_size(comm, &size);
      MPI_Comm_rank(comm, &rank);
      
      if (rank == root) {
        dest.resize(end_rank - start_rank);
        for (Iterator i=start_rank; i != end_rank; i++) {
          MPI_Irecv();
        }
      } else {
        
      }
    }
    


        
  };

} // namespace cluster

#endif // PAR_KMEDOIDS_H
