#ifndef PAR_KMEDOIDS_H
#define PAR_KMEDOIDS_H

#include <mpi.h>
#include <vector>
#include <sstream>
#include <fstream>
#include <set>

#include "Timer.h"
#include "kmedoids.h"
#include "multi_gather.h"
#include "mpi_utils.h"
#include "par_partition.h"
#include "MersenneTwister.h"

namespace cluster {

  ///
  /// This struct represents parameters for a single trial run of kmedoids.
  /// We generate a bunch of these to farm all the trials out to processes in par_clara.
  ///
  struct trial {
    size_t k;
    size_t trial;
    size_t sample_size;
  };
  
  ///
  /// This class packages iterates through all k's and trial5B numbers.  It basically packages
  /// up loop state to clean up the dispatch of clustering jobs to worker processes.
  ///
  class trial_iterator {
    const size_t max_k;         // maximum k to try
    const size_t max_trials;    // max number of trials per k
    const size_t init_size;     // initial size for samples before factoring in k, as per CLARA paper.
    const size_t max_sample;    // maximum sample size (usually number of elements in the data set)
    trial cur_trial;      // current state of the iterator.
    
    // size of the sample to cluster for particular k
    size_t get_sample_size(size_t k);
    
  public:
    trial_iterator(size_t _max_k, size_t _max_trials, size_t _init_size, size_t _max_sample);

    trial_iterator(size_t min_k, size_t _max_k, 
                   size_t _max_trials, size_t _init_size, size_t _max_sample);
    
    bool has_next();
    trial next();
  };
  

  ///
  /// Packable struct for a packable type plus its source rank.
  ///
  template <class T>
  struct source_pair {
    T element;
    int source;

    source_pair() { }
    source_pair(const T& elt, int s) : element(elt), source(s) { }

    int packed_size(MPI_Comm comm) const {
      return element.packed_size(comm) + mpi_packed_size(1, MPI_INT, comm);
    }

    void pack(void *buf, int bufsize, int *position, MPI_Comm comm) const {
      element.pack(buf, bufsize, position, comm);
      MPI_Pack(const_cast<int*>(&source), 1, MPI_INT, buf, bufsize, position, comm);
    }

    static source_pair unpack(void *buf, int bufsize, int *position, MPI_Comm comm) {
      T t = T::unpack(buf, bufsize, position, comm);
      int source;
      MPI_Unpack(buf, bufsize, position, &source, 1, MPI_INT, comm);
      return source_pair(t, source);
    }
  };

  
  /// Helper function for making source_pairs with type inference.
  template <class T>
  source_pair<T> make_source_pair(const T& elt, int source) {
    return source_pair<T>(elt, source);
  }
  
  
  class par_kmedoids : public par_partition {
  public:
    ///
    /// Constructs a parallel kmedoids object and seeds its random number generator.
    /// This is a collective operation, and needs to be called by all processes.
    ///
    /// par_kmedoids assumes that there is one object to be clustered per process.
    ///
    par_kmedoids() : par_partition() { }

    
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
    ///   medoids        Output vector where global medoids will be stored along with their source ranks.
    ///   max_trials     Max number of times to run PAM with each sampled dataset.
    ///   sample_size    defaults to 40+2*k, per Kaufman and Rousseeuw's recommendation
    ///
    /// The parallel version of CLARA will run [1..max_k] * trials instances of PAM on sets
    /// of sample_size objects distributed over all processes in the system.
    ///
    template <class T, class D>
    void par_clara(const T& object, D dmetric, size_t max_k, std::vector< source_pair<T> >& medoids,
                   size_t max_trials=5, size_t init_size=40)
    {
      Timer timer;
      
      int size, rank;
      MPI_Comm_size(comm, &size);
      MPI_Comm_rank(comm, &rank);

      seed_random_uniform(comm); // seed RN generator uniformly across ranks.

      // fix things if k is greater than the number of elements, since we can't 
      // ever find that many clusters.
      if (max_k > (size_t)size) max_k = size;

      std::vector< source_pair<T> > all_medoids[max_k][max_trials];

      timer.record("Init");
      
      trial_iterator trials(max_k, max_trials, init_size, size);
      while (trials.has_next()) {
        int my_k = -1;                      // k for local run of kmedoids
        int my_trial = -1;                  // trial id for local run of kmedoids
        std::vector<int>  my_samples;       // source ranks for each of my_objects
        std::vector<T>    my_objects;       // vector to hold local sample objects for clustering.
        multi_gather<T>   gather(comm);     // simultaneous, asynchronous local gathers for collecting samples.

        // save the iteration state so we can run through it again when we bcast the results.
        trial_iterator last_trials = trials;   

        // start gathers for each trial to aggregate samples to single worker processes.
        for (int root =0; trials.has_next() && root < size; root++) {
          trial cur_trial = trials.next();
          
          // Generate a set of indices for members of this k-medoids trial
          std::vector<int> samples;
          random_subset(size, cur_trial.sample_size, back_inserter(samples), random);
        
//           std::ostringstream msg;
//           msg << rank << ": ";
//           for (size_t i = 0; i < samples.size(); i++) {
//             msg << samples[i] << " ";
//           }
//           msg << std::endl;
//           std::cerr << msg.str();
        
          gather.start(object, samples.begin(), samples.end(), my_objects, root);
          if (rank == root) {
            my_k       = cur_trial.k;
            my_trial   = cur_trial.trial;
            samples.swap(my_samples);
          }
        }

        timer.record("StartGather");
        
        // finish all sample gathers.
        gather.finish();
        timer.record("FinishGather");
        
        // if we're a worker process (we were assigned a k and a trial number) 
        // then run PAM on the sample.
        kmedoids cluster;
        if (my_k >= 0) {
          dissimilarity_matrix mat;
          build_dissimilarity_matrix(my_objects, dmetric, mat);
          cluster.pam(mat, my_k);

          // put the local medoids we computed into the global medoids array.
          for (size_t i=0; i < cluster.medoids.size(); i++) {
            // push the medoid, along with its source rank, onto the medoids vector.
            all_medoids[my_k-1][my_trial].push_back(
              make_source_pair(my_objects[cluster.medoids[i]], 
                               my_samples[cluster.medoids[i]]));
          }
        }
        timer.record("LocalCluster");

        // once workers are done clustering, broadcast the medoids from each clustering 
        // so that all processes can figure out which cluster they're in.
        for (int root =0; last_trials.has_next() && root < size; root++) {
          trial cur_trial = last_trials.next();
          all_medoids[cur_trial.k-1][cur_trial.trial].resize(cur_trial.k);
          bcast_medoids(comm, all_medoids[cur_trial.k-1][cur_trial.trial], root);
        }

        timer.record("Broadcast");
      }
      

      std::ostringstream name;
      name << "times." << rank;
      std::ofstream file(name.str().c_str());
      timer.dump(file);
      
      
      /*
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

  protected:
    /// Random number generator, seeded the same on each process
    MTRand random;

    /// 
    /// Seeds random number generators across all processes with the same number,
    /// taken from the time in microseconds since the epoch on the process 0.
    /// 
    void seed_random_uniform(MPI_Comm comm);

    ///
    ///
    ///
    template <class T>
    void bcast_medoids(MPI_Comm comm, std::vector<T>& medoids, int root) {
      int rank;
      MPI_Comm_rank(comm, &rank);
      
      int packed_size=0;
      if (rank == root) {
        // figure out size of packed buffer
        for (size_t i=0; i < medoids.size(); i++) {
          packed_size += medoids[i].packed_size(comm);
        }
      }

      // broadcast size and allocate receive buffer.
      MPI_Bcast(&packed_size, 1, MPI_INT, root, comm);
      char packed_buffer[packed_size];

      if (rank == root) {
        // pack buffer with medoid objects
        for (size_t i=0; i < medoids.size(); i++) {
          int pos = 0;
          medoids[i].pack(packed_buffer, packed_size, &pos, comm);
        }
      }

      MPI_Bcast(&packed_buffer, packed_size, MPI_PACKED, root, comm);

      // unpack broadcast medoids.
      if (rank != root) {
        int pos = 0;
        for (size_t i=0; i < medoids.size(); i++) {
          medoids[i] = T::unpack(packed_buffer, packed_size, &pos, comm);
        }
      }
    }
  };

} // namespace cluster

#endif // PAR_KMEDOIDS_H
