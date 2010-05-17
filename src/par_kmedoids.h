#ifndef PAR_KMEDOIDS_H
#define PAR_KMEDOIDS_H
///
/// @file par_kmedoids.h
/// @brief CAPEK and XCAPEK scalable parallel clustering algorithms.
/// 

#include <mpi.h>
#include <ostream>
#include <vector>
#include <functional>

#include <boost/iterator/permutation_iterator.hpp>

#include "Timer.h"
#include "kmedoids.h"
#include "multi_gather.h"
#include "trial.h"
#include "id_pair.h"
#include "par_partition.h"
#include "stl_utils.h"
#include "bic.h"
#include "mpi_bindings.h"

namespace cluster {

  ///
  /// This class implements the CAPEK and XCAPEK scalable parallel clustering algorithms.
  /// 
  /// For more information on these algorithms, see this paper:
  /// @par
  /// Todd Gamblin, Bronis R. de Supinski, Martin Schulz, Rob Fowler, and Daniel A. Reed.
  /// <a href="http://www.cs.unc.edu/~tgamblin/pubs/scalable-cluster-ics10.pdf">
  /// <b>Clustering Performance Data Efficiently at Massive Scales</b></a>.
  /// <i>Proceedings of the International Conference on Supercomputing (ICS'10)</i>,
  /// Tsukuba, Japan, June 1-4, 2010.
  ///
  /// <b>Example usage:</b>
  /// @code
  /// // This is a theoretical point struct to be clustered.
  /// struct point {
  ///     double x, y;
  /// };
  /// 
  /// // Euclidean distance functor to use for clustering.
  /// struct euclidean_distance {
  ///     double operator()(const point& lhs, const point& rhs) {
  ///         double a = lhs.x - rhs.x;
  ///         double b = lhs.y - rhs.y;
  ///         return sqrt(a*a + b*b);
  ///     }
  /// };
  ///
  /// vector<point> points;
  /// // ... Each process puts some local points in the vector ...
  /// 
  /// par_kmedoids parkm(MPI_COMM_WORLD);
  /// vector<point> medoids;   // storage for reprsentatives
  /// 
  /// // Run xcapek(), finding a max of 20 clusters among the 2d points.
  /// parkm.xcapek(points, euclidean_distance(), 20, 2, &medoids);
  /// @endcode
  /// 
  /// After this runs, these members of <code>parkm</code> are interesting:
  /// - <code>parkm.\link par_partition::cluster_ids cluster_ids\endlink</code>: 
  ///    A vector of cluster ids for the local objects in <code>points</code> 
  /// - <code>parkm.\link par_partition::medoid_ids medoid_ids\endlink</code>:
  ///    A vector of object ids for all the cluster representatives
  /// 
  /// Together, these define the clustering that the algorithm found. See par_partition 
  /// for an explanation of how distributed partitions like this are represented.
  /// 
  /// The <code>medoids</code> vector contains actual copies of the representatives for each cluster.
  /// The copies correspond to the object ids in <code>parkm.medoid_ids</code>.  Supplying
  /// the <code>medoids</code> vector like this is optional, so if you don't need copies of
  /// the representatives, you can omit it from the call.
  /// 
  /// @endcode
  ///
  class par_kmedoids : public par_partition {
  public:
    ///
    /// Constructs a parallel kmedoids object and seeds its random number generator.
    /// This is a collective operation, and needs to be called by all processes.
    ///
    /// par_kmedoids assumes that there is one object to be clustered per process.
    ///
    par_kmedoids(MPI_Comm comm = MPI_COMM_WORLD);

    virtual ~par_kmedoids() { }

    /// Get the average dissimilarity of objects w/their medoids for the last run.
    double average_dissimilarity();

    /// BIC score for selected clustering
    double bic_score();

    ///
    /// Sets max_reps, Max number of times to run PAM with each sampled dataset.
    /// Default is 5, per Kaufman and Rousseeuw.
    ///
    void set_max_reps(size_t reps) { max_reps = reps; }

    ///
    /// Max number of times to run PAM with each sampled dataset.
    ///
    size_t get_max_reps() { return max_reps; }
    
    ///
    /// Sets init_size, baseline size for samples, added to 2*k.
    /// Defaults to 40, per Kaufman and Rousseeuw.
    ///
    void set_init_size(size_t size) { init_size = size; }
    
    ///
    /// Baseline size for samples, added to 2*k.
    ///
    size_t get_init_size() { return init_size; }

    /// Set tolerance for convergence.  This is relative error, not absolute error.  It will be
    /// multiplied by the mean distance before it is used to test convergence.
    /// Defaults to 1e-15; may need to be higher if there exist clusterings with very similar quality.
    void set_epsilon(double epsilon);


    ///
    /// Farms out trials of PAM to worker processes then collects medoids from all trials to all processors.
    /// Puts resulting medoids in all_medoids when done.
    ///
    template <class T, class D>
    void run_pam_trials(trial_generator& trials, const std::vector<T>& objects, D dmetric, 
                        std::vector<typename id_pair<T>::vector>& all_medoids, MPI_Comm comm)
    {
      int size, rank;
      CMPI_Comm_size(comm, &size);
      CMPI_Comm_rank(comm, &rank);
      
      for (size_t i=0; trials.has_next(); i++) {
        int my_k = -1;                        // trial id for local run of kmedoids
        int my_trial = -1;                    // trial id for local run of kmedoids
        std::vector<size_t>  my_ids;          // object ids for each of my_objects
        std::vector<T>       my_objects;      // vector to hold local sample of objects for clustering.
        multi_gather<T>      gather(comm);    // simultaneous, asynchronous local gathers for collecting samples.
        
        // start gathers for each trial to aggregate samples to single worker processes.
        for (int root=0; trials.has_next() && root < size; root++) {
          trial cur_trial = trials.next();    // generate a trial descriptor
          
          // Generate a set of indices for members of this k-medoids trial
          std::vector<size_t> sample_ids;
          boost::random_number_generator<random_t> rng(random);  // Boost adaptor for STL RNG's
          random_subset(trials.num_objects, cur_trial.sample_size, std::back_inserter(sample_ids), rng);

          // figure out where the sample objects live, ASSUME objects.size() objs per process.
          std::vector<int> sources;
          std::transform(sample_ids.begin(), sample_ids.end(), std::back_inserter(sources),
                         std::bind2nd(std::divides<size_t>(), objects.size()));
          sources.erase(std::unique(sources.begin(), sources.end()), sources.end());

          // make a permutation vector for the indices of the sampled *local* objects
          std::vector<size_t> sample_indices;
          transform(std::lower_bound(sample_ids.begin(), sample_ids.end(), objects.size() * rank),
                    std::lower_bound(sample_ids.begin(), sample_ids.end(), objects.size() * (rank + 1)),
                    std::back_inserter(sample_indices),
                    std::bind2nd(std::minus<int>(), objects.size() * rank));

          // gather trial members to the current worker (root)
          gather.start(boost::make_permutation_iterator(objects.begin(), sample_indices.begin()), 
                       boost::make_permutation_iterator(objects.begin(), sample_indices.end()),
                       sources.begin(), sources.end(), my_objects, root);
          
          // record which trial to use locally and save the medoids there.
          if (rank == root) {
            my_k     = cur_trial.k;
            my_trial = trials.count() - 1;
            my_ids.swap(sample_ids);
          }
        }
        timer.record("StartGather");
        
        // finish all sample gathers.
        gather.finish();
        timer.record("FinishGather");

        // if we're a worker process (we were assigned a k and a trial number) then run PAM on the sample.
        if (my_k >= 0) {
          kmedoids cluster;
          cluster.set_epsilon(epsilon);

          dissimilarity_matrix mat;
          build_dissimilarity_matrix(my_objects, dmetric, mat);
          cluster.pam(mat, my_k);

          // put this trial's medoids into its spot in the global medoids array.
          for (size_t m=0; m < cluster.medoid_ids.size(); m++) {
            all_medoids[my_trial].push_back(
              make_id_pair(my_objects[cluster.medoid_ids[m]], my_ids[cluster.medoid_ids[m]]));
          }
        }
        timer.record("LocalCluster");

        // once workers are done clustering, broadcast the medoids from each clustering 
        // so that all processes have the full all_medoids array.
        for (size_t trial_id = i * size; trial_id < trials.count(); trial_id++) {
          bcast_vector(comm, all_medoids[trial_id], trial_id % size);
        }

        timer.record("Broadcast");
      }
    }

    ///
    /// This is the Clustering Algorithm with Parallel Extensions to K-Medoids (CAPEK).
    /// 
    /// Assumes that objects to be clustered are fully distributed across parallel process, 
    /// with the same number of objects per process.  
    ///
    /// @tparam T     Type of objects to be clustered.
    ///               T must support the following operations:
    ///               - <code>int packed_size(MPI_Comm comm) const</code>
    ///               - <code>void pack(void *buf, int bufsize, int *position, MPI_Comm comm) const</code>
    ///               - <code>static T unpack(void *buf, int bufsize, int *position, MPI_Comm comm)</code>
    /// @tparam D     Dissimilarity metric type.  
    ///               D should be callable on (T, T) and should return a double representing 
    ///               the distance between the two T's.
    /// 
    /// @param[in]  objects   Local objects to cluster (ASSUME: currently must be same number per process!)
    /// @param[in]  dmetric   Distance metric to build dissimilarity matrices with
    /// @param[in]  k         Number of clusters to find.
    /// @param[out] medoids   Optional output vector where global medoids will 
    ///                       be stored along with their source ranks.
    ///
    /// CAPEK will run trials insatances of PAM on separate processors for each k in 1..max_k using the
    /// run_pam_trials() routine. Each trial aggregates sample_size objects distributed over all 
    /// processes in the system.
    ///
    /// @see xcapek() for a K-agnostic version of this algorithm.
    ///
    template <class T, class D>
    void capek(const std::vector<T>& objects, D dmetric, size_t k, std::vector<T> *medoids = NULL) 
    {
      int size, rank;
      CMPI_Comm_size(comm, &size);
      CMPI_Comm_rank(comm, &rank);

      seed_random_uniform(comm); // seed RN generator uniformly across ranks.

      // fix things if k is greater than the number of elements, since we can't 
      // ever find that many clusters.
      size_t num_objects = size * objects.size();
      k = std::min(num_objects, k);
      timer.record("Init");

      // do parallel work: farms out trials and broadcasts medoids from each trial to
      // all processes.  On completion, medoids from all trials are in all_medoids vector.
      std::vector<typename id_pair<T>::vector> all_medoids(max_reps);
      trial_generator trials(k, k, max_reps, init_size, num_objects);
      run_pam_trials(trials, objects, dmetric, all_medoids, comm);

      // Make two arrays to hold our closest medoids and their distance from our object
      std::vector<double> all_dissimilarities(trials.count(), 0.0);           // dissimilarity sums
      std::vector< std::vector<medoid_id> > all_cluster_ids(trials.count());  // local nearest medoid ids

      // Go through all the trials again, and for each of them, find the closest 
      // medoid to this process's objects and sum the dissimilarities
      for (size_t i=0; i < trials.count(); i++) {
        for (size_t o=0; o < objects.size(); o++) {
          object_id global_oid = rank + o;
          std::pair<double, size_t> closest = closest_medoid(objects[o], global_oid, all_medoids[i], dmetric);

          all_dissimilarities[i]  += closest.first;
          all_cluster_ids[i].push_back(closest.second);
        }
      }
      timer.record("FindMinima");

      // Sum up all the min dissimilarities.  We do a Reduce/Bcast instead of an Allreduce
      // to avoid FP error and guarantee that sums is the same across all processors.
      std::vector<double> sums(trials.count());         // destination vectors for reduction.

      CMPI_Reduce(&all_dissimilarities[0],  &sums[0], trials.count(), MPI_DOUBLE, MPI_SUM, 0, comm);
      CMPI_Bcast(&sums[0],  trials.count(), MPI_DOUBLE, 0, comm);
      timer.record("GlobalSums");

      // find minmum global dissimilarity among all trials.
      std::vector<double>::iterator min_sum = std::min_element(sums.begin(), sums.end());
      total_dissimilarity = *min_sum;
      size_t best = (min_sum - sums.begin());  // index of best trial.


      // Finally set up the partition to correspond to trial with best dissimilarity found
      medoid_ids.resize(all_medoids[best].size());
      for (size_t i = 0; i < medoid_ids.size(); i++) {
        medoid_ids[i] = all_medoids[best][i].id;
      }

      // Make an indirection vector from the unsorted to sorted medoids
      std::vector<size_t> mapping(medoid_ids.size());
      std::generate(mapping.begin(), mapping.end(), sequence());
      std::sort(mapping.begin(), mapping.end(), indexed_lt(medoid_ids));
      invert(mapping);

      // set up local cluster ids, medoids, and medoid_ids with the sorted mapping.
      for (size_t i=0; i < medoid_ids.size(); i++) {
        medoid_ids[i] = all_medoids[best][mapping[i]].id;
      }

      // swap in the cluster ids with the best BIC score.
      cluster_ids.swap(all_cluster_ids[best]);

      // if the caller wanted a copy of the medoids, copy them into the dstination array.
      if (medoids) {
        medoids->resize(medoid_ids.size());
        for (size_t i=0; i < medoid_ids.size(); i++) {
          (*medoids)[i] = all_medoids[best][mapping[i]].element;
        }
      }

      timer.record("BicScore");
    }    

    
    ///
    /// K-agnostic version of capek().
    /// This version attempts to guess the best K for the data using the 
    /// Bayesian Information Criterion (BIC) described in bic.h.  Evaluation of the 
    /// BIC is parallelized using global reduction operations.
    /// 
    /// Like capek(), this uses run_pam_trials() to farm out trials of the PAM clustering algorithm,
    /// but it requires more trials than capek().  In particular, it will run 
    /// (sum(1..max_k) * trials) total trials in parallel on MPI worker processes.
    ///
    /// @tparam T     Type of objects to be clustered.
    ///               T must support the following operations:
    ///                - <code>int packed_size(MPI_Comm comm) const</code>
    ///                - <code>void pack(void *buf, int bufsize, int *position, MPI_Comm comm) const</code>
    ///                - <code>static T unpack(void *buf, int bufsize, int *position, MPI_Comm comm)</code>
    /// @tparam D     Dissimilarity metric type.  
    ///               D should be callable on (T, T) and should return a double representing 
    ///               the distance between the two T's.
    /// 
    /// @param[in]  objects         Local objects to cluster (ASSUME: currently must be same number per process!)
    /// @param[in]  dmetric         Distance metric to build dissimilarity matrices with
    /// @param[in]  max_k           Max number of clusters to find.
    /// @param[in]  dimensionality  Dimensionality of objects, used by BIC.
    /// @param[out] medoids         Optional output vector where global medoids will be stored 
    ///                             along with their source ranks.
    ///
    /// @return
    /// The best BIC value found, that is, the BIC value of the final clustering.
    ///
    template <class T, class D>
    double xcapek(const std::vector<T>& objects, D dmetric, size_t max_k, size_t dimensionality,
                  std::vector<T> *medoids = NULL) 
    {
      int size, rank;
      CMPI_Comm_size(comm, &size);
      CMPI_Comm_rank(comm, &rank);

      seed_random_uniform(comm); // seed RN generator uniformly across ranks.

      // fix things if k is greater than the number of elements, since we can't 
      // ever find that many clusters.
      size_t num_objects = size * objects.size();
      max_k = std::min(num_objects, max_k);
      timer.record("Init");

      std::vector<typename id_pair<T>::vector> all_medoids(max_k * max_reps);
      trial_generator trials(max_k, max_reps, init_size, num_objects);
      run_pam_trials(trials, objects, dmetric, all_medoids, comm);

      // Make two arrays to hold our closest medoids and their distance from our object
      std::vector<double> all_dissimilarities(trials.count(), 0.0);           // dissimilarity sums
      std::vector< std::vector<medoid_id> > all_cluster_ids(trials.count());  // local nearest medoid ids

      std::vector<double> all_dissim2;      // dissimilarity sums squared
      std::vector<size_t> cluster_sizes;    // sizes of clusters in each trial

      // Go through all the trials again, and for each of them, find the closest 
      // medoid to this process's objects and sum the squared dissimilarities
      for (size_t i=0; i < trials.count(); i++) {
        const size_t num_medoids = all_medoids[i].size();
        for (size_t m=0; m < num_medoids; m++) {
          all_dissim2.push_back(0.0);
          cluster_sizes.push_back(0);
        }
        double *dissim2 = &all_dissim2[all_dissim2.size() - num_medoids];
        size_t *sizes   = &cluster_sizes[cluster_sizes.size() - num_medoids];
        
        for (size_t o=0; o < objects.size(); o++) {
          object_id global_oid = rank + o;
          std::pair<double, size_t> closest = closest_medoid(objects[o], global_oid, all_medoids[i], dmetric);

          all_dissimilarities[i]  += closest.first;
          dissim2[closest.second] += closest.first * closest.first;
          sizes[closest.second]   += 1;
          all_cluster_ids[i].push_back(closest.second);
        }
      }
      timer.record("FindMinima");
      

      // Sum up all the min dissimilarities.  We do a Reduce/Bcast instead of an Allreduce
      // to avoid FP error and guarantee that sums is the same across all processors.
      std::vector<double> sums(trials.count());         // destination vectors for reduction.
      std::vector<double> sums2(all_dissim2.size());
      std::vector<size_t> sizes(cluster_sizes.size());

      CMPI_Reduce(&all_dissimilarities[0],  &sums[0], trials.count(), MPI_DOUBLE, MPI_SUM, 0, comm);
      CMPI_Reduce(&all_dissim2[0], &sums2[0], sums2.size(), MPI_DOUBLE, MPI_SUM, 0, comm);
      CMPI_Bcast(&sums[0],  trials.count(), MPI_DOUBLE, 0, comm);
      CMPI_Bcast(&sums2[0], sums2.size(), MPI_DOUBLE, 0, comm);
      CMPI_Allreduce(&cluster_sizes[0], &sizes[0], sizes.size(), MPI_SIZE_T, MPI_SUM, comm);
      timer.record("GlobalSums");

      // find minmum global dissimilarity among all trials.
      std::vector<double>::iterator min_sum = std::min_element(sums.begin(), sums.end());
      total_dissimilarity = *min_sum;
      size_t best = (min_sum - sums.begin());  // index of best trial.

      // locally calculate the BIC for each trial
      size_t best_trial = 0;
      best_bic_score    = -DBL_MAX;

      size_t trial_offset = 0;  // offset into sizes array
      for (size_t i=0; i < trials.count(); i++) {
        size_t k = all_medoids[i].size();
        double cur_bic = bic(k, &sizes[trial_offset], &sums2[trial_offset], dimensionality);
        if (cur_bic > best_bic_score) {
          best_trial     = i;
          best_bic_score = cur_bic;
        }
        trial_offset += k;
      }
      
      best = best_trial;

      // Finally set up the partition to correspond to best trial found.
      medoid_ids.resize(all_medoids[best].size());
      for (size_t i = 0; i < medoid_ids.size(); i++) {
        medoid_ids[i] = all_medoids[best][i].id;
      }

      // Make an indirection vector from the unsorted to sorted medoids.
      std::vector<size_t> mapping(medoid_ids.size());
      std::generate(mapping.begin(), mapping.end(), sequence());
      std::sort(mapping.begin(), mapping.end(), indexed_lt(medoid_ids));
      invert(mapping);

      // set up local cluster ids, medoids, and medoid_ids with the sorted mapping.
      for (size_t i=0; i < medoid_ids.size(); i++) {
        medoid_ids[i] = all_medoids[best][mapping[i]].id;
      }

      // swap in the cluster ids with the best BIC score.
      cluster_ids.swap(all_cluster_ids[best]);

      // if the user wanted a copy of the medoids, copy them into the dstination array.
      if (medoids) {
        medoids->resize(medoid_ids.size());
        for (size_t i=0; i < medoid_ids.size(); i++) {
          (*medoids)[i] = all_medoids[best][mapping[i]].element;
        }
      }

      timer.record("BicScore");
      return best_bic_score;
    }    
    
    /// Get the Timer with info on the last run of either capek() or xcapek().
    const Timer& get_timer() { return timer; }

  protected:
    typedef boost::mt19937 random_t;   ///< Type for random number generator used here.
    random_t random;                   ///< Random number distribution to be used for samples
    
    double total_dissimilarity;   ///< Total dissimilarity bt/w objects and medoids for last clustering.
    double best_bic_score;        ///< BIC score for the clustering found.
    size_t init_size;             ///< baseline size for samples
    size_t max_reps;              ///< Max repetitions of trials for a particular k.
    double epsilon;               ///< Tolerance for convergence tests in kmedoids PAM runs.

    Timer timer;                  ///< Performance timer.

    /// 
    /// Seeds random number generators across all processes with the same number,
    /// taken from the time in microseconds since the epoch on the process 0.
    /// 
    void seed_random_uniform(MPI_Comm comm);

    ///
    /// This function broadcasts a vector of packable objects on one process 
    /// to all processes.  To be packable, a type must support these operations:
    ///                  T must also support the following operations:
    ///                    - int packed_size(MPI_Comm comm) const
    ///                    - void pack(void *buf, int bufsize, int *position, MPI_Comm comm) const
    ///                    - static T unpack(void *buf, int bufsize, int *position, MPI_Comm comm)
    ///
    template <class T>
    void bcast_vector(MPI_Comm comm, std::vector<T>& medoids, int root) {
      int rank;
      CMPI_Comm_rank(comm, &rank);

      // figure out size of packed buffer
      int packed_size=0;
      if (rank == root) {
        packed_size += cmpi_packed_size(1, MPI_SIZE_T, comm);   // num medoids for trial.
        for (size_t i=0; i < medoids.size(); i++) {         // size of medoids.
          packed_size += medoids[i].packed_size(comm);
        }
      }

      // broadcast size and allocate receive buffer.
      CMPI_Bcast(&packed_size, 1, MPI_INT, root, comm);
      char packed_buffer[packed_size];

      // pack buffer with medoid objects
      if (rank == root) {
        int pos = 0;

        size_t num_medoids = medoids.size();
        CMPI_Pack(&num_medoids, 1, MPI_SIZE_T, packed_buffer, packed_size, &pos, comm);

        for (size_t i=0; i < num_medoids; i++) {
          medoids[i].pack(packed_buffer, packed_size, &pos, comm);
        }
      }

      CMPI_Bcast(&packed_buffer, packed_size, MPI_PACKED, root, comm);

      // unpack broadcast medoids.
      if (rank != root) {
        int pos = 0;

        size_t num_medoids;
        CMPI_Unpack(packed_buffer, packed_size, &pos, &num_medoids, 1, MPI_SIZE_T, comm);
        medoids.resize(num_medoids);

        for (size_t i=0; i < medoids.size(); i++) {
          medoids[i] = T::unpack(packed_buffer, packed_size, &pos, comm);
        }
      }
    }

    ///
    /// Find the closest object in the medoids vector to the object passed in.
    /// Returns a pair of the closest medoid's id and its distance from the object.
    /// 
    /// @param[in] object   Object to compare to medoids.
    /// @param[in] oid      ID of the object (need this so medoids prefer themselves as their own medoids).
    /// @param[in] medoids  Vector of medoids to find the closest from.
    /// @param[in] dmetric  Distance metric to assess closeness with.
    ///
    template <typename T, typename D>
    std::pair<double, size_t> closest_medoid(
      const T& object, object_id oid, const std::vector< id_pair<T> >& medoids, D dmetric
    ) {
      double min_distance = DBL_MAX;
      size_t min_id = medoids.size();
      for (size_t m=0; m < medoids.size(); m++) {
        double d = dmetric(medoids[m].element, object);
        if (d < min_distance || medoids[m].id == oid) { // prefer actual medoid as closest
          min_distance = d;
          min_id = m;
        }
      }
      return std::make_pair(min_distance, min_id);
    }

  };

} // Namespace cluster

#endif // PAR_KMEDOIDS_H
