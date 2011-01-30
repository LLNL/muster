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
/// @file density.h
/// @author Juan Gonzalez juan.gonzalez@bsc.es
/// @author Todd Gamblin tgamblin@llnl.gov
/// @brief Implementations of regular and sampled density clustering.
///
#ifndef MUSTER_DENSITY_H_
#define MUSTER_DENSITY_H_

#include <vector>
#include <set>
#include <list>
#include <iostream>
#include <stdexcept>
#include <cfloat>

#include <boost/random.hpp>

#include "CGAL/Origin.h"
#include "CGAL/Convex_hull_d.h"

#include "dissimilarity.h"
#include "partition.h"
#include "cdbw.h"

using namespace CGAL;

namespace cluster {

  /// 
  /// Implementation of the classic clustering algorithm DBSCAN.
  /// 
  class density : public partition {
  public:
    // Ids of special clusters in partitions created by dbscan.
    static const medoid_id NOISE         = 0;     ///< special id for noise cluster
    static const medoid_id FIRST_CLUSTER = 1;     ///< id of first real cluster

    ///
    /// Constructor.  
    /// 
    density(size_t num_objects = 0);

    /// Destructor does nothing for now.
    virtual ~density();

    /// 
    /// DBSCAN clustering, described by Ester et. al. in the paper "A Density-Based 
    /// Algorithm for Discovering Clusters in Large Spatial Databases with Noise"
    ///
    /// @tparam T    Type of objects to be clustered.
    /// @tparam D    Dissimilarity metric type.  D should be callable 
    ///              on (T, T) and should return a double.
    ///
    /// @param epsilon          maximum distance to perform the distance searches
    /// @param min_points       minimun number of points to consider a region as a cluster
    ///
    /// 
    template <class T, class D>
    void dbscan(const std::vector<T>& objects, D dmetric, double epsilon, size_t min_points) {
      epsilon_    = epsilon;
      min_points_ = min_points;

      // put an arbitrary representative in for the noise cluster
      medoid_ids.resize(1, 0);

      // init cluster ids to UNCLASSIFIED_ before we do the density clustering
      cluster_ids.clear();
      cluster_ids.resize(objects.size(), UNCLASSIFIED);

      // go through the data set and expand each point that's not yet classified.
      current_cluster_id_ = FIRST_CLUSTER;
      for (size_t i = 0; i < objects.size(); i++) {
        if (cluster_ids[i] == UNCLASSIFIED) {
          if (expand_cluster(objects, dmetric, i)) {
            medoid_ids.push_back(i);
            current_cluster_id_++;
          }
        }
      }

      // put a real noise object in as the representative for the noise cluster.
      size_t i=0;
      while (i < objects.size() && cluster_ids[i] != NOISE) i++;
      if (i < objects.size()) {
        medoid_ids[NOISE] = i;
      } else {
        // TODO: what to do if there is no noise?
      }
    }
    
    
    template <class D, class K>
    void sdbscan(const std::vector<typename K::Point_d>& objects, D dmetric, double epsilon, size_t min_points) {
      // Just run plain density once if sample size is larger than dataset.
      if (objects.size() <= sample_size_) {
        dissimilarity_matrix mat;
        dbscan(objects, dmetric, epsilon, min_points);
        return;
      }

      // medoids and clusters for best partition so far.
      double best_cdbw = 0;
      partition best_partition;

      //run KMedoids on a sampled subset max_reps times
      for (size_t i = 0; i < reps_; i++) {
        // Take a random sample of objects, store sample in a vector
        std::vector<size_t> sample_to_full;
        std::vector<typename K::Point_d> sample;
        algorithm_r(objects.size(), sample_size, back_inserter(sample_to_full), rng);
        for (size_t i=0; i < sample_size_; i++) {
          sample.push_back(objects[sample_to_full[i]]);
        }

        // run dbscan on the subset
        density subcall;
        subcall.dbscan(sample, dmetric, epsilon, min_points);
        subcall.remove_cluster(NOISE);
        
        // make convex hulls out of subset clusters
        std::vector< Convex_hull_d<K> > hulls[subcall.num_clusters()];
        for (size_t i=0; i < subcall.size(); i++) {
          if (subcall.cluster_ids[i] != UNCLASSIFIED) {
            hulls[subcall.cluster_ids[i]].insert(sample[i]);
          }
        }
        
        // classify full data set based on hulls
        medoid_ids.clear();
        cluster_ids.resize(objects.size());
        for (size_t i=0; i < objects.size(); i++) {
          cluster_ids[i] = index_of_containing_hull(hulls, objects[i]);
        }

        // now measure the cdbw of the set.
        CDbw<typename K::Point_d> cdbw_computer(*this, objects);
        double cdbw = cdbw_computer.compute(10);
        if (cdbw > best_cdbw) {
          swap(best_partition);
          best_cdbw = cdbw;
        }
      }
    }

    void set_sample_size(size_t sample_size);  ///< Sets the sample size for sdbscan.
    size_t sample_size();                      ///< @return sample size used by sdbscan.

    void set_reps(size_t sample_size);         ///< Sets the repetition count for sdbscan.
    size_t reps();                             ///< @return repetitions executed by sdbscan.
    
  protected:
    double epsilon_;              ///< maximum distance to perform the distance searches
    size_t min_points_;           ///< minimun number of points to consider a region as a cluster
    size_t current_cluster_id_;   ///< Next cluster id to assign
    size_t sample_size_;          ///< Sample size for sdbscan
    size_t reps_;                 ///< repetitions for sdbscan

    typedef boost::mt19937 random_type;                /// Type for RNG used in this algorithm
    random_type random;                                /// Randomness source for this algorithm
    
    /// Adaptor for STL algorithms.
    typedef boost::random_number_generator<random_type, unsigned long> rng_type;
    rng_type rng;


    template <class T, class D>
    bool expand_cluster(const std::vector<T>& objects, D dmetric, size_t current_object) {
      std::list<size_t> seed_list = epsilon_range_query(objects, dmetric, current_object);
      std::list<size_t>::iterator seed_list_iterator;

      if (seed_list.size() < min_points_) {
        cluster_ids[current_object] = NOISE;
        return false;
      }

      /* Assign current cluster id to current object neighborhood */
      seed_list_iterator = seed_list.begin();
      while (seed_list_iterator != seed_list.end()) {
        size_t current_seed = (*seed_list_iterator);
        cluster_ids[current_seed] = current_cluster_id_;

        if (current_seed == current_object) {
          seed_list_iterator = seed_list.erase(seed_list_iterator);
        } else {
          seed_list_iterator++;
        }
      }

      /* Expand the search to every seed */
      for (seed_list_iterator  = seed_list.begin();
           seed_list_iterator != seed_list.end();
           ++seed_list_iterator) {
          std::list<size_t>           neighbour_seed_list;
          std::list<size_t>::iterator neighbour_seed_list_iterator;
        
          size_t current_neighbour = (*seed_list_iterator);

          neighbour_seed_list = epsilon_range_query(objects, dmetric, current_neighbour);

          if (neighbour_seed_list.size() >= min_points_) {
            for (neighbour_seed_list_iterator  = neighbour_seed_list.begin();
                 neighbour_seed_list_iterator != neighbour_seed_list.end();
                 neighbour_seed_list_iterator++) {
              size_t current_neighbour_neighbour = (*neighbour_seed_list_iterator);

              if (cluster_ids[current_neighbour_neighbour] == UNCLASSIFIED ||
                  cluster_ids[current_neighbour_neighbour] == NOISE) {
                if (cluster_ids[current_neighbour_neighbour] == UNCLASSIFIED) {
                  seed_list.push_back(current_neighbour_neighbour);
                }
                cluster_ids[current_neighbour_neighbour] = current_cluster_id_;
              }
            }
          }
          neighbour_seed_list.clear();
        }

      return true;
    }

    template <class T, class D>
    std::list<size_t> epsilon_range_query(const std::vector<T>& objects, D dmetric, size_t current_object) {
      std::list<size_t> result;

      for (size_t i = 0; i < objects.size(); i++) {
        if ( i == current_object) {
          result.push_back(i);
          continue;
        }

        /* Check the distance between current_object and i-th object */
        if (dmetric(objects[current_object], objects[i]) < epsilon_) {
          result.push_back(i);
        }
      }

      return result;
    }

    /// Finds which hull in a vector of hulls contains a point.
    template<class K>
    medoid_id index_of_containing_hull(const std::vector< Convex_hull_d<K> >& hulls, typename K::Point_d p) {
      for (medoid_id i=0; i < (medoid_id)hulls.size(); i++) {
        if (hulls[i].boundary_side(p) != ON_UNBOUNDED_SIDE) {
          return i;
        }
      }
      return UNCLASSIFIED;
    }

  }; // class density

} // namespace cluster

#endif //MUSTER_DENSITY_H_
