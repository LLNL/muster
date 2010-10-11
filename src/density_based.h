//////////////////////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2010, Lawrence Livermore National Security, LLC.  
// Produced at the Lawrence Livermore National Laboratory  
// Written by Juan Gonzalez, juan.gonzalez@bsc.es

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


#ifndef K_MEDOIDS_H
#define K_MEDOIDS_H
///
/// @file kmedoids.h
/// @brief Implementations of the classic clustering algorithms PAM and CLARA, from 
/// <i>Finding Groups in Data</i>, by Kaufman and Rousseeuw.
///

#include <vector>
#include <set>
#include <list>
#include <iostream>
#include <stdexcept>
#include <cfloat>

#include <boost/random.hpp>

#include "random.h"
#include "dissimilarity.h"
#include "partition.h"
#include "bic.h"

using std::cout;
using std::endl;

namespace cluster {

  /// 
  /// Implementation of the classic clustering algorithm DBSCAN.
  /// 
  class density_based : public partition {
  public:
    ///
    /// Constructor.  
    /// 
    /// 
    density_based(size_t num_objects = 0);
    ~density_based();

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
    void dbscan(const std::vector<T>& objects, D dmetric, double _epsilon, size_t _min_points )
    {
      epsilon    = _epsilon;
      min_points = _min_points;

      for (size_t i = 0; i < objects.size(); i++)
      {
        cluster_ids.push_back(cluster::unclassified);
      }

      for (size_t i = 0; i < objects.size(); i++)
      {
        if (cluster_ids[i] == cluster::unclassified)
        {
          if (expand_cluster(objects, dmetric, i))
          {
            medoid_ids.push_back(i);
            current_cluster_id++;
            total_clusters++;
          }
        }
      }
    }

    protected:
    typedef boost::mt19937 random_type;                /// Type for RNG used in this algorithm
    random_type random;                                /// Randomness source for this algorithm
    
    /// Adaptor for STL algorithms.
    typedef boost::random_number_generator<random_type, unsigned long> rng_type;
    rng_type rng;

    double epsilon;                          /// maximum distance to perform the distance searches
    size_t min_points;                       /// minimun number of points to consider a region as a cluster

    size_t current_cluster_id;
    size_t total_clusters;

    template <class T, class D>
    bool expand_cluster(const std::vector<T>& objects, D dmetric, size_t current_object)
    {
      std::list<size_t> seed_list = epsilon_range_query(objects, dmetric, current_object);
      std::list<size_t>::iterator seed_list_iterator;

      if (seed_list.size() < min_points)
      {
        cluster_ids[current_object] = cluster::noise;
        return false;
      }

      /* Assign current cluster id to current object neighborhood */
      seed_list_iterator = seed_list.begin();
      while (seed_list_iterator != seed_list.end())
      {
        size_t current_seed = (*seed_list_iterator);
        cluster_ids[current_seed] = current_cluster_id;

        if (current_seed == current_object)
        {
          seed_list_iterator = seed_list.erase(seed_list_iterator);
        }
        else
        {
          seed_list_iterator++;
        }
      }

      /* Expand the search to every seed */
      for (seed_list_iterator  = seed_list.begin();
           seed_list_iterator != seed_list.end();
           ++seed_list_iterator)
      {
        std::list<size_t>           neighbour_seed_list;
        std::list<size_t>::iterator neighbour_seed_list_iterator;
        
        size_t current_neighbour = (*seed_list_iterator);

        neighbour_seed_list = epsilon_range_query(objects, dmetric, current_neighbour);

        if (neighbour_seed_list.size() >= min_points)
        {
          for (neighbour_seed_list_iterator  = neighbour_seed_list.begin();
               neighbour_seed_list_iterator != neighbour_seed_list.end();
               neighbour_seed_list_iterator++)
          {
            size_t current_neighbour_neighbour = (*neighbour_seed_list_iterator);

            if (cluster_ids[current_neighbour_neighbour] == cluster::unclassified ||
                cluster_ids[current_neighbour_neighbour] == cluster::noise)
            {
              if (cluster_ids[current_neighbour_neighbour] == cluster::unclassified)
              {
                seed_list.push_back(current_neighbour_neighbour);
              }
              cluster_ids[current_neighbour_neighbour] = current_cluster_id;
            }
          }
        }
        neighbour_seed_list.clear();
      }

      return true;
    }

    template <class T, class D>
    std::list<size_t> epsilon_range_query(const std::vector<T>& objects, D dmetric, size_t current_object)
    {
      std::list<size_t> result;

      double distance;
      
      for (size_t i = 0; i < objects.size(); i++)
      {
        if ( i == current_object)
        {
          result.push_back(i);
          continue;
        }

        /* Check the distance between current_object and i-th object */
        if (dmetric(objects[current_object], objects[i]) < epsilon)
        {
          result.push_back(i);
        }

      }

      return result;
    }
  };

} // namespace cluster

#endif //K_MEDOIDS_H
