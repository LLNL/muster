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

#ifndef CDBW_H
#define CDBW_H

///
/// @file cdbw.h
/// @brief Implementation of the CDbw criterion
///
/// CDbw criterion is a validity approach for density-based clustering algorithms.
/// It is based on a selection of multiple representatives for each cluster, and
/// evaluating different density aspect using these representatives.
///
/// For more on this technique and the approach, see this paper:
/// @par
/// Maria Halkidi and Michalis Vazirgiannis.  <a href="http://www.db-net.aueb.gr/index.php/corporate/content/download/439/2032/version/3/file/PATREC_revised.pdf">
/// <b>A Density-based Cluster Validity Approach using Multi-representatives</b></a>.
/// <i> Pattern Recognition Letters, 2008</i>, 
/// 

#include <boost/numeric/ublas/matrix.hpp>
#include <vector>
using boost::numeric::ublas::matrix;

#include <ANN/ANN.h>

#include "point.h"
#include "density_based.h"


namespace cluster
{

  class dbscan_cluster {
  private:
    std::vector<point>&  points;
    medoid_id            cluster_id;
    std::vector<size_t>  cluster_points;
    point                centroid;
    std::vector<size_t>  representatives;
    double               stdev;

  public:

    dbscan_cluster(size_t& cluster_id, std::vector<point>& points):
      cluster_id(cluster_id),points(points) { };
    
    void add_point(size_t point_id) {
      cluster_points.push_back(point_id);
    }

    void compute_data(void);

    void choose_representatives(size_t r);

    std::vector<size_t>& get_representatives(void) { return representatives; };

    size_t closest_representative(point& p);

    std::vector<point>   shrunk_representatives(double s);

    size_t size(void) { return cluster_points.size(); };
      
    double get_stdev(void)    { return stdev; };
    point& get_centroid(void) { return centroid; };

    void operator=(const dbscan_cluster& other) {
      points = other.points;
      cluster_id = other.cluster_id;
      cluster_points = other.cluster_points;
      centroid = other.centroid;
      representatives = other.representatives;
      stdev = other.stdev;
    }
  };

  ///
  /// Implementation of CDbw criterion. See paper 
  ///
  
  class CDbw {
  private:
    partition& p;
    std::vector<point>& points;
    
    std::vector<dbscan_cluster> clusters;

    size_t r; /// representatives
      
    boost::numeric::ublas::matrix<std::vector<std::pair<size_t, size_t> > > RCRs;

    ANNpointArray       ann_data_points;
    ANNkd_tree*         kd_tree;

    double              intra_cluster_density_change;
      
  public:
      
    CDbw(partition& p, std::vector<point>& points): p(p), points(points) {
      create_clusters();
    };

    double compute(size_t r);

  private:

    void create_clusters(void);

    /* RCRs computation */
    void compute_rcrs(void);

    std::vector< std::pair<size_t, size_t> > compute_rcrs_i_j(medoid_id i, medoid_id j);

    /* Separation */
    double separation(void);

    double inter_cluster_density(void);

    double density_between_clusters(medoid_id i, medoid_id j);

    double cardinality(point u, double radix, medoid_id i, medoid_id j);

    double distance_between_clusters(medoid_id i, medoid_id j);

    /* Compactness */
    double compactness_and_intra_density_changes(void);

    double intra_cluster_density(double s);

    double density(double s);

    double cardinality(point u, double radix, medoid_id i);

    /* Cohesion */
    double cohesion(double compact);

    /* Range query */
    std::vector<size_t> range_query(point u, double radix);
  };

};

#endif // CDBW_H
