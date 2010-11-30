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
/// @file cdbw.h
/// @author Juan Gonzalez juan.gonzalez@bsc.es
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
#ifndef CDBW_H
#define CDBW_H

#include <boost/numeric/ublas/matrix.hpp>
#include <vector>
using boost::numeric::ublas::matrix;

#include <ANN/ANN.h>

#include "point.h"
#include "density.h"


namespace cluster {

  class dbscan_cluster {
  private:
    std::vector<point>&  points_;
    medoid_id            cluster_id_;
    std::vector<size_t>  cluster_points_;
    point                centroid_;
    std::vector<size_t>  representatives_;
    double               stdev_;

  public:

    dbscan_cluster(size_t& cluster_id, std::vector<point>& points);

    dbscan_cluster(const dbscan_cluster& other);

    ~dbscan_cluster();
    
    void add_point(size_t point_id);

    void compute_data();

    void choose_representatives(size_t r);

    size_t closest_representative(point& p);

    std::vector<point> shrunk_representatives(double s);

    size_t size() { return cluster_points_.size(); }

    std::vector<size_t>& representatives()   { return representatives_; }
    double stdev()                           { return stdev_; }
    point& centroid()                        { return centroid_; }

    void operator=(const dbscan_cluster& other);
  };

  ///
  /// Implementation of CDbw criterion. See paper 
  ///
  
  class CDbw {
  private:
    partition& p_;
    std::vector<point>& points_;
    std::vector<dbscan_cluster> clusters_;
    size_t r_;                              ///< Representatives
      
    boost::numeric::ublas::matrix<std::vector<std::pair<size_t, size_t> > > RCRs_;

    ANNpointArray       ann_data_points_;
    ANNkd_tree*         kd_tree_;
    double              intra_cluster_density_change_;
      
  public:
      
    CDbw(partition& p, std::vector<point>& points);

    ~CDbw();

    double compute(size_t r);

  private:

    void create_clusters();

    /* RCRs computation */
    void compute_rcrs();

    std::vector< std::pair<size_t, size_t> > compute_rcrs_i_j(medoid_id i, medoid_id j);

    /* Separation */
    double separation();

    double inter_cluster_density();

    double density_between_clusters(medoid_id i, medoid_id j);

    double cardinality(point u, double radix, medoid_id i, medoid_id j);

    double distance_between_clusters(medoid_id i, medoid_id j);

    /* Compactness */
    double compactness_and_intra_density_changes();

    double intra_cluster_density(double s);

    double density(double s);

    double cardinality(point u, double radix, medoid_id i);

    /* Cohesion */
    double cohesion(double compact);

    /* Range query */
    std::vector<size_t> range_query(point u, double radix);
  };

} // namespace cluster

#endif // CDBW_H
