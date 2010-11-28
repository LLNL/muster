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

#include "cdbw.h"

#include <cmath>
#include <deque>
using namespace std;

namespace cluster {

  /*****************************************************************************
   * dbscan_clusters
   ****************************************************************************/

  ///
  /// Compute the centroid and the standard deviation for each cluster data
  ///
  void dbscan_cluster::compute_data(void) {
    if (cluster_points.size() == 0)
      return;
    
    /* Compute centroid */
    centroid = points[cluster_points[0]];

    if (cluster_points.size() > 1) {
      for (size_t i = 1; i < cluster_points.size(); i++) {
        centroid += points[cluster_points[i]];
      }

      centroid /= cluster_points.size();
    }

    /* Stdev */
    double sum_squares = 0.0;
    for (size_t i = 0; i < cluster_points.size(); i++) {
      sum_squares += pow(centroid.distance(points[cluster_points[i]]), 2);
    }

    // cout << "cl_" << cluster_id << " sum_squares = " << sum_squares << endl;
    
    stdev = sqrt(sum_squares/(cluster_points.size()-1));
  }

  ///
  /// Choose the representatives of each cluster using the 'fartherst-first' method
  ///
  void dbscan_cluster::choose_representatives(size_t r) {
    double max_distance = std::numeric_limits<double>::min();
    size_t representative_id;
    size_t representative_position;

    deque<bool> point_used (cluster_points.size(), false);

    if (r >= cluster_points.size()) {
      representatives = cluster_points;
      return;
    }
    
    point& last_representative = centroid;
    
    representatives = vector<size_t>(r, 0);
  
    for (size_t i = 0; i < r; ++i) {
      max_distance = std::numeric_limits<double>::min();
      representative_id = 0;
      
      for (size_t j = 0; j < cluster_points.size(); ++j) {
        // cout << "checking point " << cluster_points[j] << "( " << point_used[j] << ")" << endl;
        if (!point_used[j]) {
          double current_distance = last_representative.distance(points[cluster_points[j]]);
          
          if (current_distance > max_distance) {
            max_distance            = current_distance;
            representative_id       = cluster_points[j];
            representative_position = j;
          }
        }
      }

      representatives[i]                  = representative_id;
      point_used[representative_position] = true;
      last_representative                 = points[representative_id];
    }
    
    return;
  }
  
  ///
  /// Returns the closest representative of the given point
  ///
  size_t dbscan_cluster::closest_representative(point& p) {
    double min_distance = p.distance(points[representatives[0]]);
    size_t result       = representatives[0];
    
    for (size_t i = 1; i < representatives.size(); i++) {
      double current_distance = p.distance(points[representatives[i]]);
      
      if (current_distance < min_distance) {
        min_distance = current_distance;
        result       = representatives[i];
      }
    }

    return result;
  }

  ///
  /// Returns a vector of the representative points shrunken through the centroid with an 's' factor
  ///
  vector<point> dbscan_cluster::shrunk_representatives(double s) {
    vector<point> result;
    
    for (size_t i = 0; i < representatives.size(); ++i) {
      point shrunken_point = points[representatives[i]];

      shrunken_point.x = shrunken_point.x + s*(centroid.x - shrunken_point.x);
      shrunken_point.y = shrunken_point.y + s*(centroid.y - shrunken_point.y);

      result.push_back(shrunken_point);
    }
    return result;
  }

  /*****************************************************************************
   * CDbw
   ****************************************************************************/

  ///
  /// Compute the CDbw criterion
  ///
  double CDbw::compute(size_t r) {
    double result;

    double sep, coh, compact, sc;

    // Criterion is not defined when number of clusters is 1
    if (p.num_clusters() == 1) {
      if (std::numeric_limits<double>::has_quiet_NaN) {
        return std::numeric_limits<double>::quiet_NaN();

      } else if (std::numeric_limits<double>::has_infinity) {
        return std::numeric_limits<double>::infinity();
          
      } else {
        return std::numeric_limits<double>::max();
      }
    }
    
    this->r = r;

    this->RCRs = matrix<vector<pair<size_t, size_t> > >(p.num_clusters(), p.num_clusters());

  
    for (size_t i = 0; i < p.num_clusters(); i++) {
      clusters[i].choose_representatives(r);
    }

    compute_rcrs();
  
    sep         = separation();
    compact     = compactness_and_intra_density_changes();
    coh         = cohesion(compact);
    sc          = sep*compact;

    result      = coh*sc;

    // DEBUG
    cout << "-- CDbw --" << endl;
    cout << "Separation = " << sep << endl;
    cout << "Compactness = " << compact << endl;
    cout << "Cohesion = " << coh << endl;
    cout << "Separation wrt Compactness = " << sc << endl;

    return result;
  }

  ///
  /// Create the 'dbscan_cluster' objects using, so as to manipulate the cluster data easily
  ///
  void CDbw::create_clusters(void) {
    for (size_t i = 0; i < p.num_clusters(); i++) {
      clusters.push_back(dbscan_cluster(i, points));
    }

    ann_data_points = annAllocPts(p.size(), 2);
  
    for (size_t i = 0; i < p.size(); i++) {
      if (get_cluster(i) != density::UNCLASSIFIED &&
          get_cluster(i) != density::NOISE) {
        // -2 -> first cluster value == 2
        clusters[(size_t) (get_cluster(i)-2)].add_point(i);
      }
      
      ANNpoint current_point = annAllocPt(2);
      
      current_point[0] = points[i].x;
      current_point[1] = points[i].y;
      
      ann_data_points[i] = current_point;
    }
    
    for (size_t i = 0; i < clusters.size(); i++) {
      clusters[i].compute_data();
      
      /*
        cout << "cl_" << i << " - Size = " << clusters[i].size();
        cout << " Centroid = " << clusters[i].get_centroid();
        cout << " stdev = " << clusters[i].get_stdev() << endl; */
    }
  
    kd_tree = new ANNkd_tree(ann_data_points, p.size(), 2);
  }

  ///
  /// Compute RCRs for each pair of clusters
  ///
  void CDbw::compute_rcrs(void) {
    for (size_t i = 0; i < clusters.size(); ++i) {
      for (size_t j = 0; j < clusters.size(); ++j) {
        if (i != j)
          RCRs(i,j) = compute_rcrs_i_j (i,j);
      }
    }
  }

  ///
  /// Compute the RCRs for a single pair of clusters
  ///
  vector<pair<size_t, size_t> > CDbw::compute_rcrs_i_j(medoid_id c_i, medoid_id c_j) {
    vector<pair<size_t, size_t> > result;

    vector<size_t>& reps_i = clusters[c_i].get_representatives();
    vector<size_t>& reps_j = clusters[c_j].get_representatives();

    vector<pair<size_t, size_t> > rcs_i_j;
    vector<pair<size_t, size_t> > rcs_j_i;

    for (size_t v = 0; v < reps_i.size(); ++v) {
      rcs_i_j.push_back(make_pair(reps_i[v],
                                  clusters[c_j].closest_representative(points[reps_i[v]])));
    }

    for (size_t v = 0; v < reps_j.size(); ++v) {
      rcs_j_i.push_back(make_pair(reps_j[v],
                                  clusters[c_i].closest_representative(points[reps_j[v]])));
    }

    // Now we have the rc's for both clusters, next, we have to compute the rcr's
    for (size_t vi = 0; vi < reps_i.size(); ++vi) {
      for (size_t vj = 0; vj < reps_j.size(); ++vj) {
        if (rcs_i_j[vi].first  == rcs_j_i[vj].second &&
            rcs_i_j[vi].second == rcs_j_i[vj].first) {
          result.push_back(rcs_i_j[vi]); // RCR pair found!
          continue;
        }
      }
    }

    /*
      cout << "-- Compute RCR btw Cluster " << c_i << " and Cluster " << c_j << " --" << endl;
      cout << "RCRs set size = " << result.size() << endl; */
  
    return result;
  }

  ///
  /// Compute the separation of the current partition C
  ///
  double CDbw::separation(void) {
    double result;
    double min_distance_btw_clusters_i_j;
    double sum_min_distance_btw_clusters = 0.0;
    double inter_cluster_density_val = inter_cluster_density();

    for (medoid_id i = 0; i < clusters.size(); ++i) {
      min_distance_btw_clusters_i_j = std::numeric_limits<double>::max();
    
      for (medoid_id j = 0; j < clusters.size(); ++j) {
        if (i != j) {
          double distance_btw_clusters = distance_between_clusters(i, j);
        
          if (distance_btw_clusters < min_distance_btw_clusters_i_j) {
            min_distance_btw_clusters_i_j = distance_btw_clusters;
          }
        }
      }

      sum_min_distance_btw_clusters += min_distance_btw_clusters_i_j;
    }

    /*
      cout << "-- Separation --" << endl;
      cout << "Sum. Min. Distance btw Clusters = " << sum_min_distance_btw_clusters << endl;
      cout << "Inter cluster density = " << inter_cluster_density_val << endl;
    */
  
    result = (sum_min_distance_btw_clusters/clusters.size())/(1 + inter_cluster_density_val);
  
    return result;
  }

  ///
  /// Compute the inter-cluster density of the current partition C
  ///
  double CDbw::inter_cluster_density(void) {
    double max_density_btw_clusters_i_j;
    double sum_max_density_btw_clusters = 0.0;

    for (medoid_id i = 0; i < clusters.size(); ++i) {
      max_density_btw_clusters_i_j = std::numeric_limits<double>::min();
      
      for (medoid_id j = 0; j < clusters.size(); ++j) {
        if (i != j) {
          double density_btw_clusters = density_between_clusters(i,j);

          if (density_btw_clusters > max_density_btw_clusters_i_j) {
            max_density_btw_clusters_i_j = density_btw_clusters;
          }
        }
      }

      sum_max_density_btw_clusters += max_density_btw_clusters_i_j;
    }

    /*
      cout << "-- Inter cluster density --" << endl;
      cout << "Max. Density btw Clusters = " << sum_max_density_btw_clusters << endl;
    */
  
    return (sum_max_density_btw_clusters/clusters.size());
  }

  ///
  /// Compute the density between a pair of clusters c_i and c_j
  ///
  double CDbw::density_between_clusters(medoid_id c_i, medoid_id c_j) {
    vector<pair<size_t, size_t> >& RCR_i_j = RCRs ((size_t) c_i, (size_t) c_j);

    double sum_densities = 0.0;
    double avg_stdev     = sqrt((pow(clusters[c_i].get_stdev(),2.0) + pow(clusters[c_j].get_stdev(),2.0))/2);

    for (size_t p = 0; p < RCR_i_j.size(); ++p) {
      point  u;
      double distance_vi_vj = points[RCR_i_j[p].first].distance(points[RCR_i_j[p].second]);
      double cardinality_u;
    
      u  = points[RCR_i_j[p].first];
      u += points[RCR_i_j[p].second];
      u /= 2.0;


    
      cardinality_u = cardinality(u, avg_stdev, c_i, c_j);

      /*
        cout << "cl_" << c_i << " & cl_" << c_j << " ";
        cout << "vij = " << points[RCR_i_j[p].first] << ", vji = " << points[RCR_i_j[p].second] << ", u = " << u;
        cout << " avg_stdev = " << avg_stdev;
        cout << " cardinality = " << cardinality_u << endl;
      */

      sum_densities += ((distance_vi_vj/(2*avg_stdev))*cardinality_u);
    }

    return (sum_densities/RCR_i_j.size());
  }

  ///
  /// Compute the cardinality of set defined by u as center point and radix, where the points
  /// belong to clusters c_i and c_j
  ///
  double CDbw::cardinality(point u, double radix, medoid_id c_i, medoid_id c_j) {
    size_t number_of_points = 0;
    // Perform the radix search (across all points), using libANN
    vector<size_t> neighbourhood_u = range_query(u, radix);
  
    // Count just those points that belong to clusters c_i and c_j
    for (size_t x = 0; x < neighbourhood_u.size(); ++x) {
      if ((get_cluster(neighbourhood_u[x])-2) == c_i ||
          (get_cluster(neighbourhood_u[x])-2) == c_j) {
        ++number_of_points;
      }
    }

    return (number_of_points/(clusters[c_i].size()+clusters[c_j].size()));
  }

  ///
  /// Computes the distance between clusters c_i and c_j
  ///
  double CDbw::distance_between_clusters(medoid_id c_i, medoid_id c_j) {
    double sum_distances = 0.0;
    vector<pair<size_t, size_t> >& RCR_i_j = RCRs (c_i, c_j);

    for (size_t p = 0; p < RCR_i_j.size(); ++p) {
      sum_distances += points[RCR_i_j[p].first].distance(points[RCR_i_j[p].second]);
    }

    /*
      cout << "-- Distance btw clusters " << c_i << " and " << c_j << " --" << endl;
      cout << "Sum distances = " << sum_distances << endl;
      cout << "RCR_i_j size = " << RCR_i_j.size() << endl;
    */
  
    return (sum_distances/RCR_i_j.size());
  }

  ///
  /// Computes the compactness of the partition and also the intra-cluster density changes, so
  /// as to save intra-cluster density computations
  ///
  double CDbw::compactness_and_intra_density_changes(void) {
    double         sum_intra_densities = 0.0;
    vector<double> intra_densities_values;
  
    for (double s = 0.1; s <= 0.8; s += 0.1) {
      double current_intra_density = intra_cluster_density(s);
      sum_intra_densities += current_intra_density;
    
      intra_densities_values.push_back(current_intra_density);
    }

    intra_cluster_density_change = 0.0;
  
    for (size_t i = 1; i < intra_densities_values.size(); ++i) {
      intra_cluster_density_change += fabs(intra_densities_values[i] - intra_densities_values[i-1]);
    }

    intra_cluster_density_change = intra_cluster_density_change/7;

    // cout << "intra_cluster_density_change = " << intra_cluster_density_change << endl;
  
    return (sum_intra_densities/8); // 8 -> number of different 's'
  }

  double CDbw::intra_cluster_density(double s) {
    double avg_stdev = 0.0;
    double result;

    for (size_t i = 0; i < clusters.size(); ++i) {
      avg_stdev += pow (clusters[i].get_stdev(), 2.0);
    }

    avg_stdev = sqrt(avg_stdev/clusters.size());

    result = density(s)/(clusters.size()*avg_stdev);
  
    return result;
  }

  double CDbw::density(double s) {
    double result;
    double sum_cardinalities = 0.0;

    for (size_t i = 0; i < clusters.size(); ++i) {
      vector<point> shrunk_rep = clusters[i].shrunk_representatives(s);

      for (size_t j = 0; j < shrunk_rep.size(); ++j) {
        sum_cardinalities += cardinality(shrunk_rep[j], clusters[i].get_stdev(), i);
      }
    }

    // cout << "sum_cardinalities = " << sum_cardinalities << endl;
  
    result = sum_cardinalities/r;

    /* DEBUG */
    // cout << "Density s_" << s << " = " << result << endl;

    return result;
  }

  double CDbw::cardinality(point u, double radix, medoid_id c_i) {
    size_t number_of_points = 0;
    double result;
    // Perform the radix search (across all points), using libANN
    vector<size_t> neighbourhood_u = range_query(u, radix);

    // cout << "neighbourhood size = " << neighbourhood_u.size() << endl;
  
    // Count just those points that belong to cluster c_i
    for (size_t x = 0; x < neighbourhood_u.size(); ++x) {
      if ((get_cluster(neighbourhood_u[x])-2) == c_i) {
        ++number_of_points;
      }
    }

    // cout << "number of points = " << number_of_points << endl;

    result = number_of_points*1.0/clusters[c_i].size();

    // cout << "cardinality = " << result << endl;

    return result;
  }

  /* Cohesion */
  double CDbw::cohesion(double compact) {
    double result = compact/(1+intra_cluster_density_change);

    return result;
  }

  ///
  /// Performs a range query in the KD Tree, centered on point u, with the radix indicated
  ///
  vector<size_t> CDbw::range_query(point u, double radix) {
    ANNpoint       ann_query_point = annAllocPt(2);
    ANNidxArray    ann_results;
    size_t         results_size;
    vector<size_t> result;

    ann_query_point[0] = u.x;
    ann_query_point[1] = u.y;

    results_size = kd_tree->annkFRSearch(ann_query_point, pow(radix, 2.0), 0);
  
    ann_results = new ANNidx[results_size];

    kd_tree->annkFRSearch(ann_query_point,
                          pow(radix, 2.0),
                          results_size,
                          ann_results);

    for (size_t i = 0; i < results_size; ++i) {
      result.push_back(ann_results[i]);
    }

    return result;
  }

  medoid_id CDbw::get_cluster(object_id i) {
    if (i < p.size()) {
      return p.cluster_ids[i];
    } else {
      return density::UNCLASSIFIED;
    }
  }

} // namespace cluster
