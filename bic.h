#ifndef BAYESIAN_INFORMATION_CRITERION_H
#define BAYESIAN_INFORMATION_CRITERION_H

#include <stdint.h>
#include <numeric>
#include <iostream>
#include <cmath>
#include <vector>
#include "dissimilarity.h"
#include "partition.h"

namespace cluster {
  
  ///
  /// Evaluates the bayesian information criterion for a clustering.
  /// 
  /// Parameters:
  ///    k                   Number of clusters in the clustering.  Same as k from k-means or k-medoids.
  ///    cluster_sizes       Start of range of k sizes.  
  ///                          *cluster_sizes..*(cluster_sizes + k) == respective sizes of clusters 1 to k
  ///    sum2_dissimilarity  Sum of squared dissimilarities of each object w.r.t. its closest medoid.
  ///    dimensionality      Dimensionality of clustered data.  E.g.: 2 for 2-dimensional points.
  ///
  /// This implementation is based on "X-Means: Extending K-means with Efficient Estimation of the 
  /// Number of Clusters" by Dan Pelleg and Andrew Moore.
  ///
  template <typename SizeIterator, typename DissimIterator>
  double bic(size_t k, SizeIterator cluster_sizes, DissimIterator sum2_dissim, size_t dimensionality) {
    // figure out total size of data set and put it in R.
    const double R      = std::accumulate(cluster_sizes, cluster_sizes + k, 0);
    const double M      = dimensionality;    // Shorthand for model dimensionality
    const double logR   = log(R);
    const double log2pi = log(2 * M_PI);
    const double pj     = (k-1) + M*k + 1;   // free parameter count
    const double s2     = std::accumulate(sum2_dissim, sum2_dissim + k, 0.0) / (R - k);

    // apply criterion formula from xmeans paper.
    double criterion = 0;
    for (size_t i=0; i < k; i++) {
      const double Rn = *(cluster_sizes + i);
      criterion += 
        - (Rn * log2pi) / 2.0
        - (Rn * M * log(s2)) / 2.0 
        - (Rn - 1) / 2.0 
        + Rn * log(Rn) 
        - Rn * logR;
    }
    criterion -= (pj/2.0 * logR);
    
    return criterion;
  }


  ///
  /// BIC metric for partitions that have already been constructed.
  ///
  template <typename D>
  double old_bic(const partition& p, D distance, size_t M) {
    size_t k = p.num_clusters();

    size_t sizes[k];     // cluster sizes
    for (size_t i=0; i < k; i++) {
      sizes[i] = p.size(i);
    }

    std::vector<double> sum2_dissim(k);  // sum of squared dissimilarity
    for (size_t i=0; i < k; i++) {
      sum2_dissim[i] = total_squared_dissimilarity(p, distance, i);
    }

    return bic(k, sizes, sum2_dissim.begin(), M);
  }



  ///
  /// This is the direct calculation of the BIC from individual point probabilities.
  /// Use this to check that the partial implementation is correct.
  ///
  template <typename D>
  double bic(const partition& p, D distance, size_t M) {
    double R = p.size();
    size_t k = p.num_clusters();

    // calculate variance.
    double s2 = total_squared_dissimilarity(p, distance) / (R - k);
    double s  = sqrt(s2);
    double sM = pow(s, (double)M);
    
    size_t sizes[k];
    for (size_t i=0; i < k; i++) {
      sizes[i] = p.size(i);
    }
    
    double root2pi = sqrt(2 * M_PI);
    double lD = 0;
    for (size_t i=0; i < p.size(); i++) {
      double d  = distance(i, p.medoid_ids[p.cluster_ids[i]]);
      double Ri = sizes[p.cluster_ids[i]];
      lD += 
        + log(1.0 / (root2pi * sM))
        - (1 / (2 * s2)) * d * d
        + log(Ri / R);
    }

    const size_t pj = (k-1) + M*k + 1;   // free parameter count

/* Added BGLFE because no prototype exists for log(size_t). Only log(float) and
   log(long double). See math.h. Works fine on the backend
*/
#ifdef BGLFE
    return lD - pj/2 * log( (double) R);
#else
    return lD - pj/2 * log(R);
#endif

  }



};

#endif // BAYESIAN_INFORMATION_CRITERION_H
