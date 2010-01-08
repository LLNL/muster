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
    const size_t R = std::accumulate(cluster_sizes, cluster_sizes + k, 0);

    // Shorthand for model dimensionality
    const size_t M = dimensionality;

    const double logR   = log(R);
    const double log2pi = log(2 * M_PI);
    const size_t pj     = k + M*k;          // free parameter count
    
    // apply criterion formula from paper.
    double criterion = 0;
    for (size_t i=0; i < k; i++) {
      const size_t Rn = *(cluster_sizes + i);
      const double s2 = *(sum2_dissim + i);
      criterion += 
        - (Rn / 2.0 * log2pi) 
        - (Rn * M * log(s2)) / 2.0 
        - (Rn - k) / 2.0 
        + Rn * log(Rn) 
        - Rn * logR;
    }
    criterion -= pj/2.0 * logR;
    
    return criterion;
  }
  
  ///
  /// BIC metric for partitions that have already been constructed.
  ///
  template <typename D>
  double bic(const partition& p, D distance, size_t M) {
    size_t k = p.num_clusters();
    std::vector<size_t> sizes(k);

    for (size_t i=0; i < k; i++) {
      sizes[i] = p.size(i);
    }

    std::vector<double> sum2_dissim(k, 0.0);
    for (size_t i=0; i < p.size(); i++) {
      double dissim = distance(i, p.medoid_ids[p.cluster_ids[i]]);
      sum2_dissim[p.cluster_ids[i]] += dissim * dissim;
    }

    return bic(k, sizes.begin(), sum2_dissim.begin(), M);
  }



};

#endif // BAYESIAN_INFORMATION_CRITERION_H
