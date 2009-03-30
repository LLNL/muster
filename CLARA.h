#ifndef CLARA_H
#define CLARA_H

#include <algorithm>
#include <vector>

#include "KMedoids.h"

/**
 * Implementation of CLARA clustering algorithm, as per
 * R. Ng and J. Han, Efficient and Effective Clustering Methods for Spatial Data Mining.
 * 
 * by Todd Gamblin, 2004
 */
class CLARA : public KMedoids {
private:
    /** Paper reports that 5 runs of KMedoids yields good results. */
    static const unsigned ITERATIONS = 5;

    double dissimilarity;

    /**
     * Finds average dissimilarity of all clusters.
     */
    double averageDissimilarity();

    /**
     * Returns the sample size for the dataset CLARA is running on.
     * Paper reports that 40+2k items is a good number.
     */
    unsigned getSampleSize() {
        return 40+2*k;
    }
    
public:
    /** Constructs a new CLARA instance */
    CLARA(ClusterDataSet& data, unsigned k);

    /** 
     * Runs CLARA algorithm to find clusters.
     */
    virtual void findClusters();
    
    /**
     * Gets the current calculated value for total dissimilarity in 
     * the clustering.
     */
    double getDissimilarity() { 
        return dissimilarity; 
    }
};

#endif //CLARA_H
