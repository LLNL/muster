#ifndef K_MEDOIDS_H
#define K_MEDOIDS_H


#include <vector>
#include <set>
#include <iostream>
using namespace std;

#include "ClusterDataSet.h"
#include "matrix.h"

/**
 * K-Mediod clustering method.  This selects q-grams to be medoids
 * at random and uses their distance function to cluster them.
 */
class KMedoids {
public:
    typedef vector< set<unsigned> > clusterList;
    typedef vector< set<unsigned> >::iterator clIterator;
    typedef set<unsigned> cluster;
    typedef set<unsigned>::iterator cIterator;

protected:

    /** More descriptive types for indices referred to in code. */
    typedef unsigned medoid_id;
    typedef unsigned object_id;
    
    /** Data set we'll do clustering on. */
    ClusterDataSet& data;

    /** Number of clusters we're searching for here. */
    const unsigned k;

    /** 
     * Medoids for each cluster.  
     * Maps medoid id (aka cluster id) to object id. 
     */
    vector<object_id> medoids;

    /**
     * Ids of clusters for each of the objects.  
     * Maps object id to cluster id. 
     */
    vector<medoid_id> clusterIds;

    /**
     * Puts randomly medoids from the sample set into 
     * medoids vector.
     */
    void assignRandomMedoids();

    /**
     * True iff object is one of our medoids.
     */
    bool isMedoid(unsigned object);

    /**
     * Finds the cost of swapping object Oh with medoid j.
     * PRE: object[oi] is the medoid with index clusterIds[oi] in medoids.
     */
    inline double cost(unsigned oi, unsigned oh, unsigned oj);

    /**
     * Total cost of swapping object h with medoid i.
     */
    inline double totalCost(medoid_id i, object_id h);


    /**
     * Assign each object to the cluster with the closest medoid.
     */
    void assignObjectsToClusters();

public:

    /**
     * Constructs a new instance of a KMedoids algorithm.
     * @param dataSet data set to cluster on
     * @param k number of medoids to search for
     * @param distance matrix containing distances between objects
     */
    KMedoids(
        ClusterDataSet& dataSet, 
        const unsigned k
    );

    virtual ~KMedoids();

    /** 
     * Run the K-Medoids clustering algorithm with the provided number of clusters.
     * @param k the number of clusters to find.
     * 
     * This method is not thread-safe.  Also, the vector passed in the constructor
     * should not change while this is running.
     */
    virtual void findClusters();

    /**
     * Returns a read-only reference to the dataset in 
     * this KMedoids instance.
     */
    const ClusterDataSet& getData();

    /** 
     * Returns a read-only view of the medoids we found. 
     * PRE: {@link #findClusters()} has been called.
     */
    const vector<unsigned>& getMedoids();
    
    /**
     * Returns clusters in the form of a newly allocated 
     * vector of sets of object indices.
     * The caller will need to delete the clusterLis.
     */
    virtual clusterList *getClustering();

    /** 
     * Writes out medoids' indices in this KMedoids' dataset's 
     * superset to the provided vector.
     */
    void getSuperSetMedoids(vector<unsigned>& dest);

    /**
     * Prints clustering out to the provided object stream.
     */
    void printClustering(ostream& out = cout);
    void printMedoids(ostream& out = cout);
    void printClusterIds(ostream& out = cout);

};

#endif //K_MEDOIDS_H

