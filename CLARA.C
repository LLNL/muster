#include "CLARA.h"

#include <cstdlib>

#include "Dissimilarity.h"

CLARA::CLARA(ClusterDataSet& _data, unsigned _k) :
    KMedoids(_data, _k)
{
    //nothing to say here.
}


double CLARA::averageDissimilarity() {
    double totalDissimilarity = 0.0;
    for (unsigned i=0; i < data.size(); i++) {
        totalDissimilarity = data.distance(i, medoids[clusterIds[i]]);
    }
    totalDissimilarity /= (double)data.size();

    return totalDissimilarity;
}

void CLARA::findClusters() {
    //run plain KMedoids once if sampling won't gain us anything
    if (data.size() <= getSampleSize()) {
        KMedoids::findClusters();
        return;
    }

    //TODO: could try to reduce copying here.
    vector<unsigned> bestMedoids(k);
    vector<unsigned> bestClustering(clusterIds.size());
    double bestDissimilarity = DISSIMILARITY_MAX;

    //run KMedoids on a sampled subset ITERATIONS times
    for (unsigned i = 0; i < ITERATIONS; i++) {
        //create a subset of the data
        ClusterDataSet *sampledSubset = data.randomSubset(getSampleSize());

        //create a subinstance of KMedoids and run it
        KMedoids *subcall = new KMedoids(*sampledSubset, k);
        subcall->findClusters();

        //copy medoids from subcall and assign all objects to clusters
        sampledSubset->sub2super(subcall->getMedoids().begin(), medoids.begin(), k);
        assignObjectsToClusters();
        
        //compute average dissimilarity and compare to best so far.
        dissimilarity = averageDissimilarity();
        if (dissimilarity < bestDissimilarity) {
            bestMedoids.swap(medoids);
            bestClustering.swap(clusterIds);
            bestDissimilarity = dissimilarity;
        }

        //clean up a bit before next iteration
        delete sampledSubset;
        delete subcall;
    }
    //POST: best medoids/clustering in bestMedoids and bestClustering

    //use the best clustering we found
    medoids.swap(bestMedoids);
    clusterIds.swap(bestClustering);
}
