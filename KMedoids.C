#include "KMedoids.h"

#include <algorithm>
#include <stdexcept>
#include <cassert>


KMedoids::KMedoids(
    ClusterDataSet& _data,
    const unsigned _k
) :
    data(_data),
    k(_k)
{
    medoids.resize(k);
    clusterIds.resize(data.size());

    if (k > data.size()) {
        throw logic_error("Attempt to instantiate KMedoids with more clusters than data.");
    }
}


KMedoids::~KMedoids() {
    //nothing necessary.
}


void KMedoids::assignRandomMedoids() {
    //find k unique objects randomly
    set<unsigned> uniqueMedoids;
    while (uniqueMedoids.size() < k) {
        uniqueMedoids.insert(data.getRandomIndex());
    }

    //copy unique medoids into the medoids vector.
    copy(uniqueMedoids.begin(), uniqueMedoids.end(), medoids.begin());
    assert(medoids.size() == k);
    
    //point the medoid objects at the right clusters
    for (medoid_id mi = 0; mi < k; mi++) {
        clusterIds[medoids[mi]] = mi;
    }
}


bool KMedoids::isMedoid(unsigned objectId) {
    //check if the object is its own cluster's medoid.
    return medoids[clusterIds[objectId]] == objectId;
}


/**
 * Goes through all objects and assigns each to the 
 * cluster with the closest medoid.
 */
void KMedoids::assignObjectsToClusters() {
    for (object_id i=0; i < data.size(); i++) {
        if (!isMedoid(i)) {
            clusterIds[i] = data.nearest(i, medoids);
        }
    }
}


/**
 * Finds the cost w/respect to object j of swapping object oh with medoid i.
 */
double KMedoids::cost(medoid_id mi, object_id oh, object_id oj) {
    //oi is the object id of medoid i, mi is the medoid id.
    object_id oi = medoids[mi];

    double cost;
    if (clusterIds[oj] == mi) {
        medoid_id mj2 = data.nearest(oj, medoids, mi);
        object_id oj2 = medoids[mj2];
        if (data.distance(oj, oj2) < data.distance(oj, oh)) {
            cost = data.distance(oj, oj2) - data.distance(oj, oi);
        } else {
            cost = data.distance(oj, oh) - data.distance(oj, oi);            
        }

    } else {
        medoid_id mj2 = clusterIds[oj];
        object_id oj2 = medoids[mj2];
        if (data.distance(oj, oj2) < data.distance(oj, oh)) {
            cost = 0.0;
        } else {
            cost = data.distance(oj, oh) - data.distance(oj, oj2);
        }
    }

    return cost;
}



/**
 * Total cost of swapping object h with medoid i.
 * Sums costs of this exchagne for all objects j.
 */
double KMedoids::totalCost(medoid_id i, object_id h) {
    double sum =0;

    //pair each medoid w/all non-medoids and see what 
    for (object_id j = 0; j < data.size(); j++) {
        //skip medoids and self
        if (isMedoid(j) || j == h) continue;  

        //add cost of swapping object h with object j
        sum += cost(i, h, j);
    }
    return sum;
}



void KMedoids::findClusters() {
    //randomly pick initial medoids
    assignRandomMedoids();
    
    if (k == 1) {
        //just bail out here if they only want one cluster.
        assignObjectsToClusters();

    //otherwise, find the clusters
    } else while (true) {
        //put objects in cluster w/nearest medoid
        assignObjectsToClusters();

        //vars to keep track of minimum
        double minTotalCost = DISSIMILARITY_MAX;
        medoid_id minMedoid = 0;
        object_id minObject = 0;

        //iterate over each medoid
        for (medoid_id i=0; i < k; i++) {
            //iterate over all non-medoid objects
            for (object_id h = 0; h < data.size(); h++) {
                if (isMedoid(h)) continue;

                //see if the total cost of swapping i & h was less than min
                double curCost = totalCost(i,h);
                if (curCost < minTotalCost) {
                    minTotalCost = curCost;
                    minMedoid = i;
                    minObject = h;
                }
            }
        }

        if (minTotalCost < 0) {
            //replace a medoid if it gains us something
            medoids[minMedoid] = minObject;
            clusterIds[minObject] = minMedoid;

        } else {
            //otherwise we're done swapping, so stop.
            break;
        }
    }
    
    //associate each object w/the correct cluster one last time.
    assignObjectsToClusters();
}



const ClusterDataSet& KMedoids::getData() {
    return data;
}

const vector<unsigned>& KMedoids::getMedoids() {
    return medoids;
}


void KMedoids::getSuperSetMedoids(vector<unsigned>& dest) {
    unsigned index = 0;
    for (vector<unsigned>::iterator i = medoids.begin(); i != medoids.end(); i++) {
        dest[index++] = data.getSuperSetIndex(*i);
    }
}


KMedoids::clusterList *KMedoids::getClustering() {
    //make a new list of empty clusters
    cluster c;
    clusterList *clusters = new clusterList(k, c);

    //insert each object into the appropriate cluster.
    for (unsigned object=0; object < data.size(); object++) {
        (*clusters)[clusterIds[object]].insert(object);
    }
    return clusters;
}


void KMedoids::printClustering(ostream& out) {
    out << "id\tmembers" << endl;
    
    clusterList &clusters = *getClustering();
    for (unsigned i=0; i < clusters.size(); i++) {
        out << i << "\t";
        cluster& c = clusters[i];
        for (cIterator object=c.begin(); object != c.end(); object++) {
            out << *object << " ";
        }
        out << endl;
    }
}


void KMedoids::printMedoids(ostream& out) {
    for (medoid_id i=0; i < medoids.size(); i++) {
        cout << i << " ";
    }
    cout << endl;
    for (medoid_id i=0; i < medoids.size(); i++) {
        cout << medoids[i] << " ";
    }
    cout << endl;
}

void KMedoids::printClusterIds(ostream& out) {
    for (object_id i=0; i < clusterIds.size(); i++) {
        cout << "    " << i << " -> " << clusterIds[i] << endl;
    }
}
