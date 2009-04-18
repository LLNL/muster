#include <set>
#include <vector>
#include <iostream>
#include <stdexcept>
using namespace std;

#include "ClusterDataSet.h"


ClusterDataSet::ClusterDataSet(const matrix<double> *dissimilarity, int size):
    numElements(size),
    distances(dissimilarity),
    superSet(NULL),
    superSetMapping(NULL)    
{
    ; //nothing to do
}


ClusterDataSet::ClusterDataSet(
    const ClusterDataSet *super, const matrix<double> *dissimilarity, 
    vector<unsigned> *mapping, int size
):
    numElements(size),
    distances(dissimilarity),
    superSet(super),
    superSetMapping(mapping)
{
    ; // nothing to do.
}

/**
 * Frees up the distance matrix.
 */
ClusterDataSet::~ClusterDataSet() {
    delete distances;
    if (superSetMapping) delete superSetMapping;
}




const ClusterDataSet *ClusterDataSet::getSuperSet() const {
    return superSet;
}


unsigned ClusterDataSet::getSuperSetIndex(const unsigned i) const {
    return (superSetMapping == NULL) ? i : (*superSetMapping)[i];
}


/**
 * Returns another ClusterDataSet containing a randomly chosen
 * set of <size> of our objects.
 */
ClusterDataSet* ClusterDataSet::randomSubset(unsigned size) const {
    if (size > numElements) {
        throw invalid_argument("size passed to randomSubset must be <= size()");
    }

    //Mapping from subset to superset indices, for the subset.
    vector<unsigned> *mapping = new vector<unsigned>();

    // This is Knuth's algorithm R for taking a sample of indices from
    // 0 to numElements.  We sample size indices from this (the superset)
    // and put them in the subset's mapping.
    insert_iterator<vector<unsigned> > out = inserter(*mapping, mapping->begin());
    unsigned first = 0;
    unsigned remaining = numElements;
    unsigned m = size;
    
    while (m > 0) {
      if (rnGenerator(remaining) < (int)m) {
            *out = first;
            ++out;
            --m;
	    }
        --remaining;
        ++first;
	}

    //use a submatrix of our precomputed distance matrix into the new data set
    matrix<double> *newDistances = distances->getSubmatrix(
        mapping->begin(), mapping->end(),
        mapping->begin(), mapping->end()
    );

    return new ClusterDataSet(this, newDistances, mapping, size);
}


unsigned ClusterDataSet::nearest(
    unsigned testObject, 
    vector<unsigned>& objects, 
    unsigned exclude
) const {
    double minDistance = DISSIMILARITY_MAX;
    unsigned nearest = 0;
    
    for (unsigned i=0; i < objects.size(); i++) {
        if (i == exclude) continue;

        double curDistance = distance(testObject, objects[i]);
        if (curDistance < minDistance) {
            minDistance = curDistance;
            nearest = i;
        }
    }
    return nearest;
}





