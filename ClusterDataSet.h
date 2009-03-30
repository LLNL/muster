#ifndef CLUSTER_DATA_SET_H
#define CLUSTER_DATA_SET_H


#include "matrix.h"
#include "Dissimilarity.h"
#include "RNGenerator.h"


/**
 * Wrapper for a dissimilarity matrix with methods used by clustering
 * algortihms.
 * <p>
 * Contains methods for taking random subsets.
 * <p>
 * Class also provides facilities for taking a subset of a ClusterDataSet,
 * and for maintaining mappings from objects in the subset to objects
 * in the superset.
 */
class ClusterDataSet {
private:

    /** 
     * For random samples... allows clustering algos to 
     * randomly sample this data. This auto-seeds from clock.
     */
    mutable RNGenerator rnGenerator;
    
    /** Number of objects in the data set. */
    const unsigned numElements;

    /** n x n Distance matrix bt/w different objects. */
    const matrix<double> *distances;

    /**
     * Superset of this data set, if one exists.
     */
    const ClusterDataSet *const superSet;

    /** 
     * Mapping from indices in this data set to indices in superset. 
     * Size of this vector is n.
     */
    vector<unsigned> *const superSetMapping;


    /**
     * Internal constructor for creating subsets.
     */
    ClusterDataSet(
        const ClusterDataSet *superSet, 
        const matrix<double> *dissimilarity, 
        vector<unsigned> *superSetMapping, 
        int size
    );

    /**
     * Constructs a data set given a dissimilarity matrix and a size.  Used to 
     * create standalone data sets.
     */
    ClusterDataSet(const matrix<double> *dissimilarity, int size);

    
    //disallow direct assignment/copying
    ClusterDataSet(const ClusterDataSet& other);
    ClusterDataSet& operator=(ClusterDataSet& other);

public:

    /**
     * Given a matrix of objects and dissimilarity measure between them, builds
     * a ClusterDataSet by precomputing all distances between objects.
     */ 
    template <class T>
    static ClusterDataSet *buildFromObjects(const vector<T>& objects, const Dissimilarity<T>& diss) {
        matrix<double> *distances = new matrix<double>(objects.size(), objects.size());
        for (unsigned i=0; i < objects.size(); i++) {
            for (unsigned j=0; j < objects.size(); j++) {
                distances->get(i, j) = diss.getDissimilarity(objects[i], objects[j]);
            }
        }

        return new ClusterDataSet(distances, objects.size());
    }


    /**
     * Deallocates things.
     */
    ~ClusterDataSet();
    

    /**
     * Returns a reference to the superset of this data set, or
     * null if none.
     */
    const ClusterDataSet *getSuperSet() const;

    
    /**
     * Gets the index in the superset corresponding to the specified
     * index in this data set.
     */
    unsigned getSuperSetIndex(const unsigned i) const;

    /**
     * Translates n subset indices to superset indices.
     * @param src input iterator to read subset indices from
     * @param dest output iterator to write equivalent superset indices to
     * @param n number of indices to translate.
     */
    template<class InputIterator, class OutputIterator>
    void sub2super(InputIterator src, OutputIterator dest, unsigned n) {
        for (unsigned i=0; i < n; i++) {
            *dest = getSuperSetIndex(*src);
            src++;
            dest++;
        }
    }


    /**
     * Returns another ClusterDataSet referencing this data
     * for a randomly chosen set of 
     */
    ClusterDataSet* randomSubset(unsigned size) const;


    /**
     * Gets an index of a random object in this data set.
     */
    unsigned getRandomIndex() const {
        return rnGenerator(numElements);
    }


    /**
     * Gets the precomputed distance between objects i and j.
     * @param i index of an object in the data set
     * @param j index of an object in the data set
     */
    double distance(const unsigned i, const unsigned j) const {
        return distances->get(i,j);
    }

    
    /**
     * Returns the index into objects of the object nearest the
     * specified testObject, according to our distance matrix.
     * testObject and contents of objects are indices into the objects
     * array.
     *
     * @param objects array of object indices to search for nearness
     * @param testObject object to test against contents of objects
     * @param exclude optional index in objects to exclude.  
     *        Defaults to INT_MAX (which usually means include everything).
     */
    unsigned nearest(
        unsigned testObject, 
        vector<unsigned>& indexArray, 
        unsigned exclude = UINT_MAX
    ) const;


    /**
     * Returns the number of objects in this cluster data set.
     */
    unsigned size() const {
        return numElements;
    }


    void printDistancesOn(ostream& out = cout) {
        distances->printOn(out);
    }

};


#endif  //CLUSTER_DATA_SET_H
