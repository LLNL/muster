

#include <vector>
#include <iostream>
using namespace std;

#include "Point.h"
#include "Dissimilarity.h"
#include "ClusterDataSet.h"
#include "KMedoids.h"
#include "CLARA.h"


struct PointDistance : public Dissimilarity<Point*> {
    virtual double getDissimilarity(Point* left, Point* right) const {
        return left->distance(right);
    }
};

/**
 * Test using Points instead of QGrams to make sure clustering 
 * algorithms work.
 */
int main(int argc, char **argv) { 

    //vector of test points
    vector<Point*> points;

    //put some clusters in the vector.  These are 5-pt crosses
    points.push_back(new Point(0,10));
    points.push_back(new Point(0,5));
    points.push_back(new Point(0,0));


    points.push_back(new Point(0,1));
    points.push_back(new Point(1,0));
    points.push_back(new Point(1,1));
    points.push_back(new Point(1,2));
    points.push_back(new Point(2,1));
    
    points.push_back(new Point(5,1));
    points.push_back(new Point(6,0));
    points.push_back(new Point(6,1));
    points.push_back(new Point(6,2));
    points.push_back(new Point(7,1));
    
    points.push_back(new Point(0,8));
    points.push_back(new Point(1,7));
    points.push_back(new Point(1,8));
    points.push_back(new Point(1,9));
    points.push_back(new Point(2,8));

    PointDistance diss;
    ClusterDataSet *data = ClusterDataSet::buildFromObjects(points, diss);

    cerr << "Constrcted DataSet." << endl;

    /*
    for (int i=0; i < 15; i++) {
        for (int j=0; j < 15; j++) {
            cout << data.distance(i,j) << " == " << data.distance(j,i) << endl;
        }

    }
    */

    for (int k = 1; k < 18; k++) {
        KMedoids kmedoids(*data, k);
        cerr << "Finding Clusters with k=" << k << endl;;
        kmedoids.findClusters();
        cerr << "done." << endl;
        kmedoids.printClustering();
    }

    CLARA clara(*data, 2);
    cerr << "Finding Clusters with CLARA...";
    clara.findClusters();
    cerr << "done." << endl;
    clara.printClustering();
}
