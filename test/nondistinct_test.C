#include <iostream>
#include <iomanip>

#include "kmedoids.h"
#include "matrix_utils.h"
#include "point.h"
#include "Timer.h"

using namespace cluster;
using namespace std;


int main(int argc, char **argv) {
  size_t max_k = 10;
  if (argc > 1) {
    max_k = atoi(argv[1]);
  }

  // make a set of points with some duplicate elements
  std::vector<point> points;
  for (size_t i=0; i < 16; i++) {
    size_t type = i % 3;
    switch (type) {
    case 0:
      points.push_back(point(1,1));
      break;
    case 1:
      points.push_back(point(3,3));
      break;
    case 2:
      points.push_back(point(4,4));
      break;
    }
  }

  dissimilarity_matrix distance;
  build_dissimilarity_matrix(points, point_distance(), distance);

  kmedoids cluster;
  for (size_t k=1; k <= 10; k++) {
    cluster.pam(distance, k);

    vector<object_id> medoids(cluster.medoid_ids);
    sort(medoids.begin(), medoids.end());
    size_t num_unique = unique(medoids.begin(), medoids.end()) - medoids.begin();
    if (num_unique != cluster.medoid_ids.size()) {
      cerr << "Error: medoids are not distinct:" << endl;
      cerr << "[";
      copy(cluster.medoid_ids.begin(), cluster.medoid_ids.end(), 
           ostream_iterator<object_id>(cerr, " "));
      cerr << "]" << endl;

      cout << "FAILED" << endl;
      exit(1);
    }

    // ensure that each cluster has at least one element
    cluster_list clist;
    cluster.to_cluster_list(clist);
    for (size_t i=0; i < clist.size(); i++) {
      if (!clist[i].size()) {
        cerr << "Error: clustering contains empty clusters." << endl;
        cerr << cluster << endl;
        cout << "FAILED" << endl;
        exit(1);
      }
    }
  }

  cout << "PASSED" << endl;
  exit(0);
}
