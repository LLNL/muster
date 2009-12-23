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

  // make a set of points to cluster in a diagonal line from the origin.
  std::vector<point> points;
  for (size_t i=0; i < max_k; i++) {
    points.push_back(point(i,i));
  }
  
  dissimilarity_matrix distance;
  build_dissimilarity_matrix(points, point_distance(), distance);

  kmedoids cluster;
  for (size_t k=1; k <= max_k; k++) {
    Timer timer;
    cluster.pam(distance, k);
    timer.record("Cluster");
    
    cerr << k << setw(20) << timer["Cluster"]/1e9 << endl;
  }
}
