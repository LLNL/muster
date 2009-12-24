#include <iostream>
#include <iomanip>

#include "kmedoids.h"
#include "matrix_utils.h"
#include "point.h"
#include "Timer.h"

using namespace cluster;
using namespace std;


int main(int argc, char **argv) {
  vector<point> points;
  parse_points("(1, 1) (1, 2)"
               "(0,  1) (2,  1) (19, 1) (19,  2) (19,  0) (18,  1) (38,  0) (39,  1)"
               "(58, 1) (58, 0) (79, 2) (101, 1) (100, 1) (124, 2) (124, 0)",
               points);
  const object_id initial_medoids[] = {4, 12, 8};//22, 25};
  const size_t k = 3;

  kmedoids cluster;
  dissimilarity_matrix distance;
  build_dissimilarity_matrix(points, point_distance(), distance);

  cerr << "Starting PAM" << endl;
  cluster.pam(distance, k, initial_medoids);
  cerr << "Finished."    << endl;

  draw("Converged", points, cluster);
}
