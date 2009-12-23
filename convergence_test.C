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
  parse_points("( 1,  2) ( 2,  1) ( 5,  1) ( 5,  0) (10,  2) (16,  2) (17,  1) (23,  2)"
               "(23,  0) (24,  1) (31,  2) (30,  1) (40,  1) (40,  2) (40,  0) (50,  2)"
               "(49,  1) (51,  1) (61,  1) (62,  1) (73,  1) (86,  1) (86,  0) (85,  1)"
               "(100, 1) (100, 0) (115, 1) (115, 2) (116, 1) (131, 2) (148, 1) (148, 2)"
               "(147, 1) (149, 1) (166, 1) (166, 2) (166, 0) (165, 1) (185, 0) (186, 1)"
               "(205, 1) (205, 0) (226, 2) (248, 1) (247, 1) (271, 2) (271, 0) (295, 1)"
               "(295, 0) (320, 1) (320, 0) (319, 1) (321, 1) (346, 2) (372, 1) (374, 1)", 
               points);
  const object_id initial_medoids[] = {1, 14, 25, 27, 34, 42, 52, 55};

  kmedoids cluster;
  dissimilarity_matrix distance;
  build_dissimilarity_matrix(points, point_distance(), distance);

  cerr << "Starting PAM" << endl;
  cluster.pam(distance, 8, initial_medoids);
  cerr << "Finished."    << endl;
}
