#include <iostream>
#include <iomanip>

#include "kmedoids.h"
#include "matrix_utils.h"
#include "point.h"
#include "Timer.h"

using namespace cluster;
using namespace std;


int main(int argc, char **argv) {
  /*
  parse_points("( 1,  2) ( 2,  1) ( 5,  1) ( 5,  0) (10,  2) (16,  2) (17,  1) (23,  2)"
               "(23,  0) (24,  1) (31,  2) (30,  1) (40,  1) (40,  2) (40,  0) (50,  2)"
               "(49,  1) (51,  1) (61,  1) (62,  1) (73,  1) (86,  1) (86,  0) (85,  1)"
               "(100, 1) (100, 0) (115, 1) (115, 2) (116, 1) (131, 2) (148, 1) (148, 2)"
               "(147, 1) (149, 1) (166, 1) (166, 2) (166, 0) (165, 1) (185, 0) (186, 1)"
               "(205, 1) (205, 0) (226, 2) (248, 1) (247, 1) (271, 2) (271, 0) (295, 1)"
               "(295, 0) (320, 1) (320, 0) (319, 1) (321, 1) (346, 2) (372, 1) (374, 1)", 
               points);
  const object_id initial_medoids[] = {1, 14, 25, 27, 34, 42, 52, 55};
  */

  vector<point> points;
  parse_points("(1, 1) (1, 2)"
               "(0,  1) (2,  1) (19, 1) (19,  2) (19,  0) (18,  1) (38,  0) (39,  1)"
               "(58, 1) (58, 0) (79, 2) (101, 1) (100, 1) (124, 2) (124, 0)",
               points);
  const object_id initial_medoids[] = {4, 12, 8};//22, 25};
  const size_t k = 3;

  kmedoids cluster;
  cluster.points = &points;

  dissimilarity_matrix distance;
  build_dissimilarity_matrix(points, point_distance(), distance);

  cerr << "Starting PAM" << endl;
  cluster.pam(distance, k, initial_medoids);
  cerr << "Finished."    << endl;
}
