#include <iostream>
#include <iomanip>

#include "kmedoids.h"
#include "matrix_utils.h"
#include "point.h"
#include "Timer.h"

using namespace cluster;
using namespace std;

struct test {
  const char *points;
  size_t k;
};

int main(int argc, char **argv) {
  const test tests[] = {
    {"(1,  1) (1,  2) (0,  1) (2,  1) (19,  1) (19,  2) (19,  0) (18,  1) (38,  0)"
     "(39, 1) (58, 1) (58, 0) (79, 2) (101, 1) (100, 1) (124, 2) (124, 0)",
     3
     },
    {"( 1, 1) ( 1, 2) ( 0, 1) ( 2, 1) ( 5, 1) ( 5, 2) ( 5, 0) ( 4, 1) ( 6, 1) (10, 1)"
     "(10, 2) (10, 0) ( 9, 1) (16, 1) (16, 2) (16, 0) (15, 1) (17, 1) (23, 2) (23, 0)"
     "(31, 1) (31, 0) (30, 1) (32, 1) (40, 1) (40, 2) (39, 1) (41, 1) (50, 1) (50, 0)"
     "(49, 1) (51, 1) (61, 2) (61, 0) (62, 1) (72, 1) (74, 1) (86, 1) (86, 2) (86, 0)"
     "(87, 1) (100, 1) (100, 2) (100, 0) (99, 1) (115, 1) (115, 2) (115, 0)",
     4},
    {"( 1, 1) ( 1, 2) ( 1, 0) ( 0, 1) ( 2, 1) ( 5, 1) ( 5, 2) ( 5, 0) ( 6, 1) (10, 1) (10, 2) (10, 0)"
     "( 9, 1) (16, 1) (16, 2) (16, 0) (15, 1) (17, 1) (23, 2) (22, 1) (24, 1) (31, 1) (31, 2) (40, 0)"
     "(39, 1) (50, 1) (50, 2) (51, 1) (61, 1) (61, 2) (61, 0) (60, 1) (73, 1) (73, 2) (73, 0) (72, 1)"
     " (86, 1) (100, 1) (100, 2) (100, 0) (101, 1) (115, 2) (115, 0) (114, 1)",
     6},
    {"( 1, 1) ( 1, 2) ( 1, 0) ( 0, 1) ( 2, 1) ( 5, 1) ( 5, 2) ( 5, 0) ( 4, 1) ( 6, 1) (10, 1) (10, 2)"
     "(10, 0) ( 9, 1) (11, 1) (16, 1) (16, 2) (16, 0) (15, 1) (17, 1) (23, 1) (23, 2) (23, 0) (22, 1)"
     "(24, 1) (31, 1) (31, 2) (31, 0) (30, 1) (32, 1) (40, 1) (40, 2) (40, 0) (39, 1) (41, 1) (50, 1)"
     "(50, 2) (50, 0) (49, 1) (51, 1) (61, 1) (61, 2) (61, 0) (60, 1) (62, 1) (73, 1) (73, 2) (73, 0)"
     "(72, 1) (74, 1) (86, 1) (86, 2) (86, 0) (85, 1) (87, 1) (100, 1) (100, 2) (100, 0) (99, 1) "
     "(101, 1) (115, 1) (115, 2) (115, 0) (114, 1)", 
     6}
  };
  const size_t num_tests = sizeof(tests) / sizeof(test);
  

  for (size_t i=0; i < num_tests; i++) {
    vector<point> points;
    parse_points(tests[i].points, points);
    const size_t k = tests[i].k;

    kmedoids cluster;
    dissimilarity_matrix distance;
    build_dissimilarity_matrix(points, point_distance(), distance);
    
    cerr << "Starting Test " << i << " with k=" << k <<endl;
    cluster.pam(distance, k);
    cerr << "Finished."    << endl;

    ostringstream label;
    label << "Trial " << i;
    draw(label.str(), points, cluster);

  }
}
