#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

#include "point.h"
#include "bic.h"
#include "par_kmedoids.h"

using namespace cluster;
using namespace std;


/// Test using points instead of QGrams to make sure clustering 
/// algorithms work.
int main(int argc, char **argv) { 
  //put 5-pt crosses inthe vector, offset by increasing distances
  vector<point> points;
  point ref(1,1);
  point stencil[] = {
    point( 0, 0),
    point( 0, 1), 
    point( 0,-1), 
    point(-1, 0), 
    point( 1, 0)
  };
  size_t stencil_size = sizeof(stencil) / sizeof(point);

  size_t num_objects = 64;
  if (argc > 1) {
    num_objects = strtol(argv[1], NULL, 0);
  }

  size_t max_clusters = num_objects / stencil_size + 5; // go 5 over to test BIC

  for (size_t i=0; points.size() < num_objects; i++) {
    for (size_t s=0; s < stencil_size && points.size() < num_objects; s++) {
      point p = ref + stencil[s];
      points.push_back(p);
    }
    ref += point(i+4, 0);
  }

  cerr << "num_objects  = " << num_objects << endl;
  cerr << "max_clusters = " << max_clusters << endl;

  cerr << points.size() << " points: ";
  copy(points.begin(), points.end(), ostream_iterator<point>(cerr, " "));
  cerr << endl;

  dissimilarity_matrix distance;
  build_dissimilarity_matrix(points, point_distance(), distance);

  cerr << "matrix size is " << distance.size1() << "x" << distance.size2() << endl;

  kmedoids km;
  for (size_t k=1; k <= max_clusters; k++) {
    km.pam(distance, k);

    ostringstream pam_msg;
    pam_msg << "PAM"
            << ", " << km.medoid_ids.size() << " clusters"
            << ", Avg. dissimilarity: " << km.average_dissimilarity()
            << ", BIC: " << bic(km, matrix_distance(distance), 2);
    
    draw(pam_msg.str(), points, km);
    cout << endl;
  }
}
