#include <vector>
#include <iostream>
#include <cstdlib>
using namespace std;

#include "point.h"
#include "dissimilarity.h"
#include "kmedoids.h"
#include "color.h"
#include "point.h"
using namespace cluster;



/// Test using points instead of QGrams to make sure clustering 
/// algorithms work.
int main(int argc, char **argv) { 
  //vector of test points
  vector<point> points;
  
  size_t clusters = 5;
  if (argc > 1) {
    clusters = strtol(argv[1], NULL, 0);    
  }

  size_t max_k = clusters;
  if (argc > 2) {
    max_k = strtol(argv[2], NULL, 0);    
  }

  //put 5-pt crosses inthe vector, offset by increasing distances
  point ref(1,1);
  point stencil[] = {
    point( 0, 0),
    point( 0, 1), 
    point( 0,-1), 
    point(-1, 0), 
    point( 1, 0)
  };
  size_t stencil_size = sizeof(stencil) / sizeof(point);

  for (size_t i=0; i < clusters; i++) {
    for (size_t s=0; s < stencil_size; s++) {
      point p = ref + stencil[s];
      points.push_back(p);
    }
    ref += point(i+4, 0);
  }

  cout << "Testing with " << points.size() << " points for 1 to " << clusters << " clusters." << endl;

  dissimilarity_matrix dmatrix;
  build_dissimilarity_matrix(points, point_distance(), dmatrix);

  kmedoids km, clara, xkm, xclara;
  for (size_t k = 1; k <= max_k; k++) {
    km.pam(dmatrix, k);
    clara.clara(points, point_distance(), k);
    xkm.xpam(dmatrix, k, 2);
    xclara.xclara(points, point_distance(), k, 2);

    cout << "k: " << k << ", Mirkin distance: " << mirkin_distance(km, clara) << endl;


    draw("PAM", points, km);
    draw("CLARA", points, clara);
    draw("XPAM", points, xkm);
    draw("XCLARA", points, xclara);
    cout << endl;
  }
}
