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
  MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

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

  size_t objects_per_process = 10;
  if (argc > 1) {
    objects_per_process = strtol(argv[1], NULL, 0);    
  }

  size_t num_objects = size * objects_per_process;
  size_t max_clusters = num_objects / stencil_size + 5; // go 5 over to test BIC
  if (argc > 2) {
    max_clusters = strtol(argv[2], NULL, 0);
  }

  vector<point> my_points;
  for (size_t i=0; points.size() < num_objects; i++) {
    for (size_t s=0; s < stencil_size && points.size() < num_objects; s++) {
      point p = ref + stencil[s];
      if (rank == (int)(points.size() / objects_per_process)) {
        my_points.push_back(p);
      }
      points.push_back(p);
    }
    ref += point(i+4, 0);
  }

  dissimilarity_matrix distance;
  build_dissimilarity_matrix(points, point_distance(), distance);

  kmedoids km;
  par_kmedoids parkm;
  for (size_t k=1; k <= max_clusters; k++) {
    km.pam(distance, k);

    vector<point> medoids;
    parkm.xclara(my_points, point_distance(), k, 2, &medoids);
    cluster::partition local_partition;
    parkm.gather(local_partition);

    if (rank == 0) {
      cout << "k: " << k 
           << ", Mirkin distance: " << setprecision(3) << mirkin_distance(km, local_partition) 
           << endl;

      ostringstream pam_msg;
      pam_msg << "PAM"
              << ", " << km.medoid_ids.size() << " clusters"
              << ", Avg. dissimilarity: " << km.average_dissimilarity()
              << ", BIC: " << bic(km, matrix_distance(distance), 2);

      ostringstream parkm_msg;
      parkm_msg << "Parallel CLARA" 
                << ", " << local_partition.medoid_ids.size() << " clusters"
                << ", Avg. dissimilarity: " << parkm.average_dissimilarity()
                << ", BIC: " << parkm.bic_score();

      draw(pam_msg.str(), points, km);
      draw(parkm_msg.str(), points, local_partition);
      cout << endl;
      parkm.get_timer().write();
      cout << endl;
    }
  }
}
