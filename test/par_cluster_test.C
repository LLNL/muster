#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

#include "point.h"
#include "bic.h"
#include "par_kmedoids.h"

using namespace cluster;
using namespace std;

void usage() {
  cerr << "Usage: par-cluster-test [-hxv] [-n num-points] [-k clusters] [-i initial-size] [-r reps]" << endl;
  cerr << "  Compare parallel clustering with sequential clustering." << endl;
  cerr << "Options:" << endl;
  cerr << "  -h         Show this message." << endl;
  cerr << "  -x         Use BIC-scored versions of PAM and CLARA." << endl;
  cerr << "  -v         Verbose mode.  Draws actual clusterings and outputs timings." << endl;
  cerr << "  -n         Number of points per process." << endl;
  cerr << "               Default is 1." << endl;
  cerr << "  -k         Max number of clusters to search for." << endl;
  cerr << "               Default is number of processes * points per process." << endl;
  cerr << "  -i         Initial sample size in clara (before 2*k is added)." << endl;
  cerr << "               Default is 40." << endl;
  cerr << "  -r         Number of repeated trials per k in clara." << endl;
  cerr << "               Default is 5." << endl;
  exit(1);
}

size_t objects_per_process = 1;
size_t num_clusters = 0;
size_t init_size = 40;
size_t max_reps = 5;
bool use_bic = false;
bool verbose = false;

/// Uses getopt to read in arguments.
void get_args(int *argc, char ***argv, int rank) {
  int c;
  char *err;

  while ((c = getopt(*argc, *argv, "hxn:i:r:k:v")) != -1) {
    switch (c) {
    case 'h':
      if (rank == 0) usage();
      exit(1);
      break;
    case 'x':
      use_bic = true;
      break;
    case 'n':
      objects_per_process = strtol(optarg, &err, 0);
      if (*err) usage();
      break;
    case 'i':
      init_size = strtol(optarg, &err, 0);
     if (*err) usage();
      break;
    case 'r':
      max_reps = strtol(optarg, &err, 0);
      if (*err) usage();
      break;
    case 'k':
      num_clusters = strtol(optarg, &err, 0);
      if (*err) usage();
      break;
    case 'v':
      verbose = true;
      break;
    default:
      if (rank == 0) usage();
      exit(1);
      break;
    }
  }

  // adjust params
  *argc -= optind;
  *argv += optind;
}


vector<point> points;
kmedoids km;
par_kmedoids parkm(MPI_COMM_WORLD);


void print_cluster_info(const cluster::partition& gathered, const dissimilarity_matrix& distance) {
  cout << "seq k: " << km.medoid_ids.size()
       << ", par k: " << gathered.medoid_ids.size()
       << ", Mirkin distance: " << setprecision(3) << mirkin_distance(km, gathered)
       << endl;
  
  if (verbose) {
    ostringstream pam_msg;
    pam_msg << "PAM"
            << ", " << km.medoid_ids.size() << " clusters"
            << ", Avg. dissimilarity: " << km.average_dissimilarity()
            << ", BIC: " << bic(km, matrix_distance(distance), 2);
    
    ostringstream parkm_msg;
    parkm_msg << "Parallel CLARA" 
              << ", " << gathered.medoid_ids.size() << " clusters"
              << ", Avg. dissimilarity: " << parkm.average_dissimilarity()
              << ", BIC: " << parkm.bic_score();
    
    draw(pam_msg.str(), points, km);
    draw(parkm_msg.str(), points, gathered);
    cout << endl;
    parkm.get_timer().write();
    cout << endl;
  }
}


/// Test using points instead of QGrams to make sure clustering 
/// algorithms work.
int main(int argc, char **argv) { 
  MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  get_args(&argc, &argv, rank);
  
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


  size_t num_objects = size * objects_per_process;
  if (num_clusters == 0) {
    num_clusters = num_objects / stencil_size + 5;
  }

  num_clusters = min(num_clusters, num_objects);

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

  parkm.set_init_size(init_size);
  parkm.set_max_reps(max_reps);

  for (size_t k=1; k <= num_clusters; k++) {
    vector<point> medoids;
    cluster::partition gathered;


    km.pam(distance, k);
    parkm.clara(my_points, point_distance(), k, &medoids);
    parkm.gather(gathered);

    if (rank == 0) {
      cout << "Max k = " << k << endl;
      cout << "  No BIC:    ";
      print_cluster_info(gathered, distance);
    }

    km.xpam(distance, k, 2);
    parkm.xclara(my_points, point_distance(), k, 2, &medoids);
    parkm.gather(gathered);

    if (rank == 0) {
      cout << "  Using BIC: ";
      print_cluster_info(gathered, distance);
      cout << endl;
    }
  }
}

