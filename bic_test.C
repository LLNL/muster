#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <sys/time.h>

#include "point.h"
#include "io_utils.h"
#include "wavelet.h"
#include "kmedoids.h"
#include "par_kmedoids.h"
#include "bic.h"
#include "matrix_utils.h"
#include "MersenneTwister.h"

using namespace cluster;
using namespace std;
using boost::numeric::ublas::matrix;

void usage() {
  cerr << "Usage: par-bic-test [-htvi] [-p num-points] [-c clusters]" << endl;
  cerr << "  Parallel test case for clustering points and evaluating their BIC scores." << endl;
  cerr << "Options:" << endl;
  cerr << "  -h         Show this message." << endl;
  cerr << "  -t         Output timing info to file." << endl;
  cerr << "  -v         Validate with sequential clustering and output Mirkin distance." << endl;
  cerr << "  -p         Number of points generated per process." << endl;
  cerr << "               Default is 1." << endl;
  cerr << "  -c         Max number of clusters to search for." << endl;
  cerr << "               Default is 1." << endl;
  exit(1);
}

bool timing = false;
bool validate = false;
size_t num_points = 128;
size_t max_clusters = 10;
const size_t dimensions = 2;

const double xscale = 4.0;
const double yscale = 1.0;

/// Uses getopt to read in arguments.
void get_args(int *argc, char ***argv, int rank) {
  int c;
  char *err;

  while ((c = getopt(*argc, *argv, "htvl:n:s:c:")) != -1) {
    switch (c) {
    case 'h':
      if (rank == 0) usage();
      exit(1);
      break;
    case 't':
      timing = true;
      break;
    case 'v':
      validate = true;
      break;
    case 'p':
      num_points = strtol(optarg, &err, 0);
      if (*err) usage();
      break;
    case 'c':
      max_clusters = strtol(optarg, &err, 0);
      if (*err) usage();
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


static point create_point(MTRand& rand, int x, int y, int xr, int yr, double xscale, double yscale) {
  xr = (int)(xr * xscale);
  yr = (int)(yr * yscale);
  x  = (int)(x  * xscale);
  y  = (int)(y  * yscale);

  int px = x + rand.randInt(xr*2) - xr;
  int py = y + rand.randInt(yr*2) - yr;
  return point(px,py);
}

vector<point> points;
dissimilarity_matrix dissimilarity;

/// Callback for printing out clusterings and their BIC scores.
static void callback(const cluster::partition& km, double bic) {
    cout << endl;
    cout << "k:   " << km.num_clusters() << endl;
    cout << "BIC: " << bic << endl;
    cout << "D:   " << total_dissimilarity(km, matrix_distance(dissimilarity)) << endl;
    cout << "D2:  " << total_squared_dissimilarity(km, matrix_distance(dissimilarity)) << endl;

    //cout << km << endl;

    draw("Clustering", points, km);
}


int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  get_args(&argc, &argv, rank);

  // make sure max_clusters is valid.
  max_clusters = min(max_clusters, size * num_points);
  
  uint32_t seed = 0;
  struct timeval time;
  gettimeofday(&time, 0);
  seed = time.tv_sec * time.tv_usec;

  MTRand rand(seed + rank);  // Seed points differently on each rank.


  for (size_t r=0; r < num_points; r++) {
    size_t type = (rank + r) % 3;
    switch (type) {
    case 0:
      points.push_back(create_point(rand, 2, 4, 2, 2, xscale, yscale));
      break;
    case 1:
      points.push_back(create_point(rand, 12, 12, 2, 2, xscale, yscale));
      break;
    case 2:
      points.push_back(create_point(rand, 12, 3, 2, 2, xscale, yscale));
      break;
    }
  }


  build_dissimilarity_matrix(points, point_distance(), dissimilarity);
  
  kmedoids km;
  km.set_xpam_callback(callback);
  double best_bic = km.xpam(dissimilarity, max_clusters, dimensions);


  cout << "========== Best ==========" << endl;
  cout << "k:   " << km.num_clusters() << endl;
  cout << "BIC: " << best_bic << endl;
  cout << "D:   " << total_dissimilarity(km, matrix_distance(dissimilarity)) << endl;
  cout << "D2:  " << total_squared_dissimilarity(km, matrix_distance(dissimilarity)) << endl;
  //cout << km << endl;
  draw("Clustering", points, km);


  ofstream point_file("points.xy");
  copy(points.begin(), points.end(), ostream_iterator<point>(point_file, "\n"));

}
