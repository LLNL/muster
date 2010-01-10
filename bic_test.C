#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iterator>
#include <sys/time.h>

#include "point.h"
#include "io_utils.h"
#include "wavelet.h"
#include "kmedoids.h"
#include "par_kmedoids.h"
#include "bic.h"
#include "matrix_utils.h"
#include "gaussian.h"

using namespace cluster;
using namespace std;
using boost::numeric::ublas::matrix;


void usage() {
  cerr << "Usage: par-bic-test [-htvi] [-n num-points] [-c clusters] [-s scale]" << endl;
  cerr << "  Parallel test case for clustering points and evaluating their BIC scores." << endl;
  cerr << "Options:" << endl;
  cerr << "  -h         Show this message." << endl;
  cerr << "  -t         Output timing info to file." << endl;
  cerr << "  -v         Validate with sequential clustering and output Mirkin distance." << endl;
  cerr << "  -n         Number of points generated." << endl;
  cerr << "               Default is 128." << endl;
  cerr << "  -c         Max number of clusters to search for." << endl;
  cerr << "               Default is 1." << endl;
  cerr << "  -s         Scale values by this amount.." << endl;
  exit(1);
}

bool timing = false;
bool validate = false;
size_t num_points = 128;
size_t max_clusters = 10;
const size_t dimensions = 1;

double scale = 1.0;

/// Uses getopt to read in arguments.
void get_args(int *argc, char ***argv, int rank) {
  int c;
  char *err;

  while ((c = getopt(*argc, *argv, "htvl:n:s:c:x:y:")) != -1) {
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
    case 'n':
      num_points = strtol(optarg, &err, 0);
      if (*err) usage();
      break;
    case 'c':
      max_clusters = strtol(optarg, &err, 0);
      if (*err) usage();
      break;
    case 's':
      scale = strtod(optarg, &err);
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

const double xscale = 2.5;  // even out based on screen character size.
const double yscale = 1.0;

vector<point> centroids;
vector<point> points;
dissimilarity_matrix dissimilarity;

/// Callback for printing out clusterings and their BIC scores.
static void print_cluster_info(const cluster::partition& km, double bic) {
    cout << "k:       " << km.num_clusters() << endl;
    cout << "BIC:     " << bic << endl;
    /*
    cout << "medoids: ";
    for (size_t i=0; i < km.medoid_ids.size(); i++) {
      cout << points[km.medoid_ids[i]] << " ";
    }
    cout << endl;
    */
    cout << "old_bic: " << old_bic(km, matrix_distance(dissimilarity), dimensions) << endl;
    //cout << "D:       " << total_dissimilarity(km, matrix_distance(dissimilarity)) << endl;
    //cout << "D2:      " << total_squared_dissimilarity(km, matrix_distance(dissimilarity)) << endl;

    //cout << km << endl;

    draw("Clustering", points, km);
    cout << endl;
}


int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  get_args(&argc, &argv, rank);

  // make sure max_clusters is valid.
  max_clusters = min(max_clusters, size * num_points);

  const double stddev = 0.3; // size of clusters

  centroids.push_back(point(2, 4));
  centroids.push_back(point(6, 1));
  centroids.push_back(point(15, 2));
  centroids.push_back(point(4, 1));
  centroids.push_back(point(4, 8));
  centroids.push_back(point(10, 7));

  vector<gaussian_generator_2d> clusters;
  for (size_t i=0; i < centroids.size(); i++) {
    clusters.push_back(
      gaussian_generator_2d(centroids[i].x  * scale, centroids[i].y * scale, stddev * scale, xscale, yscale));
  }

  for (size_t r=0; r < num_points; r++) {
    size_t type = (rank + r) % clusters.size();
    point p;
    do {
      p = clusters[type].next_point();
    } while (p.x < 0 || p.y < 0);

    points.push_back(p);
  }


  build_dissimilarity_matrix(points, point_distance(), dissimilarity);
  
  kmedoids km;
  km.set_xcallback(print_cluster_info);
  //double best_bic = km.xpam(dissimilarity, max_clusters, dimensions);
  double best_bic = km.xclara(points, point_distance(), max_clusters, dimensions);
  //double best_bic = 0;
  //km.clara(points, point_distance(), max_clusters, dimensions);


  cout << "========== Best ==========" << endl;
  print_cluster_info(km, best_bic);
}
