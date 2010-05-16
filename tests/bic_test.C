#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iterator>
#include <sys/time.h>

#include "point.h"
#include "kmedoids.h"
#include "bic.h"
#include "spherical_clustering_generator.h"

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
  cerr << "  -s         Scale values by this amount." << endl;
  cerr << "  -d         Verbose debug output." << endl;
  exit(1);
}

bool debug = false;
bool timing = false;
bool validate = false;
size_t num_points = 128;
size_t max_clusters = 10;
const size_t dimensions = 2;

double scale = 1.0;

/// Uses getopt to read in arguments.
void get_args(int *argc, char ***argv) {
  int c;
  char *err;

  while ((c = getopt(*argc, *argv, "htdvl:n:s:c:x:y:")) != -1) {
    switch (c) {
    case 'h':
      usage();
      exit(1);
      break;
    case 't':
      timing = true;
      break;
    case 'd':
      debug = true;
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
      usage();
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

vector<point> points;
dissimilarity_matrix dissimilarity;


/// Callback for printing out clusterings and their BIC scores.
static void print_cluster_info(const cluster::partition& km, double bic) {
  cout << "k:        " << km.num_clusters() << endl;
  cout << "BIC:      " << bic << endl;
  cout << "real bic: " << cluster::bic(km, matrix_distance(dissimilarity), dimensions) << endl;
  cout << "old bic:  " << old_bic(km, matrix_distance(dissimilarity), dimensions) << endl;  
  if (debug) {
    cout << "medoids: ";
    for (size_t i=0; i < km.medoid_ids.size(); i++) {
      cout << points[km.medoid_ids[i]] << " ";
    }
    cout << endl;
  }

  if (debug) {
    cout << "D:       " << total_dissimilarity(km, matrix_distance(dissimilarity)) << endl;
    cout << "D2:      " << total_squared_dissimilarity(km, matrix_distance(dissimilarity)) << endl;
    cout << km << endl;
  }

  draw("Clustering", points, km);
  cout << endl;
}


static void print_cluster_info_if_6(const cluster::partition& km, double bic) {
  if (km.num_clusters() == 6) {
    print_cluster_info(km, bic);
  }
}


int main(int argc, char **argv) {
  get_args(&argc, &argv);

  // make sure max_clusters is valid.
  max_clusters = min(max_clusters, num_points);

  spherical_clustering_generator cgen;
  cgen.set_default_stddev(0.3);
  cgen.set_xscale(scale * xscale);
  cgen.set_yscale(scale * yscale);

  cgen.add_cluster(point(2,4));
  cgen.add_cluster(point(6, 1));
  cgen.add_cluster(point(15, 2));
  cgen.add_cluster(point(4, 1));
  cgen.add_cluster(point(4, 8));
  cgen.add_cluster(point(10, 7));

  for (size_t r=0; r < num_points; r++) {
    point p;
    do {
      p = cgen.next_point();
    } while (p.x < 0 || p.y < 0);
    points.push_back(p);
  }

  kmedoids km, clara;
  build_dissimilarity_matrix(points, point_distance(), dissimilarity);
  if (debug) {
    km.set_xcallback(print_cluster_info_if_6);
    clara.set_xcallback(print_cluster_info_if_6);
  }

  double xpam_bic   = km.xpam(dissimilarity, max_clusters, dimensions);
  double xclara_bic = clara.xclara(points, point_distance(), max_clusters, dimensions);

  cout << "========== Best XPAM ==========" << endl;
  print_cluster_info(km, xpam_bic);

  cout << "========== Best XCLARA ==========" << endl;
  print_cluster_info(clara, xclara_bic);
}
