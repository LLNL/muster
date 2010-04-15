#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <sys/time.h>

#include "point.h"
#include "kmedoids.h"
#include "par_kmedoids.h"
#include "bic.h"
#include "matrix_utils.h"
#include "spherical_clustering_generator.h"

using namespace cluster;
using namespace std;
using boost::numeric::ublas::matrix;

void usage() {
  cerr << "Usage: par-bic-test [-htvdx] [-p points-per-process] [-c clusters] [-s scale] [-t trials]" << endl;
  cerr << "  Parallel test case for clustering points and evaluating their BIC scores." << endl;
  cerr << "Options:" << endl;
  cerr << "  -h         Show this message." << endl;
  cerr << "  -t         Output timing info to file." << endl;
  cerr << "  -v         Validate with sequential clustering and output Mirkin distance." << endl;
  cerr << "  -d         Verbose debug output, actually draw clusterings." << endl;
  cerr << "  -s         Scale point values by a factor." << endl;
  cerr << "               Default is 1.0" << endl;
  cerr << "  -p         Number of points generated per process." << endl;
  cerr << "               Default is 1." << endl;
  cerr << "  -k         Max number of clusters to search for (or number to search for if no BIC)." << endl;
  cerr << "               Default is 1." << endl;
  cerr << "  -t         Number of trials of parallel algorithms to run for speed measurement." << endl;
  cerr << "               Default is 10." << endl;
  cerr << "  -i         Number of times to iterate through complete test with validation." << endl;
  cerr << "               Default is 10." << endl;
  exit(1);
}

bool timing = false;
bool validate = false;
bool debug = false;
size_t points_per_process = 1;
size_t num_clusters = 10;
size_t speed_trials = 10;
size_t iterations = 10;
const size_t dimensions = 2;

double scale = 1.0;
const double xscale = 2.5;
const double yscale = 1.0;

/// Uses getopt to read in arguments.
void get_args(int *argc, char ***argv, int rank) {
  int c;
  char *err;

  while ((c = getopt(*argc, *argv, "htvds:i:p:k:")) != -1) {
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
    case 'd':
      debug = true;
      break;
    case 's':
      scale = strtod(optarg, &err);
      if (*err) usage();
      break;
    case 'i':
      iterations = strtol(optarg, &err, 0);
      if (*err) usage();
      break;
    case 'p':
      points_per_process = strtol(optarg, &err, 0);
      if (*err) usage();
      break;
    case 'k':
      num_clusters = strtol(optarg, &err, 0);
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

/// Globals for sequential clustering info, needed by print_cluster_info
vector<point> all_points;
dissimilarity_matrix dissimilarity;

/// Callback for printing out clusterings and their BIC scores.
static void print_cluster_info(const cluster::partition& km, double bic, string name="") {
  cout << name << " ";
  for (size_t i=0; i < name.size(); i++) cout << "=";
  cout << endl;

  cout << "  k:           " << km.num_clusters() << endl;
  cout << "  BIC:         " << bic << endl;
  cout << "  actual  bic: " << cluster::bic(km, matrix_distance(dissimilarity), dimensions) << endl;
  cout << "  old bic:     " << old_bic(km, matrix_distance(dissimilarity), dimensions) << endl;  
    
  if (debug) {
    cout << "medoids: ";
    for (size_t i=0; i < km.medoid_ids.size(); i++) {
      cout << all_points[km.medoid_ids[i]] << " ";
    }
    cout << endl;
  }

  if (debug) {
    cout << "D:       " << total_dissimilarity(km, matrix_distance(dissimilarity)) << endl;
    cout << "D2:      " << total_squared_dissimilarity(km, matrix_distance(dissimilarity)) << endl;
    cout << km << endl;
  }

  draw("Clustering", all_points, km);
  cout << endl;
}


int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // get args and make sure num_clusters is valid.
  get_args(&argc, &argv, rank);
  size_t num_points = size * points_per_process;
  num_clusters = min(num_clusters, num_points);

  // create generators for a bunch of clusters.
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

  // do the generating sequentially, since boost isn't that great for
  // parallel random number generation.
  matrix<double> full_data(num_points, dimensions);
  for (size_t r=0; r < num_points; r++) {
    point p;
    do {
      p = cgen.next_point();
    } while (p.x < 0 || p.y < 0);

    // copy all_points into matrix of raw point data, for easy (kludgy) MPI transfer.
    all_points.push_back(p);
    full_data(r,0) = p.x;
    full_data(r,1) = p.y;
  }

  // scatter points to processes
  matrix<double> data(points_per_process, dimensions);
  MPI_Scatter(&full_data(0,0),  data.size1() * data.size2(), MPI_DOUBLE,
              &data(0,0),       data.size1() * data.size2(), MPI_DOUBLE,
              0, MPI_COMM_WORLD);

  // stick received points in local point vector.
  vector<point> points(points_per_process);
  for (size_t i=0; i < points.size(); i++) {
    points[i] = point(data(i,0), data(i,1));
  }


  double total_mirkin = 0;
  for (size_t iter=0; iter < iterations; iter++) {
    par_kmedoids parkm;

    //
    // First we do speed trials of the parallel clustering algorithm.
    //
    long long start = get_time_ns();    
    double best_bic = 0;
    for (size_t i=0; i < speed_trials; i++) {
      best_bic = parkm.xclara(points, point_distance(), num_clusters, dimensions);
      //parkm.clara(points, point_distance(), num_clusters);
    }
    double total = get_time_ns() - start;
    double avg = total / speed_trials;

    
    if (rank == 0) {
      cout << size << " processes" << endl;
      cout << "TOTAL:    " << total  / 1e9 << endl;
      cout << "AVERAGE:  " << avg / 1e9 << endl;
  
      if (timing) {
        ostringstream timing_filename;
        timing_filename << "times-" << size;
        ofstream timing(timing_filename.str().c_str());
        parkm.get_timer().write(timing);
      }
    }
    

    if (validate) {
      cluster::partition parallel;
      parkm.gather(parallel, 0);


      if (rank == 0) {
        cerr << points.size() << " local points." << endl;
        cerr << all_points.size() << " points." << endl;

        // do sequential xpam clustering to compare.
        kmedoids km;
        build_dissimilarity_matrix(all_points, point_distance(), dissimilarity);
        double best_bic = km.xpam(dissimilarity, num_clusters, dimensions);

        print_cluster_info(km, best_bic, "XPAM");
        print_cluster_info(parallel, 0, "PAR XCLARA");

        cout << "k: " << num_clusters << endl;

        double mirkin = mirkin_distance(km, parallel);
        total_mirkin += mirkin;

        cout << "  mirkin(xpam, par_xclara):    " << mirkin << endl;
      }
    }
  }

  if (validate && rank == 0) {
    cout << endl;
    cout << "Average Mirkin Distance: " << total_mirkin / iterations << endl;
  }
}
