#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <sys/time.h>

#include "timing.h"
#include "point.h"
#include "bic.h"
#include "par_kmedoids.h"
#include "MersenneTwister.h"

using namespace cluster;
using namespace std;


void usage() {
  cerr << "Usage: par-cluster-speed-test [-htvi] [-n num-points] [-c clusters] [-s scale]" << endl;
  cerr << "  Compare parallel clustering with sequential clustering." << endl;
  cerr << "Options:" << endl;
  cerr << "  -h         Show this message." << endl;
  cerr << "  -x         Use BIC-scored versions of PAM and CLARA." << endl;
  cerr << "  -t         Save details timing info in a file." << endl;
  cerr << "  -n         Number of points per process." << endl;
  cerr << "               Default is 1." << endl;
  cerr << "  -i         Initial sample size in clara (before 2*k is added)." << endl;
  cerr << "               Default is 40." << endl;
  cerr << "  -r         Number of trials per k in clara." << endl;
  cerr << "               Default is 5." << endl;
  cerr << "  -k         Max number of clusters to search for." << endl;
  cerr << "               Default is 10." << endl;
  exit(1);
}

size_t objects_per_process = 1;
size_t num_clusters = 10;
size_t init_size = 40;
size_t max_reps = 5;
bool use_bic = false;
bool timing = false;

/// Uses getopt to read in arguments.
void get_args(int *argc, char ***argv, int rank) {
  int c;
  char *err;

  while ((c = getopt(*argc, *argv, "htxn:i:r:k:")) != -1) {
    switch (c) {
    case 'h':
      if (rank == 0) usage();
      exit(1);
      break;
    case 'x':
      use_bic = true;
      break;
    case 't':
      timing = true;
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




int main(int argc, char **argv) { 
  MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  get_args(&argc, &argv, rank);

  uint32_t seed = 0;
  if (rank == 0) {
    struct timeval time;
    gettimeofday(&time, 0);
    seed = time.tv_sec * time.tv_usec;
  }

  MTRand rand(seed + rank);  // generator to make points.  Seed differently on each rank.
  vector<point> points;      // vector of local points  
  
  // generate randomly distributed, zero-centered points.
  for (size_t i=0; i < objects_per_process; i++) {
    int x = rand.randInt(5000) - 2500;
    int y = rand.randInt(5000) - 2500;
    points.push_back(point(x,y));
  }

  par_kmedoids parkm;
  parkm.set_init_size(init_size);
  parkm.set_max_reps(max_reps);

  // trials of whole algorithm, to account for any variability
  const size_t trials = 10;

  long long start = get_time_ns();
  for (size_t i=0; i < trials; i++) {
    parkm.clara(points, point_distance(), num_clusters);
  }

  double total = get_time_ns() - start;
  double avg = total / trials;
  
  if (rank == 0) {
    if (timing) {
      ostringstream timing_filename;
      timing_filename << "times-" << size;
      ofstream timing_file(timing_filename.str().c_str());
      parkm.get_timer().write(timing_file);
    }

    cout << "PROCS:   " << size << " " << endl;
    cout << "TOTAL:   " << total / 1e9 << endl;
    cout << "AVERAGE: " << avg   / 1e9 << endl;
  }
}
