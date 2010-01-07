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


/// Test using points instead of QGrams to make sure clustering 
/// algorithms work.
int main(int argc, char **argv) { 
  MPI_Init(&argc, &argv);

  size_t objects_per_process = 1;
  if (argc > 1) {
    objects_per_process = strtol(argv[1], NULL, 0);    
  }
  
  size_t max_clusters = 10;
  if (argc > 2) {
    max_clusters = strtol(argv[2], NULL, 0);
  }

  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

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
  const size_t trials = 10;

  long long start = get_time_ns();
  for (size_t i=0; i < trials; i++) {
    parkm.xclara(points, point_distance(), max_clusters, 2);
  }

  double total = get_time_ns() - start;
  double avg = total / trials;
  
  if (rank == 0) {
    ostringstream timing_filename;
    timing_filename << "times-" << size;
    ofstream timing(timing_filename.str().c_str());
    parkm.get_timer().write(timing);

    cout << size << " processes" << endl;
    cout << "TOTAL:   " << total / 1e9 << endl;
    cout << "AVERAGE: " << avg   / 1e9 << endl;
  }
}
