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
  cerr << "Usage: par-bic-test [-htvi] [-p points-per-process] [-c clusters]" << endl;
  cerr << "  Parallel test case for clustering points and evaluating their BIC scores." << endl;
  cerr << "Options:" << endl;
  cerr << "  -h         Show this message." << endl;
  cerr << "  -t         Output timing info to file." << endl;
  cerr << "  -v         Validate with sequential clustering and output Mirkin distance." << endl;
  cerr << "  -i         Number of times to iterate through test." << endl;
  cerr << "               Default is 10 tries." << endl;
  cerr << "  -p         Number of points generated per process." << endl;
  cerr << "               Default is 1." << endl;
  cerr << "  -c         Max number of clusters to search for." << endl;
  cerr << "               Default is 1." << endl;
  exit(1);
}

bool timing = false;
bool validate = false;
size_t points_per_process = 1;
size_t max_clusters = 10;
size_t iterations = 10;
const size_t dimensions = 2;

const double xscale = 4.0;
const double yscale = 1.0;

/// Uses getopt to read in arguments.
void get_args(int *argc, char ***argv, int rank) {
  int c;
  char *err;

  while ((c = getopt(*argc, *argv, "htvi:l:n:s:c:")) != -1) {
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
    case 'i':
      iterations = strtol(optarg, &err, 0);
      if (*err) usage();
      break;
    case 'p':
      points_per_process = strtol(optarg, &err, 0);
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



int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  get_args(&argc, &argv, rank);

  // make sure max_clusters is valid.
  max_clusters = min(max_clusters, size*points_per_process);
  
  uint32_t seed = 0;
  struct timeval time;
  gettimeofday(&time, 0);
  seed = time.tv_sec * time.tv_usec;

  MTRand rand(seed + rank);  // Seed points differently on each rank.
  matrix<int> data(points_per_process, dimensions);     // vector of local "traces"

  for (size_t r=0; r < data.size1(); r++) {
    point p;

    size_t type = (rank + r) % 3;
    switch (type) {
    case 0:
      p = create_point(rand, 2, 4, 2, 2, xscale, yscale);
      break;
    case 1:
      p = create_point(rand, 12, 12, 2, 2, xscale, yscale);
      break;
    case 2:
      p = create_point(rand, 12, 3, 2, 2, xscale, yscale);
      break;
    }
    data(r,0) = p.x;
    data(r,1) = p.y;
  }

  vector<point> points(data.size1());
  for (size_t i=0; i < data.size1(); i++) {
    points[i] = point(data(i,0), data(i,1));
  }

  double total_mirkin = 0;
  for (size_t iter=0; iter < iterations; iter++) {
    par_kmedoids parkm;
    
    const size_t trials = 10;
    long long start = get_time_ns();
    
    for (size_t i=0; i < trials; i++) {
      parkm.xclara(points, point_distance(), max_clusters, dimensions);
    }

    double total = get_time_ns() - start;
    double avg = total / trials;
    if (rank == 0) {
      cout << size << " processes" << endl;
      cout << "TOTAL:   " << total / 1e9 << endl;
      cout << "AVERAGE: " << avg   / 1e9 << endl;
  
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


      matrix<int> full_data(size * points_per_process, dimensions);
      MPI_Gather(&data(0,0),      points_per_process * dimensions, MPI_INT,
                 &full_data(0,0), points_per_process * dimensions, MPI_INT,
                 0, MPI_COMM_WORLD);
      
      if (rank == 0) {
        //ostringstream fn;
        //fn << "full." << iter;
        //ofstream full(fn.str().c_str());
        //output(full_data, full);

        cerr << "pushing" << endl;

        // compare parallel clustering with local clustering.
        vector<point> all_points;
        for (size_t i=0; i < full_data.size1(); i++) {
          all_points.push_back(point(full_data(i,0), full_data(i,1)));
        }

        cerr << "done pushing" << endl;
        
        dissimilarity_matrix distance;
        build_dissimilarity_matrix(all_points, point_distance(), distance);
        

        cerr << "done dissim" << endl;

        kmedoids km;
        double best_bic = km.xpam(distance, max_clusters, dimensions);
        cout << endl;
        cout << "Seq k:   " << km.num_clusters() << endl;
        cout << "Seq BIC: " << best_bic << endl;
        cout << "Seq D:   " << total_dissimilarity(km, matrix_distance(distance)) << endl;
        cout << "Seq D2:  " << total_squared_dissimilarity(km, matrix_distance(distance)) << endl;

        cerr << "done km" << endl;

        cout << km << endl;
        cerr << "done out km" << endl;
        cerr << "size: " << all_points.size() << endl;

        copy(all_points.begin(), all_points.end(), ostream_iterator<point>(cerr, " "));


        draw("SEQ", all_points, km);

        cerr << "done draw" << endl;

        cout << endl;
        cout << "Par k:   " << parallel.num_clusters() << endl;
        cout << "Par BIC: " << bic(parallel, matrix_distance(distance), dimensions) << endl;
        cout << "Par D:   " << total_dissimilarity(parallel, matrix_distance(distance)) << endl;
        cout << "Par D2:  " << total_squared_dissimilarity(parallel, matrix_distance(distance)) << endl;

        //cout << parallel << endl;
        draw("PAR", all_points, parallel);

        cout << endl;

        kmedoids par_clone;
        par_clone.pam(distance, parallel.num_clusters());
        cout << "Clone k:   " << par_clone.num_clusters() << endl;
        cout << "Clone BIC: " << bic(par_clone, matrix_distance(distance), dimensions) << endl;
        cout << "Clone D:   " << total_dissimilarity(par_clone, matrix_distance(distance)) << endl;
        cout << "Clone D2:  " << total_squared_dissimilarity(par_clone, matrix_distance(distance)) << endl;
        //cout << par_clone << endl;
        draw("CLONE", all_points, par_clone);
        cout << endl;

        
        double mirkin = mirkin_distance(parallel, km);
        total_mirkin += mirkin;
        cout << "Mirkin Distance:       " << mirkin << endl;

        double clone_dist = mirkin_distance(par_clone, parallel);
        cout << "Clone Mirkin Distance: " << clone_dist << endl;

      }
    }
  }

  if (validate && rank == 0) {
    cout << endl;
    cout << "Average Mirkin Distance: " << total_mirkin / iterations << endl;
  }
}
