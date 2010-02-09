#include <mpi.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <iterator>
#include <iomanip>
#include <sstream>

#include "multi_gather.h"
#include "point.h"
#include "random.h"
#include "Timer.h"
#include "MersenneTwister.h"

using namespace std;
using namespace cluster;


template <class OutputIterator>
void generate_points_for_rank(int rank, OutputIterator out) {
  for (int i=0; i < (rank % 10 + 1); i++) {
    int p = rank * 10 + i;
    *out++ = point(p,p);
  }
}



int main(int argc, char **argv) {
  Timer timer;
  MPI_Init(&argc, &argv);
  
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  bool verbose = false;
  bool timing  = false;

  for (int i=1; i < argc; i++) {
    if (string(argv[i]) == "-v") {
      verbose = true;
    }

    if (string(argv[i]) == "-t") {
      timing = true;
    }
  }

  MTRand random;
  random.seed(1);  // make sure all ranks generate the same numbers.

  vector<point> points;   // local points to send
  vector<point> dest;     // destination vector for gathered points
  vector<int> sources;    // ranks we received from, so we can check the local points vector

  generate_points_for_rank(rank, back_inserter(points));

  timer.record("init");

  // now fire off <size> gathers, each with ~size/2 elements.
  multi_gather<point> gather(MPI_COMM_WORLD);
  for (int root=0; root < size; root++) {
    vector<int> cur_sources;
    random_subset(size, (int)ceil(sqrt((double)size)), back_inserter(cur_sources), random);
    gather.start(points.begin(), points.end(), cur_sources.begin(), cur_sources.end(), dest, root);

    if (rank == root) {
      // record sources so we can check later.
      cur_sources.swap(sources);
    }
  }

  timer.record("start_gathers");
  
  if (verbose) {
    cerr << rank << " gathering from [";
    for (size_t i=0; i < sources.size(); i++) cerr << setw(3) << sources[i] << " ";
    cerr << "]" << endl;
    timer.record("verbose");
  }

  if (verbose) cerr << rank << " calling multi_gather::finish()." << endl;

  gather.finish();
  timer.record("finish_gathers");

  if (verbose) cerr << rank << " finished gather." << endl;
  
  int passed = 1;
  size_t index = 0;
  for (size_t i=0; passed && i < sources.size(); i++) {
    vector<point> expected;
    generate_points_for_rank(sources[i], back_inserter(expected));
    
    for (size_t j=0; j < expected.size(); j++) {
      if (index >= dest.size() || dest[index] != expected[j]) {
        passed = 0;
      }
      index++;
    }
  }
  timer.record("check");

  if (verbose) {
    ostringstream msg;

    msg << rank << " Expected: ";
    for (size_t i=0; i < sources.size(); i++) {
      generate_points_for_rank(sources[i], ostream_iterator<point>(msg, " "));
    }
    msg << endl;

    msg << rank << " Found:   ";
    for (size_t i=0; i < dest.size(); i++) {
      msg << " " << dest[i];
    }
    msg << endl;
    cerr << msg.str();
    timer.record("verbose");
  }

  int num_passed = 0;
  MPI_Allreduce(&passed, &num_passed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Finalize();

  bool all_passed = (num_passed == size);
  if (rank == 0) {
    cerr << (all_passed ? "PASSED" : "FAILED") << endl;
  }

  if (timing && rank == 0) timer.dump(cerr);

  return all_passed ? 0 : 1;
}
