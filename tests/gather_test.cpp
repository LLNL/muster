#include <mpi.h>
#include <algorithm>
#include <vector>
#include <iterator>
#include <string>
#include "point.h"
#include "packable_vector.h"
#include "gather.h"
#include "Timer.h"

using namespace cluster;
using namespace std;

static bool verbose = false;
static bool debug = false;

void verify(const vector< packable_vector<point> >& all_points, int rank) {
  if (debug) cout << "Points on rank " << rank << ":" << endl;
  for (size_t i = 0; i < all_points.size(); i++) {
    vector<point>& packables(*all_points[i]._packables);

    for (size_t j = 0; j < packables.size(); j++) {
      if (debug) cout << packables[j] << " ";
      if (packables[j].x != i || packables[j].y != i) {
        cerr << "FAILED: Wrong value from rank " << i << ": " << packables[j] << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
    }
  }
  if (debug) cout << endl;
}


int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  // parse args
  for (int i=0; i < argc; i++) {
    if (string(argv[i]) == string("-v")) {
      verbose = true;
    }
    if (string(argv[i]) == string("-d")) {
      debug = true;
    }
  }

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  // Local vector of points -- we're going to gather a bunch of these.
  vector<point> my_points;
  my_points.push_back(point(rank, rank));
  my_points.push_back(point(rank, rank));

  // Make a packable vector handle out of our local items, and put them all in an
  // array of packable vectors for the default root
  vector< packable_vector<point> > all_points;
  gather(packable_vector<point>(&my_points, false), all_points, MPI_COMM_WORLD);

  // verify that everything in this vector is from the correct rank for all roots
  Timer timer;
  for (int root = 1; root < size; root++) {
    all_points.clear();
    gather(packable_vector<point>(&my_points, false), all_points, MPI_COMM_WORLD, root);
    timer.record("gather");

    if (rank == root) {
      verify(all_points, root);
    }
    timer.record("verify");
  }

  if (rank == 0 && verbose) {
    cout << "Average Gather time: " << timer["gather"] / size / 1e9 << " sec" << endl;
    cout << "PASSED" << endl;
  }


  // verify that everything in this vector is from the correct rank for all roots
  for (int root = 1; root < size; root++) {
    all_points.clear();
    allgather(packable_vector<point>(&my_points, false), all_points, MPI_COMM_WORLD, root);
    timer.record("allgather");

    verify(all_points, root);
    timer.record("verify");
  }

  if (rank == 0 && verbose) {
    cout << "Average AllGather time: " << timer["allgather"] / size / 1e9 << " sec" << endl;
    cout << "PASSED" << endl;
  }


  MPI_Finalize();
  return 0;
}
