#include <iostream>

#include "par_kmedoids.h"
#include "point.h"

using namespace cluster;
using namespace std;


int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  
  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  

  // local data
  point p(rank, rank);

  cout << p << endl;

  // parallel clusterer for data
  par_kmedoids cluster;
  vector< source_pair<point> > medoids;
  cluster.par_clara(p, point_distance(), 10, medoids);

  MPI_Finalize();
}
