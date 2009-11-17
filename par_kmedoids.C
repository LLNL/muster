#include "par_kmedoids.h"

#include <cstdlib>
#include <sys/time.h>
using namespace std;

namespace cluster {

  void par_kmedoids::seed_random(MPI_Comm comm) {
    int rank;
    MPI_Comm_rank(comm, &rank);

    int seed = 0;
    if (rank == 0) {
      struct timeval time;
      gettimeofday(&time, 0);
    
      seed = time.tv_sec * time.tv_usec;
    }

    MPI_Bcast(&seed, 1, MPI_INT, 0, comm);
    random.seed(seed);
  }

} // namespace cluster
