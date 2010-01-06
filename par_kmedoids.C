#include "par_kmedoids.h"

#include <cstdlib>
#include <stdint.h>
#include <sys/time.h>
using namespace std;

namespace cluster {



  par_kmedoids::par_kmedoids(MPI_Comm comm) 
    : par_partition(comm), 
      total_dissimilarity(numeric_limits<double>::infinity()),
      min_bic_score(0),
      init_size(40),
      max_reps(5),
      epsilon(1e-15)
  { }


  void par_kmedoids::set_epsilon(double e) {
    epsilon = e;
  }


  double par_kmedoids::average_dissimilarity() {
    int size;
    PMPI_Comm_size(comm, &size);
    return total_dissimilarity / size;
  }

  double par_kmedoids::bic_score() {
    return min_bic_score;
  }

  void par_kmedoids::seed_random_uniform(MPI_Comm comm) {
    int rank;
    PMPI_Comm_rank(comm, &rank);

    uint32_t seed = 0;
    if (rank == 0) {
      struct timeval time;
      gettimeofday(&time, 0);
    
      seed = time.tv_sec * time.tv_usec;
    }

    PMPI_Bcast(&seed, 1, MPI_INT, 0, comm);
    random.seed(seed);
  }

} // namespace cluster
