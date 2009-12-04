#include "par_partition.h"

#include <iostream>

#include "mpi_utils.h"

#define DEBUG

using namespace std;

namespace cluster {
  
  par_partition::par_partition(MPI_Comm _comm) : my_id(0), comm(_comm) { 
    medoids.push_back(0);
  }


  void par_partition::gather(partition& destination, int root) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
#ifdef DEBUG
    size_t count = medoids.size();
    MPI_Bcast(&count, 1, MPI_SIZE_T, root, comm);
    if (count != medoids.size()) {
      cerr << "Error: incorrect number of medoids on " << rank << endl;
      exit(1);
    }
    
    std::vector<object_id> bcast_medoids = medoids;
    MPI_Bcast(&bcast_medoids[0], bcast_medoids.size(), MPI_SIZE_T, root, comm);
    for (size_t i=0; i < medoids.size(); i++) {
      if (medoids[i] != bcast_medoids[i]) {
        cerr << "Error: medoids do not match on " << rank << endl;
        exit(1);
      }
    }
#endif // DEBUG
    
    if (rank == root) {
      destination.cluster_ids.resize(size);
    }

    MPI_Gather(&my_id, 1, MPI_SIZE_T, 
               &destination.cluster_ids[0], 1, MPI_SIZE_T, 
               root, comm);
  }

} // namespace cluster
