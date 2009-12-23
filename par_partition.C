#include "par_partition.h"

#include <iostream>

#include "mpi_utils.h"

#define DEBUG

using namespace std;

namespace cluster {
  
  par_partition::par_partition(MPI_Comm _comm) : comm(_comm) { }

  par_partition::~par_partition() { }


  void par_partition::gather(partition& destination, int root) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
#ifdef DEBUG
    size_t count = medoid_ids.size();
    MPI_Bcast(&count, 1, MPI_SIZE_T, root, comm);
    if (count != medoid_ids.size()) {
      cerr << "Error: incorrect number of medoids on " << rank << ": " << medoid_ids.size() 
           << ", expected " << count << endl;
      exit(1);
    }
    
    size_t object_count = cluster_ids.size();
    MPI_Bcast(&object_count, 1, MPI_SIZE_T, root, comm);    
    if (object_count != cluster_ids.size()) {
      cerr << "Error: incorrect number of objects on " << rank << ": " << cluster_ids.size() 
           << ", expected " << object_count << endl;
      exit(1);
    }

    std::vector<object_id> bcast_medoid_ids = medoid_ids;
    MPI_Bcast(&bcast_medoid_ids[0], bcast_medoid_ids.size(), MPI_SIZE_T, root, comm);
    for (size_t i=0; i < medoid_ids.size(); i++) {
      if (medoid_ids[i] != bcast_medoid_ids[i]) {
        cerr << "Error: medoids do not match on " << rank << endl;
        exit(1);
      }
    }
#endif // DEBUG
    
    if (rank == root) {
      destination.medoid_ids = medoid_ids;
      destination.cluster_ids.resize(cluster_ids.size() * size);
    }

    MPI_Gather(&cluster_ids[0],             cluster_ids.size(), MPI_SIZE_T,
               &destination.cluster_ids[0], cluster_ids.size(), MPI_SIZE_T, 
               root, comm);
  }

} // namespace cluster
