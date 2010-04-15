#include "par_partition.h"

#include <iostream>

#include "mpi_utils.h"
#include "partition.h"

//#define DEBUG

using namespace std;

namespace cluster {
  
  par_partition::par_partition(MPI_Comm _comm) : comm(_comm) { }

  par_partition::~par_partition() { }


  void par_partition::gather(partition& destination, int root) {
    int rank, size;
    PMPI_Comm_rank(comm, &rank);
    PMPI_Comm_size(comm, &size);
    
#ifdef DEBUG
    size_t count = medoid_ids.size();
    PMPI_Bcast(&count, 1, MPI_SIZE_T, root, comm);
    if (count != medoid_ids.size()) {
      cerr << "Error: incorrect number of medoids on " << rank << ": " << medoid_ids.size() 
           << ", expected " << count << endl;
      exit(1);
    }
    
    size_t object_count = cluster_ids.size();
    PMPI_Bcast(&object_count, 1, MPI_SIZE_T, root, comm);    
    if (object_count != cluster_ids.size()) {
      cerr << "Error: incorrect number of objects on " << rank << ": " << cluster_ids.size() 
           << ", expected " << object_count << endl;
      exit(1);
    }

    std::vector<object_id> bcast_medoid_ids = medoid_ids;
    PMPI_Bcast(&bcast_medoid_ids[0], bcast_medoid_ids.size(), MPI_SIZE_T, root, comm);
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

    PMPI_Gather(&cluster_ids[0],             cluster_ids.size(), MPI_SIZE_T,
                &destination.cluster_ids[0], cluster_ids.size(), MPI_SIZE_T, 
                root, comm);
  }


  void par_partition::get_sizes(std::vector<size_t>& sizes) {
    vector<size_t> local_sizes(medoid_ids.size(), 0);
    for (size_t i=0; i < cluster_ids.size(); i++) {
      local_sizes[cluster_ids[i]]++;
    }
    
    sizes.resize(medoid_ids.size());
    PMPI_Allreduce(&local_sizes[0], &sizes[0], medoid_ids.size(), MPI_SIZE_T, MPI_SUM, comm); 
  }



  std::ostream& operator<<(std::ostream& out, const par_partition& par) {
    // TODO: fix this; it's hacky.
    cluster::partition p;
    p.medoid_ids = par.medoid_ids;
    p.cluster_ids = par.cluster_ids;
    out << p;
    return out;
  }


} // namespace cluster
