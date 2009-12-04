#ifndef PAR_PARTITION_H
#define PAR_PARTITION_H

#include <mpi.h>
#include <vector>

#include "partition.h"


namespace cluster {
  
  /// Parallel partition object.  
  struct par_partition {
    /// Gives the rank of the process owning the ith medoid.
    /// medoids[i] == rank of process owning ith medoid
    std::vector<object_id> medoids;
    
    /// Identity of medoid this process is a member of.
    medoid_id my_id;
    
    /// Communicator, the processes of which this partition divides
    MPI_Comm comm;

    /// Construct a parallel partition for the communicator supplied
    /// Partition starts off with everyone in one cluster with medoid 0.
    par_partition(MPI_Comm comm = MPI_COMM_WORLD);

    /// Collective operation.  Gathers my_id from all processes into a 
    /// local partition object. If size of system is large, then this method
    /// will not scale.
    void gather(partition& local, int root=0);
  };
  
} // namespace cluster

#endif // PAR_PARTITION_H


