#ifndef PAR_PARTITION_H
#define PAR_PARTITION_H
///
/// @file par_partition.h
/// @brief Distributed representation of a partitioning of a data set.
///

#include <mpi.h>
#include <vector>
#include <ostream>

#include "partition.h"


namespace cluster {
  
  ///
  /// par_partition represents a partitioning of a distributed data set.
  /// It is analogous to partition, but its object_ids are distributed across
  /// the ranks of the communicator it is instantiated with.  Each process is assumed
  /// to "own" some set of objects in the data set this describes, and each process's
  /// par_partition object contains medoid_ids only for its own objects. Thus, the cluster_ids
  /// array will contain different object_ids on different processes within the "same"
  /// par_partition object.
  /// 
  /// While cluster_ids will vary, the "same" par_partition object on the same 
  /// communicator will have the same medoid_ids.  If this is not true, then the medoid_ids
  /// won't make any sense between processes.  Partitioning algorithms that use a 
  /// par_partition for output should preserve this property.
  /// 
  /// You can convert a par_partition to a partition on a single process using 
  /// the gather() method.  This is a collective operation.  It is not scalable, in that 
  /// it will aggregate ids from <i>every</i> process in the communicator to <i>one</i> process.
  /// However, it's useful for small systems and debugging.
  /// 
  /// @see partition, the non-distributed equivalent of this class.
  ///
  struct par_partition {
    /// Gives the object id for the ith medoid.  This object may not be local.
    std::vector<object_id> medoid_ids;
    
    /// Global cluster ids for local objects.  These are indices into medoid_ids.
    /// The object id of the medoid of the ith local object is medoid_ids[cluster_ids[i]].
    std::vector<object_id> cluster_ids;
    
    /// Communicator, the processes of which this partition divides
    MPI_Comm comm;

    /// Construct a parallel partition for the communicator supplied
    /// Partition starts off with everyone in one cluster with medoid 0.
    par_partition(MPI_Comm comm = MPI_COMM_WORLD);

    /// Virtual destructor for inheritance.
    virtual ~par_partition();

    /// Scalably get the sizes of all the clusters in this partition.
    /// POST: sizes is valid on all processes
    void get_sizes(std::vector<size_t>& sizes);

    /// Collective operation.  Gathers my_id from all processes into a 
    /// local partition object. If size of system is large, then this method
    /// will not scale.
    void gather(partition& local, int root=0);
  };

  ///
  /// Right now this just uses parition::operator<<() by making a 
  /// partition with this par_partition's cluster_ids and medoid_ids vectors
  /// and outputting it.
  /// 
  /// @todo fix the implementation; it's a little hacky.  
  ///
  std::ostream& operator<<(std::ostream& out, const par_partition& par);

} // namespace cluster

#endif // PAR_PARTITION_H


