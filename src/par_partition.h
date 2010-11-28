//////////////////////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2010, Lawrence Livermore National Security, LLC.  
// Produced at the Lawrence Livermore National Laboratory  
// LLNL-CODE-433662
// All rights reserved.  
//
// This file is part of Muster. For details, see http://github.com/tgamblin/muster. 
// Please also read the LICENSE file for further information.
//
// Redistribution and use in source and binary forms, with or without modification, are
// permitted provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright notice, this list of
//    conditions and the disclaimer below.
//  * Redistributions in binary form must reproduce the above copyright notice, this list of
//    conditions and the disclaimer (as noted below) in the documentation and/or other materials
//    provided with the distribution.
//  * Neither the name of the LLNS/LLNL nor the names of its contributors may be used to endorse
//    or promote products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
// OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
// LAWRENCE LIVERMORE NATIONAL SECURITY, LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
// ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//////////////////////////////////////////////////////////////////////////////////////////////////

///
/// @file par_partition.h
/// @author Todd Gamblin tgamblin@llnl.gov
/// @brief Distributed representation of a partitioning of a data set.
///
#ifndef PAR_PARTITION_H
#define PAR_PARTITION_H

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


