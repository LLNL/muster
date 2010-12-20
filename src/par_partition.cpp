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
/// @file par_partition.cpp
/// @author Todd Gamblin tgamblin@llnl.gov
///
#include "par_partition.h"

#include <iostream>

#include "mpi_utils.h"
#include "partition.h"
#include "mpi_bindings.h"

//#define DEBUG

using namespace std;

namespace cluster {
  
  par_partition::par_partition(MPI_Comm _comm) : comm(_comm) { }

  par_partition::~par_partition() { }


  void par_partition::gather(partition& destination, int root) {
    int rank, size;
    CMPI_Comm_rank(comm, &rank);
    CMPI_Comm_size(comm, &size);
    
#ifdef DEBUG
    size_t count = medoid_ids.size();
    CMPI_Bcast(&count, 1, MPI_SIZE_T, root, comm);
    if (count != medoid_ids.size()) {
      cerr << "Error: incorrect number of medoids on " << rank << ": " << medoid_ids.size() 
           << ", expected " << count << endl;
      exit(1);
    }
    
    size_t object_count = cluster_ids.size();
    CMPI_Bcast(&object_count, 1, MPI_SIZE_T, root, comm);    
    if (object_count != cluster_ids.size()) {
      cerr << "Error: incorrect number of objects on " << rank << ": " << cluster_ids.size() 
           << ", expected " << object_count << endl;
      exit(1);
    }

    std::vector<object_id> bcast_medoid_ids = medoid_ids;
    CMPI_Bcast(&bcast_medoid_ids[0], bcast_medoid_ids.size(), MPI_SIZE_T, root, comm);
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

    CMPI_Gather(&cluster_ids[0],             cluster_ids.size(), MPI_SIZE_T,
                &destination.cluster_ids[0], cluster_ids.size(), MPI_SIZE_T, 
                root, comm);
  }


  void par_partition::get_sizes(std::vector<size_t>& sizes) {
    vector<size_t> local_sizes(medoid_ids.size(), 0);
    for (size_t i=0; i < cluster_ids.size(); i++) {
      local_sizes[cluster_ids[i]]++;
    }
    
    sizes.resize(medoid_ids.size());
    CMPI_Allreduce(&local_sizes[0], &sizes[0], medoid_ids.size(), MPI_SIZE_T, MPI_SUM, comm); 
  }



  std::ostream& operator<<(std::ostream& out, const par_partition& par) {
    cluster::partition p;
    p.medoid_ids = par.medoid_ids;
    p.cluster_ids = par.cluster_ids;
    out << p;
    return out;
  }


} // namespace cluster
