#ifndef MUSTER_GATHER_H
#define MUSTER_GATHER_H

#include <mpi.h>
#include <numeric>
#include <algorithm>
#include "mpi_bindings.h"
#include "mpi_utils.h"
#include "binomial.h"

namespace cluster {

  ///
  /// Binomial gather of char buffers into a single agglomerated clump of buffers
  ///
  template <class T>
  void gather(const T& src, std::vector<T>& dest, MPI_Comm comm, int root = 0) {
    int rank, size;
    CMPI_Comm_rank(comm, &rank);
    CMPI_Comm_size(comm, &size);

    binomial_embedding binomial(size, root);
    int parent = binomial.parent(rank);
    std::vector<int> children = binomial.children(rank);

    std::vector<int> sizes;                  // sizes of buffers to receive.
    sizes.push_back(src.packed_size(comm));  // size of local packed data

    for (size_t i=0; i < children.size(); i++) {
      // Receive sizes from all children
      int cur_size;
      CMPI_Recv(&cur_size, 1, MPI_INT, children[i], 0, comm, MPI_STATUS_IGNORE);
      sizes.push_back(cur_size);
    }

    // construct offsets to receive buffers into.
    std::vector<int> offsets;
    offsets.push_back(0);
    std::partial_sum(sizes.begin(), sizes.end(), back_inserter(offsets));

    // create a buffer to receive into, then to send to parent
    std::vector<char> sendbuf(accumulate(sizes.begin(), sizes.end(), 0));
    
    // pack local object before its children
    int pos = 0;
    src.pack(&sendbuf[0], sendbuf.size(), &pos, comm);

    for (size_t i = 0; i < children.size(); i++) {
      CMPI_Recv(&sendbuf[offsets[i+1]], sizes[i+1], MPI_PACKED, children[i], 0, comm, MPI_STATUS_IGNORE);
    }
    
    if (parent != -1) {
      // Done receiving.  Send to parent.
      int size = sendbuf.size();
      CMPI_Send(&size,                    1,    MPI_INT, parent, 0, comm);
      CMPI_Send(&sendbuf[0], sendbuf.size(), MPI_PACKED, parent, 0, comm);

    } else {
      // unpack everything at the root, which has no parent.
      int pos = 0;
      dest.resize(size);
      for (size_t i=0; i < size; i++) {
        dest[binomial.reverse_relative_rank(i)] = T::unpack(&sendbuf[0], sendbuf.size(), &pos, comm);
      }
    }
  }

} // namespace cluster

#endif // MUSTER_GATHER_H
