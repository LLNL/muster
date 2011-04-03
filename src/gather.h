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
  /// Packs and gathers a buffer full of packed representation of src's.  All packed
  /// src's are stored in the destination buffer in binomial traversal order on completion.
  /// That is, element i in dest is the value for rank binomial.reverse_relative_rank(i).
  ///
  /// This allows you to use native MPI operations like bcast on the packed buffer 
  /// once it's gathered.
  ///
  /// @see gather() for a version of this that will unpack the gathered data for you.
  ///
  template <class T>
  void gather_packed(const T& src, std::vector<char>& dest, const binomial_embedding binomial, MPI_Comm comm) {
    int rank, size;
    CMPI_Comm_rank(comm, &rank);
    CMPI_Comm_size(comm, &size);

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

    // receive from all children
    for (size_t i = 0; i < children.size(); i++) {
      CMPI_Recv(&sendbuf[offsets[i+1]], sizes[i+1], MPI_PACKED, children[i], 0, comm, MPI_STATUS_IGNORE);
    }
    
    if (parent != -1) {
      // Done receiving.  Send to parent.
      int size = sendbuf.size();
      CMPI_Send(&size,                    1,    MPI_INT, parent, 0, comm);
      CMPI_Send(&sendbuf[0], sendbuf.size(), MPI_PACKED, parent, 0, comm);

    } else {
      // put packed data in the destination.
      dest.swap(sendbuf);
    }
  }


  ///
  /// Unpacks a packed vector in binomial order into objects in rank order in the destination vector.
  ///
  template <class T>
  void unpack_binomial(const std::vector<char>& src, std::vector<T>& dest, const binomial_embedding binomial, 
                       MPI_Comm comm) {
    int pos = 0;
    dest.resize(binomial.size());
    for (size_t i=0; i < binomial.size(); i++) {
      dest[binomial.reverse_relative_rank(i)] = T::unpack(const_cast<char*>(&src[0]), src.size(), &pos, comm);
    }
  }
  

  ///
  /// Binomial gather of char buffers into a single agglomerated clump of buffers
  ///
  template <class T>
  void gather(const T& src, std::vector<T>& dest, MPI_Comm comm, int root = 0) {
    int rank, size;
    CMPI_Comm_rank(comm, &rank);
    CMPI_Comm_size(comm, &size);

    // gather everything to a packed buffer at the root.
    binomial_embedding binomial(size, root);
    std::vector<char> packed;
    gather_packed(src, packed, binomial, comm);

    // now unpack everything.
    if (rank == root) {
      unpack_binomial(packed, dest, binomial, comm);
    }
  }


  ///
  /// Allgather for variable-length data.
  ///
  template <class T>
  void allgather(const T& src, std::vector<T>& dest, MPI_Comm comm, int root = 0) {
    int rank, size;
    CMPI_Comm_rank(comm, &rank);
    CMPI_Comm_size(comm, &size);

    // gather everything to a packed buffer at the root.
    binomial_embedding binomial(size, root);
    std::vector<char> packed;
    gather_packed(src, packed, binomial, comm);

    size_t packed_size = packed.size();
    CMPI_Bcast(&packed_size, 1, MPI_SIZE_T, root, comm);

    packed.resize(packed_size);
    CMPI_Bcast(const_cast<char*>(&packed[0]), packed.size(), MPI_PACKED, root, comm);

    unpack_binomial(packed, dest, binomial, comm);
  }

} // namespace cluster

#endif // MUSTER_GATHER_H
