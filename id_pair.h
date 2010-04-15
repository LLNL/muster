#ifndef ID_PAIR_H
#define ID_PAIR_H

#include "libra-config.h"

#ifdef LIBRA_HAVE_MPI
#include <mpi.h>
#endif // LIBRA_HAVE_MPI



#include <cstdlib>
#include <ostream>

namespace cluster {

  ///
  /// Packable struct for a packable type plus its object id.
  ///
  template <class T>
  struct id_pair {
    T element;
    size_t id;

    /// Template typedef for creating vectors of id_pair<T>
    typedef std::vector< id_pair<T> > vector;

    id_pair() { }
    id_pair(const T& elt, size_t _id) : element(elt), id(_id) { }

#ifdef LIBRA_HAVE_MPI
    int packed_size(MPI_Comm comm) const {
      return element.packed_size(comm) + mpi_packed_size(1, MPI_SIZE_T, comm);
    }

    void pack(void *buf, int bufsize, int *position, MPI_Comm comm) const {
      element.pack(buf, bufsize, position, comm);
      MPI_Pack(const_cast<size_t*>(&id), 1, MPI_SIZE_T, buf, bufsize, position, comm);
    }

    static id_pair unpack(void *buf, int bufsize, int *position, MPI_Comm comm) {
      T t = T::unpack(buf, bufsize, position, comm);
      size_t id;
      MPI_Unpack(buf, bufsize, position, &id, 1, MPI_SIZE_T, comm);
      return id_pair(t, id);
    }
#endif // LIBRA_HAVE_MPI
  };
  
  /// Helper function for making id_pairs with type inference.
  template <class T>
  id_pair<T> make_id_pair(const T& elt, int id) {
    return id_pair<T>(elt, id);
  }
  
  /// print out an id_pair
  template <class T>
  std::ostream& operator<<(std::ostream& out, const id_pair<T>& p) {
    out << "<" << p.element << ", " << p.id << ">";
    return out;
  }
  
  
  
  
} // namespace cluster

#endif // ID_PAIR_H
