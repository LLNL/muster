#ifndef ID_PAIR_H
#define ID_PAIR_H
///
/// @file id_pair.h
/// @brief MPI-packable, templated struct for shipping around an MPI-packable
///        object plus the id of the process it came from.
///

#include <mpi.h>
#include "mpi_bindings.h"

#include <cstdlib>
#include <ostream>

namespace cluster {

  ///
  /// MPI-packable struct for an MPI-packable type plus its object id.
  /// 
  /// Each id_pair<T> has an element and an id for that element and supports
  /// packed_size(), pack(), and unpack() methods for transferring these 
  /// things with MPI.
  /// 
  /// @tparam T Type of contained element.  
  ///           T Must support MPI pack(), packed_size(), and unpack() methods.
  ///
  template <class T>
  struct id_pair {
    T element;    ///< The object wrapped by this id_pair.
    size_t id;    ///< Id of the rank where element came from.

    /// Template typedef for declaring vectors of id_pair<T>
    typedef std::vector< id_pair<T> > vector;

    id_pair() { }
    id_pair(const T& elt, size_t _id) : element(elt), id(_id) { }

    int packed_size(MPI_Comm comm) const {
      return element.packed_size(comm) + cmpi_packed_size(1, MPI_SIZE_T, comm);
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
  };
  
  ///
  /// Helper function for making arbitrary id_pairs with type inference.
  /// 
  template <class T>
  id_pair<T> make_id_pair(const T& elt, int id) {
    return id_pair<T>(elt, id);
  }
  
  ///
  /// Print out an id_pair as a tuple of its element and its source rank.
  /// 
  /// @tparam T inferred from the id_pair<T> this is called on.  Must support operator<<.
  /// 
  /// @param out Output stream to write the id_pair p onto
  /// @param p   An id_pair<T>.  T must support operator<<.
  /// 
  template <class T>
  std::ostream& operator<<(std::ostream& out, const id_pair<T>& p) {
    out << "<" << p.element << ", " << p.id << ">";
    return out;
  }
    
} // namespace cluster

#endif // ID_PAIR_H
