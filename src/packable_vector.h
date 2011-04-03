#ifndef MUSTER_PACKABLE_VECTOR_H
#define MUSTER_PACKABLE_VECTOR_H

#include <cstdlib>
#include <vector>
#include <boost/shared_ptr.hpp>
#include "mpi_utils.h"
#include "mpi_bindings.h"

namespace cluster {

  /// Pass this to a shared pointer if you do NOT want it to own its object
  struct null_deleter {
    void operator()(void const *) const { }
  };
  
  
  ///
  /// This class allows a vector of packable objects to be packed as though
  /// it were a packable object itself.
  ///
  template <class T>
  struct packable_vector {
    boost::shared_ptr< std::vector<T> > _packables;

    packable_vector(std::vector<T> *packables, bool owned = true) {
      if (owned) {
        _packables = boost::shared_ptr< std::vector<T> >(packables);
      } else {
        _packables = boost::shared_ptr< std::vector<T> >(packables, null_deleter());
      }
    }
    

    packable_vector(const packable_vector& other) : _packables(other._packables) { }
    packable_vector() : _packables(new std::vector<T>()) { }

    ~packable_vector() { }

    ///
    /// Assignment
    ///
    packable_vector& operator=(const packable_vector& other) {
      if (this == &other) return *this;
      _packables = other._packables;
      return *this;
    }

    ///
    /// get the number of bytes required to pack this buffer.
    ///
    int packed_size(MPI_Comm comm) const {
      // figure out size of packed buffer
      int size = 0;
      size += cmpi_packed_size(1, MPI_SIZE_T, comm);  // num packables for trial.
      for (size_t i=0; i < _packables->size(); i++) {          // size of packables.
        size += (*_packables)[i].packed_size(comm);
      }
      return size;
    }

    ///
    /// Pack onto an MPI buffer.
    ///
    void pack(void *buf, int bufsize, int *pos, MPI_Comm comm) const {
      // pack buffer with medoid objects
      size_t num_packables = _packables->size();
      CMPI_Pack(&num_packables, 1, MPI_SIZE_T, buf, bufsize, pos, comm);
      for (size_t i=0; i < num_packables; i++) {
        (*_packables)[i].pack(buf, bufsize, pos, comm);
      }
    }

    ///
    /// Unpack from an input buffer.  Note that this creates a new vector.
    ///
    static packable_vector unpack(void *buf, int bufsize, int *pos, MPI_Comm comm) {
      size_t num_packables;
      CMPI_Unpack(buf, bufsize, pos, &num_packables, 1, MPI_SIZE_T, comm);

      packable_vector vec;
      vec._packables->resize(num_packables);
      for (size_t i=0; i < vec._packables->size(); i++) {
        (*vec._packables)[i] = T::unpack(buf, bufsize, pos, comm);
      }
      return vec;
    }
  };
  
  ///
  /// Tempate function adaptor so we can have type inference when making 
  /// packable vectors.
  ///
  template <class T>
  packable_vector<T> make_packable_vector(std::vector<T> *vec, bool owned = true) {
    return packable_vector<T>(vec, owned);
  }
  
} // namespace cluster

#endif // MUSTER_PACKABLE_VECTOR_H
