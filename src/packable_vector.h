#ifndef MUSTER_PACKABLE_VECTOR_H
#define MUSTER_PACKABLE_VECTOR_H

#include <cstdlib>
#include <vector>
#include <boost/shared_ptr.hpp>

namespace muster {

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
      _packables(packables) ;
        
    }
    

    packable_vector(const packable_vector& other) : _packables(other._packables) { }
    packable_vector() : _packables(new vector<T>());

    ~packable_vector() { }

    ///
    /// Assignment
    ///
    packable_vector& operator=(const packable_vector& other) {
      if (this == &other) return *this;
      _packables = other.packables;
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
      CMPI_Pack(&num_packables, 1, MPI_SIZE_T, &buffer[0], packed_size, &pos, comm);
      for (size_t i=0; i < num_packables; i++) {
        (*_packables)[i].pack(&buffer[0], packed_size, &pos, comm);
      }
    }

    ///
    /// Unpack from an input buffer.  Note that this creates a new vector.
    ///
    static packable_vector unpack(void *buf, int bufsize, int *pos, MPI_Comm comm) {
      size_t num_packables;
      CMPI_Unpack((void*)&buffer[0], buffer.size(), &pos, &num_packables, 1, MPI_SIZE_T, comm);

      packable_vector vec;
      vec._packables->resize(num_packables);
      for (size_t i=0; i < _packables.size(); i++) {
        (*vec._packables)[i] = T::unpack((void*)&buffer[0], buffer.size(), &pos, comm);
      }
    }
  };
  
} // namespace muster

#endif // MUSTER_PACKABLE_VECTOR_H
