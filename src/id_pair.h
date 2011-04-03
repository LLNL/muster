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
/// @file id_pair.h
/// @brief MPI-packable, templated struct for shipping around an MPI-packable
///        object plus the id of the process it came from.
///
#ifndef ID_PAIR_H
#define ID_PAIR_H

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
      CMPI_Pack(const_cast<size_t*>(&id), 1, MPI_SIZE_T, buf, bufsize, position, comm);
    }

    static id_pair unpack(void *buf, int bufsize, int *position, MPI_Comm comm) {
      T t = T::unpack(buf, bufsize, position, comm);
      size_t id;
      CMPI_Unpack(buf, bufsize, position, &id, 1, MPI_SIZE_T, comm);
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
