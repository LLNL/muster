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
/// @file mpi_packable_serializer.h
/// @author Todd Gamblin tgamblin@llnl.gov
/// @brief Simple serializer for classes that already support packed_size, pack, and unpack operations.
/// 
#ifndef MUSTER_MPI_PACKABLE_SERIALIZER_H
#define MUSTER_MPI_PACKABLE_SERIALIZER_H

namespace cluster {

  ///
  /// This is a simple serializer for classes that support pack, packed_size, and unpack
  /// operations.  The methods here simply delegate to the appropriate methods on T.
  /// 
  /// @note that this does NOT work well with polymorphic classes, as it 
  /// requires that T support a static unpack method.
  /// 
  /// @tparam T Type of objects to be serialized.
  ///   T must support the following operations:
  ///   - <code>int packed_size(MPI_Comm comm) const</code>
  ///   - <code>void pack(void *buf, int bufsize, int *position, MPI_Comm comm) const</code>
  ///   - <code>static T unpack(void *buf, int bufsize, int *position, MPI_Comm comm)</code>
  ///
  template <class T>
  class mpi_packable_serializer {
  public:
    int packed_size(const T& object, MPI_Comm comm) const {
      return object.packed_size(comm);
    }

    void pack(const T& object, void *buf, int bufsize, int *position, MPI_Comm comm) const {
      object.pack(buf, bufsize, position, comm);
    }

    void unpack(T& object, void *buf, int bufsize, int *position, MPI_Comm comm) const {
      object = T::unpack(buf, bufsize, position, comm);
    }
  };

}

#endif // MUSTER_MPI_PACKABLE_SERIALIZER_H
