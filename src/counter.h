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
/// @file counter.h
/// @author Todd Gamblin tgamblin@llnl.gov
/// @brief Dummy output iterator that counts how many times it was assigned to.
///        without actually storing anything.
///
#ifndef COUNTER_ITERATOR_H
#define COUNTER_ITERATOR_H

#include <cstdlib>
#include <iterator>

namespace cluster {

  /// 
  /// Counting output iterator that records how many times an output iterator
  /// was assigned to, but ignores the value stored.  
  ///
  /// This is useful if you just want to know the size of something that an STL 
  /// algorithm would output, without actually allocating space for it.
  /// 
  /// @note
  /// Don't use this directly; Instantiate this using the counter() template function 
  /// so that you don't have to supply a type (see sample usage in file docs).
  ///
  template <class T>
  struct counter_iterator {
    typedef T                        value_type;
    typedef T*                       pointer;
    typedef T&                       reference;
    typedef size_t                   difference_type;
    typedef std::output_iterator_tag iterator_category;
    
    /// struct representation of a no-op.  Makes assignment to target do nothing.
    struct target {
      void operator=(T t) { }
    };
    
    pointer count;
    counter_iterator(value_type& c) : count(&c) { *count = 0; }
    counter_iterator(const counter_iterator& other) : count(other.count) { }
    counter_iterator& operator=(const counter_iterator& other) { count = other.count; return *this; }
    target operator*()              { return target(); }
    counter_iterator& operator++()    { (*count)++; return *this; }
    counter_iterator  operator++(int) { (*count)++; return *this; }
  };
  
  ///
  /// Adaptor for creating type-inferred counters.
  ///
  /// <b>Example Usage:</b>
  /// @code
  /// #include <algorithm>
  /// 
  /// // construct two sets
  /// set<int> s1, s2;
  /// 
  /// // insert some things so that their intersection has 2 ints.
  /// s1.insert(1); s1.insert(2); s1.insert(3);
  /// s2.insert(2); s2.insert(3); s2.insert(4);
  /// 
  /// // Compute intersection, but throw away the values and just 
  /// // store the size in the count variable.
  /// size_t count;
  /// set_intersection(s1.begin(), s1.end(), 
  ///                  s2.begin(), s2.end(), counter(count));
  ///
  /// // now count == 2, since 2 items were inserted 
  /// // by <code>set_intersection.
  ///
  /// @endcode
  ///
  template <class T>
  counter_iterator<T> counter(T& ref) {
    return counter_iterator<T>(ref);
  }

} // namespace cluster
  
#endif // COUNTER_ITERATOR_H
