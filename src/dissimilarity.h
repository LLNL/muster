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
/// @file dissimilarity.h
/// @author Todd Gamblin tgamblin@llnl.gov
/// @brief Data types and functions for dealing with dissimilarity matrices.
///
#ifndef DISSIMILARITY_MATRIX_H
#define DISSIMILARITY_MATRIX_H

#include <vector>
#include <boost/numeric/ublas/symmetric.hpp>
#include <iostream>

namespace cluster {
  ///
  /// Packed repersentation of symmetric dissimilarity matrix.
  ///
  typedef boost::numeric::ublas::symmetric_matrix<double> dissimilarity_matrix;

  ///
  /// Computes a dissimilarity matrix from a vector of objects.
  ///
  /// @param[in]  objects         Vector of any type T.
  /// @param[in]  dissimilarity   A dissimilarity measure that gives the distance between two T's.
  ///                             Needs to be callable on (T, T).
  /// @param[out] mat             Output parameter.  Dissimiliarity matrix is stored here.
  /// 
  template <class T, class D>
  void build_dissimilarity_matrix(const std::vector<T>& objects, D dissimilarity, 
                                  dissimilarity_matrix& mat) {
    if (mat.size1() != objects.size() || mat.size2() != objects.size()) {
      mat.resize(objects.size(), objects.size());
    }
    
    for (size_t i=0; i < objects.size(); i++) {
      for (size_t j=0; j <= i; j++) {
        mat(i,j) = dissimilarity(objects[i], objects[j]);
      }
    }
  }


  ///
  /// Computes a dissimilarity matrix from a subset of a vector of objects.
  ///
  /// @param objects         Vector of any type T.
  /// @param subset          Indirection vector.  Contains indices into objects for 
  ///                        elements to be compared.
  /// @param dissimilarity   A dissimilarity measure that gives the distance between two T's.
  ///                        Needs to be callable(T, T).
  /// @param mat             Output parameter.  Dissimiliarity matrix is stored here.
  template <class T, class D>
  void build_dissimilarity_matrix(const std::vector<T>& objects, const std::vector<size_t>& subset,
                                  D dissimilarity, dissimilarity_matrix& mat) {
    if (mat.size1() != subset.size() || mat.size2() != subset.size()) {
      mat.resize(subset.size(), subset.size());
    }

    for (size_t i=0; i < subset.size(); i++) {
      for (size_t j=0; j <= i; j++) {
        mat(i,j) = dissimilarity(objects[subset[i]], objects[subset[j]]);
      }
    }
  }

  
  /// Adaptor for passing a matrix by reference to template functions that take
  /// a callable distance function.  Avoids copying distance matrix.
  struct matrix_distance {
    const dissimilarity_matrix& mat;
    matrix_distance(const dissimilarity_matrix& m) : mat(m) { }
    double operator()(size_t i, size_t j) { return mat(i,j); }
  };


  /// Functor for computing distance lazily from an object array and
  /// a distance metric.  Use this for CLARA, where we don't want to
  /// precompute the entire distance matrix.
  template <class T, class D>
  struct lazy_distance_functor {
    const std::vector<T>& objects;
    D dissimilarity;

    lazy_distance_functor(const std::vector<T>& objs, D d)
      : objects(objs), dissimilarity(d) { }

    double operator()(size_t i, size_t j) {
      return dissimilarity(objects[i], objects[j]);
    }
  };

  /// Type-inferred syntactic sugar for constructing lazy_distance_functor.
  template <class T, class D>
  lazy_distance_functor<T,D> lazy_distance(const std::vector<T>& objs, D dist) {
    return lazy_distance_functor<T,D>(objs, dist);
  }

}; // namespace cluster

#endif // DISSIMILARITY_MATRIX_H
