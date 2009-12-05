#ifndef DISSIMILARITY_MATRIX_H
#define DISSIMILARITY_MATRIX_H

#include <vector>
#include <boost/numeric/ublas/symmetric.hpp>

namespace cluster {
  ///
  /// Packed repersentation of symmetric dissimilarity matrix.
  ///
  typedef boost::numeric::ublas::symmetric_matrix<double> dissimilarity_matrix;

  ///
  /// Computes a dissimilarity matrix from a vector of objects.
  ///
  /// Parameters:
  ///     objects         Vector of any type T.
  ///     dissimilarity   A dissimilarity measure that gives the distance between two T's.
  ///                     Needs to be callable(T, T).
  ///     mat             Output parameter.  Dissimiliarity matrix is stored here.
  template <class T, class D>
  void build_dissimilarity_matrix(const std::vector<T>& objects, D dissimilarity, 
                                  dissimilarity_matrix& mat) {
    mat.resize(objects.size(), objects.size());

    for (size_t i=0; i < objects.size(); i++) {
      for (size_t j=0; j <= i; j++) {
        mat(i,j) = dissimilarity(objects[i], objects[j]);
      }
    }
  }


  ///
  /// Computes a dissimilarity matrix from a subset of a vector of objects.
  ///
  /// Parameters:
  ///     objects         Vector of any type T.
  ///     subset          Indirection vector.  Contains indices into objects for 
  ///                     elements to be compared.
  ///     dissimilarity   A dissimilarity measure that gives the distance between two T's.
  ///                     Needs to be callable(T, T).
  ///     mat             Output parameter.  Dissimiliarity matrix is stored here.
  template <class T, class D>
  void build_dissimilarity_matrix(const std::vector<T>& objects, const std::vector<size_t>& subset,
                                  D dissimilarity, dissimilarity_matrix& mat) {
    mat.resize(subset.size(), subset.size());

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
