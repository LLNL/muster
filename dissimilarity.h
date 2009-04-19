#ifndef DISSIMILARITY_MATRIX_H
#define DISSIMILARITY_MATRIX_H

#include <vector>
#include <boost/numeric/ublas/matrix.hpp>

typedef boost::numeric::ublas::matrix<double> dissimilarity_matrix;

///
/// Computes a dissimilarity matrix from a vector of objects.   Omits half the comparisons 
/// by assuming dissimilarity matrix is symmetric.  Also assumes that the dissimilarity 
/// between any object and itself is 0.
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
    mat(i,i) = 0;
    for (size_t j=i+1; j < objects.size(); j++) {
      mat(i,j) = mat(j,i) = dissimilarity(objects[i], objects[j]);
    }
  }
}


///
/// Computes a dissimilarity matrix from a subset of a vector of objects.
/// Makes the same assumptions as the exhaustive version.
///
/// Parameters:
///     objects         Vector of any type T.
///     subset          Indices into objects that we need the distance between.
///     dissimilarity   A dissimilarity measure that gives the distance between two T's.
///                     Needs to be callable(T, T).
///     mat             Output parameter.  Dissimiliarity matrix is stored here.
template <class T, class D>
void build_dissimilarity_matrix(const std::vector<T>& objects, const std::vector<size_t>& subset,
                                D dissimilarity, dissimilarity_matrix& mat) {
  mat.resize(subset.size(), subset.size());

  for (size_t i=0; i < subset.size(); i++) {
    mat(i,i) = 0;
    for (size_t j=i+1; j < subset.size(); j++) {
      mat(i,j) = mat(j,i) = dissimilarity(objects[subset[i]], objects[subset[j]]);
    }
  }
}

#endif // DISSIMILARITY_MATRIX_H
