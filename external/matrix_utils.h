#ifndef MATRIX_UTILS_H
#define MATRIX_UTILS_H

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <cfloat>
#include <algorithm>
#include <numeric>
#include <stdint.h>
#include <boost/numeric/ublas/matrix.hpp>


/// True if and only if n is divisible by 2 <level> times.
bool isDivisibleBy2(size_t n, int level);


template <class Matrix>
bool in_bounds(const Matrix& mat, size_t row, size_t col) {
  return (row < mat.size1()) && (col < mat.size2());
}


bool read_matrix(const char *filename, boost::numeric::ublas::matrix<double>& mat);


template <class Matrix>
void output(const Matrix& mat, std::ostream& out = std::cout) {
  const size_t width = 12;

  for (size_t i=0; i < mat.size1(); i++) {
    if (mat.size2()) {
      out << std::setw(width) << mat(i,0);
    }

    for (size_t j=1; j < mat.size2(); j++) {
      out << " " << std::setw(width) << mat(i,j);
    }
    out << std::endl;
  }
}


struct ms_summary {
  double sum_squares;
  double orig_max, orig_min;
  double repro_max, repro_min;
  double both_max, both_min;

  ms_summary(double ss, double ma, double mi) 
    : sum_squares(ss), 
      orig_max(ma),  orig_min(mi), 
      repro_max(ma), repro_min(mi),
      both_max(ma),  both_min(mi) 
  { }
};


template <class Matrix>
ms_summary get_summary(const Matrix& orig, const Matrix& repro,
                       size_t row_start = 0, size_t row_end = std::numeric_limits<size_t>::max(),
                       size_t col_start = 0, size_t col_end = std::numeric_limits<size_t>::max()) {
  assert(orig.size1() == repro.size1());
  assert(orig.size2() == repro.size2());

  if (row_end > orig.size1()) row_end = orig.size1();
  if (col_end > orig.size2()) col_end = orig.size2();

  ms_summary summary(0, -DBL_MAX, DBL_MAX);

  for (size_t i=row_start; i < row_end; i++) {
    for (size_t j=col_start; j < col_end; j++) { 
      double diff = repro(i,j) - orig(i,j);
      summary.sum_squares += diff*diff;

      summary.orig_max =  std::max(summary.orig_max, orig(i,j));
      summary.orig_min =  std::min(summary.orig_min, orig(i,j));
      summary.repro_max = std::max(summary.repro_max, repro(i,j));
      summary.repro_min = std::min(summary.repro_min, repro(i,j));
      summary.both_max =  std::max(summary.orig_max, summary.repro_max);
      summary.both_min =  std::min(summary.orig_min, summary.repro_min);
    }
  }

  return summary;
}


template <class Matrix>
double rmse(const Matrix& orig, const Matrix& repro,
             size_t row_start = 0, size_t row_end = std::numeric_limits<size_t>::max(),
             size_t col_start = 0, size_t col_end = std::numeric_limits<size_t>::max()) {
  assert(orig.size1() == repro.size1());
  assert(orig.size2() == repro.size2());

  if (row_end > orig.size1()) row_end = orig.size1();
  if (col_end > orig.size2()) col_end = orig.size2();

  ms_summary summary = get_summary(orig, repro, row_start, row_end, col_start, col_end);
  size_t rows = (row_end - row_start);
  size_t cols = (col_end - col_start);
  return sqrt(summary.sum_squares / (rows * cols));
}


/// Normalized rms error
template <class Matrix>
double nrmse(const Matrix& orig, const Matrix& repro,
             size_t row_start = 0, size_t row_end = std::numeric_limits<size_t>::max(),
             size_t col_start = 0, size_t col_end = std::numeric_limits<size_t>::max()) {
  assert(orig.size1() == repro.size1());
  assert(orig.size2() == repro.size2());

  if (row_end > orig.size1()) row_end = orig.size1();
  if (col_end > orig.size2()) col_end = orig.size2();

  ms_summary summary = get_summary(orig, repro, row_start, row_end, col_start, col_end);
  double range = summary.orig_max - summary.orig_min;
  size_t rows = (row_end - row_start);
  size_t cols = (col_end - col_start);
  return sqrt(summary.sum_squares / (rows * cols)) / range;
}


/// Peak Signal to Noise Ratio
template <class Matrix>
double psnr(const Matrix& orig, const Matrix& repro) {
  assert(orig.size1() == repro.size1());
  assert(orig.size2() == repro.size2());

  ms_summary summary = get_summary(orig, repro);
  double range = summary.orig_max - summary.orig_min;
  return 20 * log10(range / sqrt(summary.sum_squares / (orig.size1() * orig.size2())));
}


/// Similarity (NRMSE defined to be symmetric, with max and min taken from both matrices)
template <class Matrix>
double similarity(const Matrix& orig, const Matrix& repro) {
  assert(orig.size1() == repro.size1());
  assert(orig.size2() == repro.size2());

  ms_summary summary = get_summary(orig, repro);
  double range = summary.both_max - summary.both_min;
  return sqrt(summary.sum_squares / (orig.size1() * orig.size2())) / range;
}


template <class Matrix>
void standardize(Matrix& mat) {
  size_t n = mat.size1() * mat.size2();
  double sum = 0;
  double sum2 = 0;
  for (size_t i=0; i < mat.size1(); i++) {
    for (size_t j=0; j < mat.size2(); j++) {
      sum += mat(i,j);
      sum2 += mat(i,j) * mat(i,j);
    }
  }
  double mean = sum / n;
  double stdDevInv = 1/sqrt(sum2/n - mean*mean);
  
  for (size_t i=0; i < mat.size1(); i++) {
    for (size_t j=0; j < mat.size2(); j++) {
      mat(i,j) = (mat(i,j) - mean) * stdDevInv;
    }
  }
}


/// Template class for generic, type-inferred absolute value 
template <typename T>
struct Abs {
  T value;
  Abs(double num)        : value(::fabs(num)) { }
  Abs(long num)          : value(::labs(num)) { }
  Abs(unsigned long num) : value(::labs(num)) { }
  Abs(unsigned int num)  : value(::labs(num)) { }
  Abs(long long num)     : value(::llabs(num)) { }
  Abs(int num)           : value(::abs(num)) { }
};


/// Generic, type-inferred absolute value function.
/// Applies fabs to doubles, labs to longs, etc.
template<typename T> T abs_val(T num) {
  return Abs<T>(num).value;
}



template <class Matrix>
typename Matrix::value_type sum(const Matrix& mat,
      size_t row_start = 0, size_t row_end = std::numeric_limits<size_t>::max(),
      size_t col_start = 0, size_t col_end = std::numeric_limits<size_t>::max())
{
  if (row_end > mat.size1()) row_end = mat.size1();
  if (col_end > mat.size2()) col_end = mat.size2();

  typename Matrix::value_type total = 0;
  for (size_t i = row_start; i < row_end; i++) {
    for (size_t j = col_start; j < col_end; j++) { 
      total += mat(i,j);
    }
  }
  return total;
}


template <class Matrix>
double mean_val(const Matrix& mat,
                size_t row_start = 0, size_t row_end = std::numeric_limits<size_t>::max(),
                size_t col_start = 0, size_t col_end = std::numeric_limits<size_t>::max())
{
  if (row_end > mat.size1()) row_end = mat.size1();
  if (col_end > mat.size2()) col_end = mat.size2();

  size_t row_len = (row_end - row_start);
  size_t col_len = (col_end - col_start);
  return sum(mat, row_start, row_end, col_start, col_end) 
    / ((double)row_len * col_len);
}


template <class Matrix>
typename Matrix::value_type max_val(const Matrix& mat,
          size_t row_start = 0, size_t row_end = std::numeric_limits<size_t>::max(),
          size_t col_start = 0, size_t col_end = std::numeric_limits<size_t>::max()) 
{
  if (row_end > mat.size1()) row_end = mat.size1();
  if (col_end > mat.size2()) col_end = mat.size2();
  if (row_end <= row_start || col_end <= col_start) return 0;

  typename Matrix::value_type max = mat(row_start, col_start);
  for (size_t i = row_start; i < row_end; i++) {
    for (size_t j = col_start; j < col_end; j++) { 
      max = mat(i,j) > max ? mat(i,j) : max;
    }
  }
  return max;
}


template <class Matrix>
typename Matrix::value_type min_val(const Matrix& mat,
          size_t row_start = 0, size_t row_end = std::numeric_limits<size_t>::max(),
          size_t col_start = 0, size_t col_end = std::numeric_limits<size_t>::max()) 
{
  if (row_end > mat.size1()) row_end = mat.size1();
  if (col_end > mat.size2()) col_end = mat.size2();
  if (row_end <= row_start || col_end <= col_start) return 0;

  typename Matrix::value_type min = mat(row_start, col_start);
  for (size_t i = row_start; i < row_end; i++) {
    for (size_t j = col_start; j < col_end; j++) { 
      min = mat(i,j) < min ? mat(i,j) : min;
    }
  }
  return min;
}

template <class Matrix>
typename Matrix::value_type abs_max_val(const Matrix& mat,
          size_t row_start = 0, size_t row_end = std::numeric_limits<size_t>::max(),
          size_t col_start = 0, size_t col_end = std::numeric_limits<size_t>::max()) 
{
  if (row_end > mat.size1()) row_end = mat.size1();
  if (col_end > mat.size2()) col_end = mat.size2();

  typename Matrix::value_type amax = 0;
  for (size_t i = row_start; i < row_end; i++) {
    for (size_t j = col_start; j < col_end; j++) { 
      typename Matrix::value_type tmp = abs_val(mat(i,j));
      amax = (tmp > amax) ? tmp : amax;
    }
  }
  return amax;
}

template <class Matrix, typename V>
void set_all(const Matrix& mat, V value,
          size_t row_start = 0, size_t row_end = std::numeric_limits<size_t>::max(),
          size_t col_start = 0, size_t col_end = std::numeric_limits<size_t>::max()) 
{
  if (row_end > mat.size1()) row_end = mat.size1();
  if (col_end > mat.size2()) col_end = mat.size2();

  for (size_t i = row_start; i < row_end; i++) {
    for (size_t j = col_start; j < col_end; j++) { 
      mat(i,j) = value;
    }
  }
}


template <typename Iterator1, typename Iterator2>
double manhattan_distance(Iterator1 first1, Iterator1 last1, Iterator2 first2) {
  double sum = 0.0;
  Iterator1 i = first1;
  Iterator2 j = first2;
  while (i != last1) {
    sum += abs_val(*i - *j);
    i++;
    j++;
  }
  return sum;
}


template <typename Iterator1, typename Iterator2>
double euclidean_distance(Iterator1 first1, Iterator1 last1, Iterator2 first2) {
  double sum = 0.0;
  Iterator1 i = first1;
  Iterator2 j = first2;
  while (i != last1) {
    double diff = (*i - *j);
    sum += diff * diff;
    i++;
    j++;
  }
  return sqrt(sum);
}


/// Interpolates a value for point (x,y) based on values in the Matrix.
/// Finds the nearest values by taking floor and ceil of x and y, then
/// uses bilinear interpolation to estimate the value at (x,y).
template <class Matrix>
double interp_bilinear(const Matrix& mat, double x, double y) {
  if (x < 0) {
    x = 0;
  } else if (x > (mat.size1() - 1)) {
    x = mat.size1() - 1;
  }

  if (y < 0) {
    y = 0;
  } else if (y > (mat.size2() - 1)) {
    y = mat.size2() - 1;
  }
  
  // get corners of interpolation
  size_t x1 = (size_t)floor(x);
  size_t x2 = (size_t)ceil(x);

  size_t y1 = (size_t)floor(y);
  size_t y2 = (size_t)ceil(y);

  size_t xd = (x2-x1);
  size_t yd = (y2-y1);
  
  double value;
  if (xd == 0 && yd == 0) {
    value = mat(x1,y1);

  } else if (yd == 0) {
    value = (x2 - x) * mat(x1,y1) + (x - x1) * mat(x2, y2);
    
  } else if (xd == 0) {
    value = (y2 - y) * mat(x1,y1) + (y - y1) * mat(x2, y2);
    
  } else {
    double d = xd * yd;
    value =  
      mat(x1,y1)/d * ((x2-x) * (y2-y)) +
      mat(x2,y1)/d * ((x-x1) * (y2-y)) +
      mat(x1,y2)/d * ((x2-x) * (y-y1)) +
      mat(x2,y2)/d * ((x-x1) * (y-y1));
  }
  
  return value;
}

#endif // MATRIX_UTILS_H
