#ifndef STL_UTILS
#define STL_UTILS

#include <utility>
#include <vector>
#include <iostream>

// ====================================================================================
// stl_utils.h
// This file contains some utility functions for dealing with the STL.
// ====================================================================================


/// Functor to get the first element of a pair.  Use with STL functions like transform().
struct get_first {
  template <typename P>
  typename P::first_type operator()(const P& pair) {
    return pair.first;
  }
};

/// Functor to get the second element of a pair.  Use with STL functions like transform().
struct get_second {
  template <typename P>
  typename P::second_type operator()(const P& pair) {
    return pair.second;
  }
};


template <typename Indexable>
struct indexed_lt_functor {
  const Indexable& container;
  indexed_lt_functor(const Indexable& c) : container(c) { }
  template <typename I>
  bool operator()(const I& lhs, const I& rhs) {
    return container[lhs] < container[rhs];
  }
};

template <typename Indexable>
indexed_lt_functor<Indexable> indexed_lt(const Indexable& container) {
  return indexed_lt_functor<Indexable>(container);
}


template <typename Index>
void invert(std::vector<Index>& vec) {
  std::vector<Index> inverse(vec.size());
  for (size_t i=0; i < vec.size(); i++) {
    inverse[vec[i]] = i;
  }
  inverse.swap(vec);
}


///
/// Generator object for a strided sequence of ints.
///
struct sequence {
  int value, stride;

  sequence(int _start=0, int _stride=1) 
    : value(_start), stride(_stride) { }

  int operator()() {
    int result = value;
    value += stride;
    return result;
  }
};




#endif // STL_UTILS
