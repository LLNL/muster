#ifndef COUNTER_ITERATOR_H
#define COUNTER_ITERATOR_H

#include <cstdlib>

/// Counting output iterator.  Stores (in your own size_t somewhere) how
/// many elements were output, but assignment is a no-op.  Usage:
///
/// size_t count;
/// vector<int> test;
/// test.push_back(1);
/// test.push_back(2);
/// test.push_back(3);
/// copy(test.begin(), test.end(), counter<int>(count));
/// 
/// POST: count is 3, since 3 items were inserted.
///
template <class T>
struct counter {
  typedef T                   value_type;
  typedef T*                  pointer;
  typedef T&                  reference;
  typedef size_t              difference_type;
  typedef output_iterator_tag iterator_category;

  struct target {
    void operator=(T t) { }  // no-op; assignment to target does nothing.
  };

  size_t *count;
  counter(size_t& c) : count(&c) { *count = 0; }
  counter(const counter& other) : count(other.count) { }
  counter& operator=(const counter& other) { count = other.count; return *this; }
  target operator*()       { return target(); }
  counter& operator++()    { (*count)++; return *this; }
  counter  operator++(int) { (*count)++; return *this; }
};

#endif // COUNTER_ITERATOR_H
