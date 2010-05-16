#ifndef COUNTER_ITERATOR_H
#define COUNTER_ITERATOR_H
/// 
/// @file counter.h
/// @brief Dummy output iterator that counts how many times it was assigned to.
///        without actually storing anything.
///

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
