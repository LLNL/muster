#ifndef EFFORT_RANDOM_H
#define EFFORT_RANDOM_H

#include "RNGenerator.h"

namespace effort {

  template <class OutputIterator>
  void randomSubset(size_t numElements, size_t sample_size, OutputIterator out) {
    static RNGenerator rnGenerator;

    // This is Knuth's algorithm R for taking a sample of indices from
    // 0 to numElements.  We sample size indices from this (the superset)
    // and put them in the subset's mapping.
    size_t first = 0;
    size_t remaining = numElements;
    size_t m = sample_size;
    
    while (m > 0) {
      if (rnGenerator(remaining) < m) {
        *out = first;
        ++out;
        --m;
      }
      --remaining;
      ++first;
    }
  }

} // namespace effort

#endif // EFFORT_RANDOM_H
