#ifndef EFFORT_RANDOM_H
#define EFFORT_RANDOM_H

#include "RNGenerator.h"

template <class OutputIterator, class Random>
void random_subset(size_t numElements, size_t sample_size, OutputIterator out, 
                   Random& random = RNGenerator()) {
  // This is Knuth's algorithm R for taking a sample of indices from
  // 0 to numElements.  We sample size indices from this (the superset)
  // and put them in the subset's mapping.
  size_t first = 0;
  size_t remaining = numElements;
  size_t m = sample_size;
  
  while (m > 0) {
    if (random(remaining) < m) {
      *out = first;
      ++out;
      --m;
    }
    --remaining;
    ++first;
  }
}

#endif // EFFORT_RANDOM_H
