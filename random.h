#ifndef EFFORT_RANDOM_H
#define EFFORT_RANDOM_H

///
/// This is Knuth's algorithm R for taking a sample of indices from
/// 0 to numElements.  We sample size indices from this (the superset)
/// and put them in the subset's mapping.
///
/// Parameters:
///    numElements    total number of elements to select from
///    sample_size    number of elements to select
///    out            destination for selected elements 
///    random         model of STL Random Number Generator.
///                   must be callable as random(N), returning a random number in [0,N).
///
template <class OutputIterator, class Random>
void random_subset(size_t numElements, size_t sample_size, OutputIterator out, Random& random) {
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
