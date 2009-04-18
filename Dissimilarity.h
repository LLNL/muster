#ifndef DISSIMILARITY_H
#define DISSIMILARITY_H

#include <float.h>
#define DISSIMILARITY_MAX DBL_MAX

/**
 * Dissimilarity measure for Kmedoids clustering algorithms (KMedoids, CLARA, etc).
 */
template <class T>
struct Dissimilarity {
  /**
   * Takes the dissimilarity between two objects
   */
  virtual double getDissimilarity(T left, T right) const = 0;
  virtual ~Dissimilarity() { }
};

#endif // DISSIMILARITY_H
