#ifndef GAUSSIAN_GENERATOR_H
#define GAUSSIAN_GENERATOR_H

#include <cfloat>
#include <sys/time.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

#include "point.h"

namespace cluster {
  ///
  /// Gaussian cluster with cutoff radius.
  ///
  class gaussian_generator_2d {
  public:    
    gaussian_generator_2d(int _x, int _y, double stddev, double xscale=1.0, double yscale=1.0) 
      : x(_x), y(_y), xs(xscale), ys(yscale), normal(0, stddev), radius(rand, normal) { 
      struct timeval time;
      gettimeofday(&time, 0);
      rand.seed(time.tv_sec * time.tv_usec);
    }
    
    point next_point() {
      double r = radius();
      
      double theta = rand() / (double)rand.max() * 2 * M_PI;
      double px = x + r * cos(theta);
      double py = y + r * sin(theta);
      
      return point(round(px * xs), round(py * ys));
    }

  private:
    int x, y, xs, ys;
    boost::mt19937 rand;
    
    typedef boost::normal_distribution<double> dist_type;
    dist_type normal;
    boost::variate_generator<boost::mt19937, dist_type> radius;
  };


} // namespace

#endif // GAUSSIAN_GENERATOR_H
