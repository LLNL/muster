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
    typedef boost::mt19937 rand_type;

  public:    
    gaussian_generator_2d(double _x, double _y, double stddev, double xscale=1.0, double yscale=1.0) 
      : x(_x * xscale), y(_y * yscale), xs(xscale), ys(yscale), normal(0, stddev), radius(rand, normal) { 
      rand.seed(get_time_seed());
    }

    void seed(rand_type::result_type s) {
      rand.seed(s);
    }

    point next_point() {
      double r = radius();
      double theta = rand() / (double)rand.max() * 2 * M_PI;
      double px = x + r * cos(theta);
      double py = y + r * sin(theta);
      
      return point(px * xs, py * ys);
    }

    point center() {
      return point(x,y);
    }

  private:
    double x, y;
    double xs, ys;
    rand_type rand;
    
    typedef boost::normal_distribution<double> dist_type;
    dist_type normal;
    boost::variate_generator<rand_type, dist_type> radius;
  };


} // namespace

#endif // GAUSSIAN_GENERATOR_H
