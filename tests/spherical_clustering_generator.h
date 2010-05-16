#ifndef SPHERICAL_CLUSTERING_GENERATOR_H
#define SPHERICAL_CLUSTERING_GENERATOR_H

#include <vector>
#include "point.h"
#include "gaussian.h"
#include "random.h"
#include <climits>


namespace cluster {
  
  class spherical_clustering_generator {
    std::vector<gaussian_generator_2d> generators;
    double default_stddev;
    double xscale;
    double yscale;
    size_t cur_generator;

  public:
    spherical_clustering_generator() 
      : default_stddev(1.0), xscale(1.0), yscale(1.0), cur_generator(0) { }

    void add_cluster(point center, double stddev = -1) {
      if (stddev < 0) stddev = default_stddev;
      generators.push_back(gaussian_generator_2d(center.x, center.y, stddev, xscale, yscale));
    }
    
    /// Get the next point from the next cluster in round-robin order.
    point next_point() {
      if (!generators.size()) {
        return point(0,0);
      }      
      point next = generators[cur_generator].next_point();
      cur_generator = (cur_generator + 1) % generators.size();
      return next;
    }

    /// Get the next point for a specific cluster
    point next_point(size_t i) {
      return generators[i].next_point();
    }
  
    void set_default_stddev(double dsd) {
      default_stddev = dsd;
    }
    
    void set_scale(double scale) {
      xscale = yscale = scale;
    }
    
    void set_xscale(double scale) {
      xscale = scale;
    }

    void set_yscale(double scale) {
      yscale = scale;
    }

    size_t size() const {
      return generators.size();
    }
    
  };

} // namespace cluster


#endif // SPHERICAL_CLUSTERING_GENERATOR_H
