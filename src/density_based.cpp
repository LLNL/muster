
#include "density_based.h"

#include <vector>
#include <sstream>

#include <algorithm>
#include <numeric>
#include <cassert>
#include <cstdlib>
#include <sys/time.h>
using namespace std;

#include "random.h"


namespace cluster {
  
  density_based::density_based(size_t num_objects) 
    : partition(num_objects), 
      random(get_time_seed()),
      rng(random),
      current_cluster_id(cluster::first_cluster),
      total_clusters(0)
  {
  }


  density_based::~density_based() {  }

} // namespace cluster  
