#include "partition.h"
#include <iostream>
#include <iterator>
#include <map>
#include <cassert>
#include <algorithm>
using namespace cluster;
using namespace std;

#include "counter.h"
#include "color.h"
#include "stl_utils.h"

namespace cluster {

  partition::partition(size_t num_objects) : cluster_ids(num_objects, 0) {
    if (num_objects) medoid_ids.resize(1);
  }


  partition::~partition() { }
  

  void partition::to_cluster_list(cluster_list& clusters) const {
    clusters.clear();
    clusters.resize(medoid_ids.size(), cset());
    for (unsigned object=0; object < cluster_ids.size(); object++) {
      clusters[cluster_ids[object]].insert(object);
    }
  }


  void partition::swap(partition& other) {
    medoid_ids.swap(other.medoid_ids);
    cluster_ids.swap(other.cluster_ids);
  }
  
  
  void partition::sort() {
    // first create a mapping from new ids to old ids.
    vector<size_t> mapping(medoid_ids.size());
    generate(mapping.begin(), mapping.end(), sequence());

    std::sort(mapping.begin(), mapping.end(), indexed_lt(medoid_ids));
    invert(mapping);
    std::sort(medoid_ids.begin(), medoid_ids.end());

    // translate old cluster ids to new ones.
    for (size_t i=0; i < cluster_ids.size(); i++) {
      cluster_ids[i] = mapping[cluster_ids[i]];
    }
  }


  static void write(ostream& out, const cluster_list& clusters, const vector<object_id> *medoid_ids = NULL) {
    if (medoid_ids) {
      out << "Medoids: ";
      copy(medoid_ids->begin(), medoid_ids->end(), ostream_iterator<object_id>(out, " "));
      out << endl;
    }

    for (unsigned i=0; i < clusters.size(); i++) {
      out << i << "\t";
      const cset& c = clusters[i];
      for (cset::iterator obj=c.begin(); obj != c.end(); obj++) {
        if (medoid_ids && (*medoid_ids)[i] == *obj) {
          out << Red << *obj << None << " ";
        } else {
          out << *obj << " ";          
        }
      }
      out << endl;
    }
  }


  ostream& operator<<(ostream& out, const cluster_list& clusters) {
    write(out, clusters);
    return out;
  }


  size_t partition::size(size_t i) const {
    return count(cluster_ids.begin(), cluster_ids.end(), i);
  }



  double mirkin_distance(const cluster_list& c1, const cluster_list& c2) {
    size_t c1_sum2 = 0;
    size_t n = 0;
    for (size_t i=0; i < c1.size(); i++) {
      c1_sum2 += c1[i].size() * c1[i].size();
      n += c1[i].size();
    }

    size_t c2_sum2 = 0;
    for (size_t i=0; i < c2.size(); i++) {
      c2_sum2 += c2[i].size() * c2[i].size();    
    }
  

    size_t c1c2_sum2 = 0;
    for (size_t i=0; i < c1.size(); i++) {
      for (size_t j=0; j < c2.size(); j++) {
        size_t size;
        set_intersection(c1[i].begin(), c1[i].end(), 
                         c2[j].begin(), c2[j].end(),
                         counter<unsigned>(size));
        c1c2_sum2 += size * size;
      }
    }

    return (c1_sum2 + c2_sum2 - (2 * c1c2_sum2)) / (double)(n*n);
  }


  std::ostream& operator<<(std::ostream& out, const partition& p) {
    cluster_list list;
    p.to_cluster_list(list);
    write(out, list, &p.medoid_ids);
    return out;
  }


  double mirkin_distance(partition& c1, partition& c2) {
    cluster_list l1, l2;
    c1.to_cluster_list(l1);
    c2.to_cluster_list(l2);
    return mirkin_distance(l1, l2);
  }


  void expand(cluster_list& list, size_t level) {
    if (!level) return;

    cluster_list expanded(list.size());
    for (size_t i=0; i < list.size(); i++) {
      for (cset::iterator o=list[i].begin(); o != list[i].end(); o++) {
        size_t start = (1 << level) * (*o);
        size_t end   = (1 << level) * (*o + 1);
        for (size_t ex=start; ex < end; ex++) {
          expanded[i].insert(ex);
        }
      }
    }
    expanded.swap(list);
  }

  
  void partition::write_members_with_runs(medoid_id m, ostream& out) {
    object_id o=0;
    bool first = true;
    while (o < cluster_ids.size()) {
      if (cluster_ids[o] == m) {
        object_id start = o++;
        while (o < cluster_ids.size() && cluster_ids[o] == m) o++;
        if (!first) out << " ";
        if (o == start+1) {
          out << start;
        } else {
          out << start << "-" << (o-1);
        }
        first = false;
      }
      o++;
    }
  }

} // namespace cluster
