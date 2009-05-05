#include "partition.h"
#include <iostream>
using namespace cluster;
using namespace std;

#include "counter.h"

namespace cluster {

  partition::partition(size_t num_objects) : cluster_ids(num_objects, 0) {
    if (num_objects) medoids.resize(1);
  }


  partition::~partition() { }
  

  void partition::to_cluster_list(cluster_list& clusters) const {
    clusters.clear();
    clusters.resize(medoids.size(), cset());
    for (unsigned object=0; object < cluster_ids.size(); object++) {
      clusters[cluster_ids[object]].insert(object);
    }
  }


  void partition::swap(partition& other) {
    medoids.swap(other.medoids);
    cluster_ids.swap(other.cluster_ids);
  }


  
  ostream& operator<<(ostream& out, const cluster_list& clusters) {
    out << "id\tmembers" << endl;
    for (unsigned i=0; i < clusters.size(); i++) {
      out << i << "\t";
      const cset& c = clusters[i];
      for (cset::iterator obj=c.begin(); obj != c.end(); obj++) {
        out << *obj << " ";
      }
      out << endl;
    }
    return out;
  }


  double mirkin_distance(const cluster_list& c1, const cluster_list& c2) {
    assert(c1.size() == c2.size());

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
    out << list;
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



} // namespace cluster
