#ifndef PARTITION_H
#define PARTITION_H

#include <cstddef>
#include <vector>
#include <set>
#include <ostream>
#include <stdint.h>

namespace cluster {

  typedef std::vector< std::set<size_t> > cluster_list;
  typedef std::set<size_t> cset;
  
  typedef size_t medoid_id;       /// More descriptive type for medoid index
  typedef size_t object_id;       /// More descriptive type for object index

  /// Class to represent a partitioning of a dataset into 
  /// clusters with medoids.
  struct partition {
    /// Gives the index of the object that is the ith medoid.
    /// medoids[i] == index in object array for last call to findClusters()
    std::vector<object_id> medoid_ids;

    /// Gives cluster id (index in medoids) for the ith object.
    /// clusterid[i]          == id of cluster of which object i is a member.
    /// medoids[clusterid[i]] == representative of that cluster.
    std::vector<medoid_id> cluster_ids;
    
    /// Constructor.  Can optionall supply the number of objects to be partitioned
    /// and this will start out with one cluster containing all of them.
    partition(size_t num_objects = 0);
    virtual ~partition();

    /// True if and only if object i is a medoid.
    bool is_medoid(object_id oi) const {
      return medoid_ids[cluster_ids[oi]] == oi;
    }

    /// Creates a list of std::sets from the partition info in 
    /// medoids and cluster_ids.
    void to_cluster_list(cluster_list& list) const;

    /// Fast swap with other patrition objects
    void swap(partition& other);
    
    /// puts medoids in order by their object id, and adjusts cluster_ids accordingly.
    void sort();
    
    /// Total number of objects in the partition
    size_t size() const { return cluster_ids.size(); }

    /// Total number of clusters in the partition
    size_t num_clusters() const { return medoid_ids.size(); }
    
    /// Number of objects in cluster i
    size_t size(size_t i) const;
    
    /// Write the members of cluster m out to the output stream as object_ids
    template <class OutputIterator>
    void write_members(medoid_id m, OutputIterator out) {
      for (object_id o=0; o < cluster_ids.size(); o++) {
        if (cluster_ids[o] == m) {
          *out++ = o;
        }
      }
    }

    /// Write the members of cluster m out to the output stream formatted nicely with
    /// hyphenated runs of consecutive numbers
    void write_members_with_runs(medoid_id m, std::ostream& out);

    /// writable structure returned by members() function.
    struct member_writer {    
      partition* p;
      medoid_id m;
      member_writer(partition *_p, medoid_id _m) : p(_p), m(_m) { }
    };
    member_writer members(medoid_id m) { return member_writer(this, m); }
  };

  inline std::ostream& operator<<(std::ostream& out, const partition::member_writer& mw) {
    mw.p->write_members_with_runs(mw.m, out);
    return out;
  }
  
  /// Prints out nicely formatted clustering
  std::ostream& operator<<(std::ostream& out, const cluster_list& list);

  /// For convenience
  std::ostream& operator<<(std::ostream& out, const partition& km);  

  /// Mirkin distance bt/w two clusterings.
  double mirkin_distance(const cluster_list& c1, const cluster_list& c2);

  /// Convenience overload for comparing partition objects directly
  double mirkin_distance(partition& c1, partition& c2);

  /// Expand a cluster_list by l levels.  That is, replace each index i
  /// in the cluster_list with indices in [2^l * i ... 2^l * (i+1) - 1]
  /// TODO: deprecate and delete.
  void expand(cluster_list& list, size_t level = 1);

  ///
  /// Compute the total dissimilarity between all objects and their medoids.
  ///
  template <typename D>
  double total_dissimilarity(const partition& p, D dist) {
    double dissim = 0.0;
    for (size_t i=0; i < p.cluster_ids.size(); i++) {
      dissim += dist(i, p.medoid_ids[p.cluster_ids[i]]);
    }
    return dissim;
  }

  ///
  /// Compute the total dissimilarity between all objects in a particular cluster and its medoid.
  ///
  template <typename D>
  double total_dissimilarity(const partition& p, D dist, medoid_id m) {
    double dissim = 0.0;
    for (size_t i=0; i < p.cluster_ids.size(); i++) {
      if (p.cluster_ids[i] == m) {
        dissim += dist(i, p.medoid_ids[p.cluster_ids[i]]);
      }
    }
    return dissim;
  }
  
  ///
  /// Compute the total squared dissimilarity between all objects and their medoids.
  ///
  template <typename D>
  double total_squared_dissimilarity(const partition& p, D dist) {
    double dissim2 = 0.0;
    for (size_t i=0; i < p.cluster_ids.size(); i++) {
      double d = dist(i, p.medoid_ids[p.cluster_ids[i]]);
      dissim2 += d * d;
    }
    return dissim2;
  }

  ///
  /// Compute the total squared dissimilarity between all objects in a particular 
  /// cluster and its medoid.
  ///
  template <typename D>
  double total_squared_dissimilarity(const partition& p, D dist, medoid_id m) {
    double dissim2 = 0.0;
    for (size_t i=0; i < p.cluster_ids.size(); i++) {
      if (p.cluster_ids[i] == m) {
        double d = dist(i, p.medoid_ids[p.cluster_ids[i]]);
        dissim2 += d * d;
      }
    }
    return dissim2;
  }


} // namespace cluster

#endif // PARTITION_H
