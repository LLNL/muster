#ifndef PARTITION_H
#define PARTITION_H

#include <cstddef>
#include <vector>
#include <set>
#include <ostream>

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
  };

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

} // namespace cluster

#endif // PARTITION_H
