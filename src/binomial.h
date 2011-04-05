#ifndef MUSTER_BINOMIAL_H
#define MUSTER_BINOMIAL_H

#include <cstdlib>
#include <vector>

namespace cluster {
  
  class binomial_embedding {
  private:
    const int _size;
    const int _root;
    
  public:
    /// Construct a binomial rank embedding with size nodes, rooted at root.
    binomial_embedding(int size, int root = 0);
    
    int relative_rank(int rank) const;          ///< This permutes ranks in case the root is not zero
    int reverse_relative_rank(int rank) const;  ///< Reverse rank permutation
    std::vector<int> children(int rank) const;  ///< Same as get_children, but returns vector.
    int parent(int rank) const;                 ///< Get the parent of a particular rank.
    
    int size() const { return _size; }
    int root() const { return _root; }
    
    /// This allows you to putting children into any structure that supports output iterators
    template <class OutputIterator>
    void get_children(int rank, OutputIterator o) const {
      int relrank = relative_rank(rank);
      for (int mask = 0x1; mask < _size; mask <<= 1) {
        if ((mask & relrank) != 0) {
          break;
        }
        
        int child = (relrank | mask);
        if (child < _size) {
          *o++ = (child + _root) % _size;
        }
      }
    }
  };

} // namespace cluster

#endif // MUSTER_BINOMIAL_H

