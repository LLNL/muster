#include <cstdlib>
#include <vector>
#include <iostream>
#include "binomial.h"
using namespace std;

namespace cluster {

  binomial_embedding::binomial_embedding(int size, int root) 
    : _size(size), _root(root) { }

  /// This permutes ranks in case the root is not zero
  int binomial_embedding::relative_rank(int rank) const {
    return (rank - _root + _size) % _size;
  }

  int binomial_embedding::reverse_relative_rank(int rank) const {
    return (rank + _root) % _size;
  }

  vector<int> binomial_embedding::children(int rank) const {
    vector<int> childvec;
    get_children(rank, back_inserter(childvec));
    return childvec;
  }

  int binomial_embedding::parent(int rank) const {
    int relrank = relative_rank(rank);
    for (int mask = 0x1; mask < _size; mask <<= 1) {
      if ((mask & relrank) != 0) {
        return ((relrank & (~ mask)) + _root) % _size;
      }    
    }
    return -1;
  }

} // namespace cluster
