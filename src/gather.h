#ifndef MUSTER_GATHER_H
#define MUSTER_GATHER_H

#include <mpi.h>
#include <algorithm>
#include "mpi_bindings.h"
#include "mpi_utils.h"

namespace muster {

  ///
  /// Binomial gather of char buffers into a single agglomerated clump of buffers
  /// TODO: clean this up and make it do packing automatically.
  ///
  void gather(const std::vector<char>& src, std::vector<char>& dest,
              MPI_Comm comm, int root = 0) {
    int rank, size;
    CMPI_Comm_rank(comm, &rank);
    CMPI_Comm_size(comm, &size);

    int relrank = (rank - root + size) % size;

    // start with one int for the total number of buffers packed.
    vector<int> sizes;   // sizes of buffers to receive
    vector<int> counts;  // number of sub-buffers in each
    for (int mask = 0x1; mask < size; mask <<= 1) {
      // Receive sizes from all children
      if ((mask & relrank) == 0) {
        int child = (relrank | mask);
        if (child < size) {
          child = (child + root) % size;

          // add in the size of the received buffer
          int buf[2];
          CMPI_Recv(buf, 2, MPI_INT, child, 0, comm, MPI_STATUS_IGNORE);
          sizes.push_back(buf[0]);
          counts.push_back(buf[1]);
        }
      }
    }

    // create a buffer to receive into, then to send to parent
    std::vector<char> sendbuf(sum(sizes.begin(), sizes.end()));
    
    int count = sum(counts)

    int pos = 0;
    for (int mask = 0x1; mask < size; mask <<= 1, count++) {
      if ((mask & relrank) == 0) {
        int child = (relrank | mask);
        if (child < size) {
          child = (child + root) % size;
          
        }
      } else {
        // I've received all that I'm going to.  Send my result to my parent
        int parent = ((relrank & (~ mask)) + root) % size;
        
        break;
      }
    }
  }


  /*
   * ----------------------------------
   * Binomial Tree Reduction from MPICH
   * ----------------------------------
   *
   * This algorithm is adapted from MPICH. It performs a binomial-tree reduction
   * assuming that the operations are commutative. The root process can be
   * specified in the reduce() function.
   */
  template <class T>
  void some_to_many(vector<T>int root, int numProcesses, int processRank, T *reducedObject) {
		int mask = 0x1, relrank, source, destination;
		relrank = (processRank - root + numProcesses) % numProcesses;
		while (mask < numProcesses) {
			// Receive
			if ((mask & relrank) == 0) {
				source = (relrank | mask);
				if (source < numProcesses) {
					source = (source + root) % numProcesses;
					//cout << "Proc " << processRank << " recv from " << source << endl;
					reducedObject->receive(source);
				}
			} else {
				// I've received all that I'm going to.  Send my result to my parent
				destination = ((relrank & (~ mask)) + root) % numProcesses;
				//cout << "Proc " << processRank << " send to " << destination << endl;
				reducedObject->send(destination);
				break;
			}
			mask <<= 1;
		}
	}

} // namespace muster

#endif // MUSTER_GATHER_H
