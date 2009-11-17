#ifndef MULTI_GATHER_H
#define MULTI_GATHER_H

#include <mpi.h>
#include <vector>
#include <iostream>

namespace cluster {
  
  ///
  ///   T must support the following operations:
  ///     - int packed_size(MPI_Comm comm) const
  ///     - void pack(void *buf, int bufsize, int *position, MPI_Comm comm) const
  ///     - static T unpack(void *buf, int bufsize, int *position, MPI_Comm comm)
  ///
  template <class T>
  class multi_gather {
    /// internal struct for buffering sends and recvs.
    struct buffer {
      int size;       // buffer for size of Isend or Irecv
      char *buf;      // buffer for data to be sent/recv'd
      T *dest;        // destination for unpacked data.  valid only for receives.
      
      /// constructor for receive buffers
      buffer(T *_dest) 
        : size(0), buf(NULL), dest(_dest) { }

      /// constructor for send buffers
      buffer(int _size) 
        : size(_size), buf(new char[_size]), dest(NULL) { }

      /// Destructor
      ~buffer() { 
        if (buf) delete [] buf; 
      }

      /// Whether buffer has been allocated.
      bool is_allocated() { return buf; }

      /// Allocates a buffer of <size> chars, to be called after size is received.
      void allocate() { buf = new char[size]; }

      /// This is a buffer for a send if true.  It's a buffer for a recv if false.
      bool is_send()   { return !dest; }
    };


    MPI_Comm comm;                   /// Communicator on which gather takes place
    int tag;                         /// tag for communication in multi_gathers.

    std::vector<MPI_Request> reqs;   /// Oustanding requests to be completed.    
    std::vector<buffer*> buffers;     /// Send and receive buffers for packed data in gathers.
    size_t unfinished_reqs;          /// Number of still outstanding requests
    
  public:

    multi_gather(MPI_Comm _comm, int _tag=0) 
      : comm(_comm), tag(_tag), unfinished_reqs(0) { }

    /// 
    /// Starts initial send and receive requests for this gather.  Must be followed up with a call to finish().
    ///
    template <class Iterator>
    void start(const T& src_object, Iterator first_src, Iterator last_src, std::vector<T>& dest, int root) {
      int size, rank;
      MPI_Comm_size(comm, &size);
      MPI_Comm_rank(comm, &rank);
      
      if (rank != root) {
        // stop if this rank isn't a member of the gather.
        if (find(first_src, last_src, rank) == last_src) return;

        // we're not on the root, so we need to send to the root.
        int packed_size = src_object.packed_size(comm);      // first, send size
        MPI_Send(&packed_size, 1, MPI_INT, root, tag, comm);
        
        buffers.push_back(new buffer(packed_size));              // buffer data during send
        reqs.push_back(MPI_REQUEST_NULL);       

        // send packed data along to destination.
        int pos = 0;
        src_object.pack(buffers.back()->buf, packed_size, &pos, comm);
        MPI_Isend(buffers.back()->buf, packed_size, MPI_PACKED, root, tag, comm, &reqs.back());
        unfinished_reqs++;

      } else {        // rank is root; do receives instead
        // reserve space in the destination vector for all received T's
        dest.resize(distance(first_src, last_src));

        // initiate all the receives for sizes
        Iterator src = first_src;
        size_t i = 0;
        while (src != last_src) {
          if (*src == root) {
            dest[i] = src_object;   // on root, copy local object and leave request as null.
          } else {
            // make some buffer space for the receive, record its eventual destination
            buffers.push_back(new buffer(&dest[i]));
            reqs.push_back(MPI_REQUEST_NULL);

            MPI_Irecv(&buffers.back()->size, 1, MPI_INT, *src, tag, comm, &reqs.back());
            unfinished_reqs++;
          }

          src++;
          i++;
        }
      }
    }
    

    void finish() {
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);

      while (unfinished_reqs) {
        int outcount;
        std::vector<int> indices(reqs.size());
        std::vector<MPI_Status> status(reqs.size());

        MPI_Waitsome(reqs.size(), &reqs[0], &outcount, &indices[0], &status[0]);
        for (int o=0; o < outcount; o++) {
          const int r = indices[o];   // index of received object.

          if (buffers[r]->is_send()) {
            // buffers[r] is a send, and we just need to note that it's done.
            delete buffers[r];
            unfinished_reqs--;

          } else if (!buffers[r]->is_allocated()) {
            // buffers[r] is a recv and we just received packed size.  Allocate space and recv data.
            int src = status[o].MPI_SOURCE;
            buffers[r]->allocate();
            MPI_Irecv(buffers[r]->buf, buffers[r]->size, MPI_PACKED, src, tag, comm, &reqs[r]);

          } else {
            // buffers[r] is a recv and we just received full packed data, so unpack it.
            int pos = 0;
            *buffers[r]->dest = T::unpack(buffers[r]->buf, buffers[r]->size, &pos, comm);
            delete buffers[r];
            unfinished_reqs--;
          }
        }
      }

      // clear these out before the next call to start()
      buffers.clear();
      reqs.clear();
    }
    
  }; // class multi_gather
  
} // namespace cluster  
    

#endif // MULTI_GATHER_H

