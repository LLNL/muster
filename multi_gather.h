#ifndef MULTI_GATHER_H
#define MULTI_GATHER_H

#include <mpi.h>
#include <vector>
#include <iostream>
#include "mpi_utils.h"
#ifdef AIX
#include <algorithm>
#endif

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
      int size;             /// buffer for size of Isend or Irecv
      char *buf;            /// buffer for data to be sent/recv'd
      std::vector<T> *dest; /// vector to push unpacked data onto

      /// constructor for receive buffers
      buffer(std::vector<T>& _dest) 
        : size(0), buf(NULL), dest(&_dest) { }

      /// constructor for send buffers
      buffer(int _size) 
        : size(_size), buf(new char[_size]), dest(NULL) { }

      /// Destructor
      ~buffer() { 
        if (buf) delete [] buf; 
      }

      /// Turn a send buffer into a receive buffer (for local "sends")
      void set_destination(std::vector<T>& _dest) { dest = &_dest; }

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
    template <class ObjIterator, class RankIterator>
    void start(ObjIterator  begin_obj, ObjIterator  end_obj, 
               RankIterator begin_src, RankIterator end_src, 
               std::vector<T>& dest, int root) 
    {
      int size, rank;
      MPI_Comm_size(comm, &size);
      MPI_Comm_rank(comm, &rank);

      // stop if this rank isn't a member of the gather.
      if (rank != root && find(begin_src, end_src, rank) == end_src) return;

      // determine size and pack local data.
      int packed_size = mpi_packed_size(1, MPI_SIZE_T, comm);  // num objects
      for (ObjIterator o=begin_obj; o != end_obj; o++) {       // size of each object
        packed_size += o->packed_size(comm);
      }

      if (rank != root) {
        MPI_Send(&packed_size, 1, MPI_INT, root, tag, comm);
      }

      // pack up local data into a buffer
      buffer *send_buffer = new buffer(packed_size);
      int pos = 0;
      size_t num_objects = distance(begin_obj, end_obj);
      MPI_Pack(&num_objects, 1, MPI_SIZE_T, send_buffer->buf, send_buffer->size, &pos, comm);
      for (ObjIterator o=begin_obj; o != end_obj; o++) {
        o->pack(send_buffer->buf, packed_size, &pos, comm);
      }

      if (rank != root) {
        // send packed data along to destination.
        buffers.push_back(send_buffer);     // buffer data during send
        reqs.push_back(MPI_REQUEST_NULL);
        MPI_Isend(send_buffer->buf, packed_size, MPI_PACKED, root, tag, comm, &reqs.back());
        unfinished_reqs++;

      } else {        // rank is root; do receives instead
        // initiate all the receives for sizes
        for (RankIterator src=begin_src; src != end_src; src++) {
          if (*src == root) {
            // for the root, just insert the local packed buffer onto the array of buffers.
            // don't increment unfinished reqs here b/c local buffer is already "done."
            send_buffer->set_destination(dest);
            buffers.push_back(send_buffer);
            reqs.push_back(MPI_REQUEST_NULL);

          } else {
            // make some buffer space for the receive, record its eventual destination
            buffers.push_back(new buffer(dest));
            reqs.push_back(MPI_REQUEST_NULL);
            MPI_Irecv(&buffers.back()->size, 1, MPI_INT, *src, tag, comm, &reqs.back());
            unfinished_reqs++;
          }
        }
      }
    }

    ///
    /// Starts a gather with one object instead of a range of objects.
    ///
    template <class RankIterator>
    void start(const T& obj, RankIterator begin_src, RankIterator end_src, std::vector<T>& dest, int root) {
      start(&obj, (&obj) + 1, begin_src, end_src, dest, root);
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

          if (!buffers[r]->is_send() && !buffers[r]->is_allocated()) {
            // buffers[r] is a recv and we just received packed size.  Allocate space and recv data.
            int src = status[o].MPI_SOURCE;
            buffers[r]->allocate();
            MPI_Irecv(buffers[r]->buf, buffers[r]->size, MPI_PACKED, src, tag, comm, &reqs[r]);

          } else {
            // buffers[r] is a send, or it's a receive and we just received full packed data.
            // in either case, the buffer is done, so decrement the number of unfinished reqs.
            unfinished_reqs--;
          }
        }
      }

      // Unpack all the received buffers into their destination vectors.  This preserves order
      // as unpacked data are only pushed onto the backs of destination vectors *after* everything
      // is received.  Buffers are still received in any order above, though.
      for (size_t i=0; i < buffers.size(); i++) {
        if (!buffers[i]->is_send()) {
          int pos = 0;
          size_t num_objects;

          MPI_Unpack(buffers[i]->buf, buffers[i]->size, &pos, &num_objects, 1, MPI_SIZE_T, comm);
          for (size_t o=0; o < num_objects; o++) {
            buffers[i]->dest->push_back(T::unpack(buffers[i]->buf, buffers[i]->size, &pos, comm));
          }
        }
        delete buffers[i];
      }
            
      // clear these out before the next call to start()
      buffers.clear();
      reqs.clear();
    }
    
  }; // class multi_gather
  
} // namespace cluster  
    

#endif // MULTI_GATHER_H

