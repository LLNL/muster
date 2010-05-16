#ifndef MPI_UTILS_H
#define MPI_UTILS_H

#include <iostream>
#include <cstdlib>
#include <stdint.h>
#include <mpi.h>

/// Overloaded function for getting corresponding MPI types for C types.
inline MPI_Datatype mpi_typeof(char)                       {return MPI_CHAR;}
inline MPI_Datatype mpi_typeof(signed short)               {return MPI_SHORT;}
inline MPI_Datatype mpi_typeof(signed int)                 {return MPI_INT;}
inline MPI_Datatype mpi_typeof(signed long)                {return MPI_LONG;}
inline MPI_Datatype mpi_typeof(unsigned char)              {return MPI_UNSIGNED_CHAR;}
inline MPI_Datatype mpi_typeof(unsigned short)             {return MPI_UNSIGNED_SHORT;}
inline MPI_Datatype mpi_typeof(unsigned)                   {return MPI_UNSIGNED;}
inline MPI_Datatype mpi_typeof(unsigned long)              {return MPI_UNSIGNED_LONG;}
inline MPI_Datatype mpi_typeof(signed long long)           {return MPI_LONG_LONG_INT;}
inline MPI_Datatype mpi_typeof(double)                     {return MPI_DOUBLE;}
inline MPI_Datatype mpi_typeof(long double)                {return MPI_LONG_DOUBLE;}
inline MPI_Datatype mpi_typeof(std::pair<int,int>)         {return MPI_2INT;}
inline MPI_Datatype mpi_typeof(std::pair<float,int>)       {return MPI_FLOAT_INT;}
inline MPI_Datatype mpi_typeof(std::pair<double,int>)      {return MPI_DOUBLE_INT;}
inline MPI_Datatype mpi_typeof(std::pair<long double,int>) {return MPI_LONG_DOUBLE_INT;}
inline MPI_Datatype mpi_typeof(std::pair<short,int>)       {return MPI_SHORT_INT;}


// helper routine for packing functions below.
// TODO: this really shouldn't call exit()
inline int mpi_packed_size(int count, MPI_Datatype type, MPI_Comm comm) {
  int size;
  if (PMPI_Pack_size(count, type, comm, &size)) {
    std::cerr << "Error getting size!" << std::endl;
    exit(1);
  }
  return size;
}

/// Handy datatype for stdlib datatypes
#define MPI_SIZE_T      (mpi_typeof(size_t()))
#define MPI_INTPTR_T    (mpi_typeof(intptr_t()))
#define MPI_UINTPTR_T   (mpi_typeof(uintptr_t()))

#endif // MPI_UTILS_H
