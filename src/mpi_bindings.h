//////////////////////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2010, Lawrence Livermore National Security, LLC.  
// Produced at the Lawrence Livermore National Laboratory  
// LLNL-CODE-433662
// All rights reserved.  
//
// This file is part of Muster. For details, see http://github.com/tgamblin/muster. 
// Please also read the LICENSE file for further information.
//
// Redistribution and use in source and binary forms, with or without modification, are
// permitted provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright notice, this list of
//    conditions and the disclaimer below.
//  * Redistributions in binary form must reproduce the above copyright notice, this list of
//    conditions and the disclaimer (as noted below) in the documentation and/or other materials
//    provided with the distribution.
//  * Neither the name of the LLNS/LLNL nor the names of its contributors may be used to endorse
//    or promote products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
// OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
// LAWRENCE LIVERMORE NATIONAL SECURITY, LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
// ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//////////////////////////////////////////////////////////////////////////////////////////////////

///
/// @file mpi_bindings.h
/// @author Todd Gamblin tgamblin@llnl.gov
/// @brief <code>\#defines</code> for switching between MPI and PMPI bindings.
/// 
/// User of the API can <code>\#define</code> <code>MUSTER_USE_PMPI</code> to use the 
/// PMPI bindings instead of the plain MPI bindings.  Useful
/// for including this algorithm in tools.
/// 
/// This file should contain <code>\#defines</code> for all MPI calls used in 
/// the cluster library, and needs to be kept current.
///
#ifndef MUSTER_MPI_BINDINGS_H
#define MUSTER_MPI_BINDINGS_H

#include "muster-config.h"

// External header for MPI type information
#include "mpi_utils.h"

#ifdef MUSTER_USE_PMPI

#define CMPI_Allreduce   PMPI_Allreduce
#define CMPI_Bcast       PMPI_Bcast
#define CMPI_Comm_rank   PMPI_Comm_rank
#define CMPI_Comm_size   PMPI_Comm_size
#define CMPI_Gather      PMPI_Gather
#define CMPI_Scatter     PMPI_Scatter
#define CMPI_Recv        PMPI_Recv
#define CMPI_Send        PMPI_Send
#define CMPI_Irecv       PMPI_Irecv
#define CMPI_Isend       PMPI_Isend
#define CMPI_Pack        PMPI_Pack
#define CMPI_Pack_size   PMPI_Pack_size
#define CMPI_Reduce      PMPI_Reduce
#define CMPI_Unpack      PMPI_Unpack
#define CMPI_Waitsome    PMPI_Waitsome
#define CMPI_Comm_free   PMPI_Comm_free
#define CMPI_Comm_group  PMPI_Comm_group
#define CMPI_Comm_create PMPI_Comm_create
#define CMPI_Group_incl  PMPI_Group_incl
#define CMPI_Group_free  PMPI_Group_free

#define cmpi_packed_size pmpi_packed_size

#else  // MUSTER_USE_PMPI

#define CMPI_Allreduce   MPI_Allreduce
#define CMPI_Bcast       MPI_Bcast
#define CMPI_Comm_rank   MPI_Comm_rank
#define CMPI_Comm_size   MPI_Comm_size
#define CMPI_Gather      MPI_Gather
#define CMPI_Scatter     MPI_Scatter
#define CMPI_Recv        MPI_Recv
#define CMPI_Send        MPI_Send
#define CMPI_Irecv       MPI_Irecv
#define CMPI_Isend       MPI_Isend
#define CMPI_Pack        MPI_Pack
#define CMPI_Pack_size   MPI_Pack_size
#define CMPI_Reduce      MPI_Reduce
#define CMPI_Unpack      MPI_Unpack
#define CMPI_Waitsome    MPI_Waitsome
#define CMPI_Comm_free   MPI_Comm_free
#define CMPI_Comm_group  MPI_Comm_group
#define CMPI_Comm_create MPI_Comm_create
#define CMPI_Group_incl  MPI_Group_incl
#define CMPI_Group_free  MPI_Group_free

#define cmpi_packed_size mpi_packed_size

#endif // MUSTER_USE_PMPI



#endif // MUSTER_MPI_BINDINGS_H
