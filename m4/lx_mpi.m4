#
# LX_MPI ([default-location])
#  ------------------------------------------------------------------------
# Tests for presence of MPI. Doesn't support a custom MPI location; 
# assumes you're compiling with an MPI compiler wrapper.
#
AC_DEFUN([LX_MPI], [
    AC_TRY_LINK([#include <mpi.h>],
                [int foo; MPI_Initialized(&foo);],
                [have_mpi=yes],
                [have_mpi=no])
])
