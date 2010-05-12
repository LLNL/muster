
#############################################################
# This is the main autoconf file for building Libra.
#############################################################
AC_INIT
AC_CONFIG_SRCDIR([kmedoids.C])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE(cluster,1.0)
AC_CONFIG_HEADERS([config.h])
AX_PREFIX_CONFIG_H([cluster-config.h])

# Programs we need
AC_PROG_CXX
AC_PROG_INSTALL

# now that we've specified deps, declare libtool
LT_INIT
AC_SUBST(LIBTOOL_DEPS)


# Test whether we're on BlueGene/L and define things in config.h accordingly
AC_LANG_PUSH([C++])   # Force using the user-supplied C++ compiler.
LX_DETECT_BLUEGENE
AC_LANG_POP

# We need the boost matrix class and boost random to build this library.
AX_BOOST_BASE

# Check for MPI.  If we don't have it, don't build parallel
# Assumes that you're using an MPI compiler (mpicxx, etc) that handles libs.
# This will NOT add MPI libraries to the build.
CC=${CXX}          # Use CXX instead of CC for tests
CPP="${CXX} -E"    # Use CXX instead of default CPP for tests
AC_CHECK_HEADER([mpi.h],
    [AC_DEFINE([HAVE_MPI], [1], [Define to 1 if you have <mpi.h>.])
     have_mpi=yes],
    [AC_MSG_NOTICE([Couldn't find mpi.h. Building without MPI.])
     have_mpi=no])
AM_CONDITIONAL([HAVE_MPI], [test "x$have_mpi" == xyes])

# Check for Linux timers.
AC_CHECK_LIB(rt, clock_gettime)
AC_CHECK_FUNCS(clock_gettime gettimeofday)

# If this is provided, use PMPI bindings instead of MPI bindings.  
AC_ARG_WITH([with-pmpi-bindings], 
  AS_HELP_STRING([--with-pmpi-bindings], 
                 [Optionally build the cluster library with name-shifted PMPI bindings instead of plain MPI bindings.  Default is to use MPI bindings.]),
  [if [[ "x$withval" = xyes ]]; then
      pmpi_bindings=no
   else
      pmpi_bindings=yes
   fi],
  [pmpi_bindings=no]
)
if [[ "x$pmpi_bindings" = xyes ]]; then
    AC_DEFINE([USE_PMPI], [1], [Define to 1 to use PMPI bindings instead of MPI bindings for communication.])
    AC_MSG_NOTICE([Building with PMPI bindings.])
else
    AC_MSG_NOTICE([Building with MPI bindings.])
fi

echo
echo "========================================================"
echo "==         Libra: final build configuration           =="
echo "========================================================"
echo "  MPI ............................................ $have_mpi"
echo "  PMPI Bindings .................................. $pmpi_bindings"
echo "========================================================"
echo

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
