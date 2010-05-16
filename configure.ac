
#############################################################
# This is the main autoconf file for building Libra.
#############################################################
AC_INIT
AC_CONFIG_SRCDIR([src/kmedoids.C])
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

# If this is provided, we'll install the test library and the tests along with the
# cluster library headers.
AC_ARG_WITH([tests], 
    AS_HELP_STRING([--with-tests], 
                   [If provided, we'll install the tests and the test library too.]),
    [    if [[ "x$withval" = xyes ]]; then
             install_tests=yes
         else
             install_tests=no
         fi
     ],
    [    install_tests=no   
     ]
)
AM_CONDITIONAL([INSTALL_TESTS], [test "x$install_tests" == xyes])
if [[ "x$install_tests" = xyes ]]; then
    AC_MSG_NOTICE([Will install test cases and tests library.])
fi


echo
echo "========================================================"
echo "==         Libra: final build configuration           =="
echo "========================================================"
echo "  MPI ............................................ $have_mpi"
echo "  Install Tests .................................. $install_tests"
echo "========================================================"
echo

AC_CONFIG_FILES([Makefile \
                 src/Makefile \
                 test/Makefile])
AC_OUTPUT
