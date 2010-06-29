#############################################################
# This is the main autoconf file for building Muster
#############################################################
AC_INIT
AC_CONFIG_SRCDIR([src/kmedoids.cpp])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE(muster,1.0)
AC_CONFIG_HEADERS([config.h])
AX_PREFIX_CONFIG_H([muster-config.h])

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
AC_LANG_PUSH([C++])   # Force using the user-supplied C++ compiler.
LX_MPI
AM_CONDITIONAL([HAVE_MPI], [test "x$have_mpi" == xyes])
if [[ x$have_mpi != xno ]]; then
    AC_DEFINE([HAVE_MPI], [1], [Define to 1 if you have <mpi.h>.])
else
    AC_MSG_NOTICE([Couldn't find MPI library or header. Building without MPI.])
fi
AC_LANG_POP

# Check for Linux timers.
AC_CHECK_LIB(rt, clock_gettime)
AC_CHECK_FUNCS(clock_gettime gettimeofday)

# If this is provided, we'll install the test library and the tests along with the
# muster library headers.
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


# If this is provided, we use PMPI bindings instead of MPI bindings for communication.
AC_ARG_WITH([pmpi-bindings], 
    AS_HELP_STRING([--with-pmpi-bindings], 
                   [If provided, build with PMPI bindings instead of raw MPI bindings.  Useful for including this library in tools.]),
    [    if [[ "x$withval" = xyes ]]; then
             pmpi_bindings=yes
         else
             pmpi_bindings=no
         fi
     ],
    [    pmpi_bindings=no
     ]
)
if [[ "x$pmpi_bindings" = xyes ]]; then
    AC_DEFINE([USE_PMPI], [1], [Define to 1 to use PMPI bindings instead of raw MPI bindings.])
fi

echo
echo "========================================================"
echo "==        muster: final build configuration           =="
echo "========================================================"
echo "  MPI ............................................ $have_mpi"
echo "  PMPI Bindings .................................. $pmpi_bindings"
echo "  Install Tests .................................. $install_tests"
echo "========================================================"
echo

AC_CONFIG_FILES([Makefile \
                 src/Makefile \
                 tests/Makefile \
                 external/Makefile])
AC_OUTPUT
