###############################################################
# Muster: Massively Scalable Clustering
# by Todd Gamblin
# tgamblin@llnl.gov
#
# This is the top-level CMake build file for Muster.
###############################################################
project(muster)
cmake_minimum_required(VERSION 2.8)

# Version information comes from here.  Incrememnt these when new stuff happens.
set(MUSTER_MAJOR_VERSION 1)
set(MUSTER_MINOR_VERSION 0)
set(MUSTER_PATCH_VERSION 0)
set(MUSTER_VERSION "${MUSTER_MAJOR_VERSION}.${MUSTER_MINOR_VERSION}.${MUSTER_PATCH_VERSION}")

# Build shared libraries by default.
get_property(DEFAULT_SHARED GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS)
option(BUILD_SHARED_LIBS "Build shared libraries?" ${DEFAULT_SHARED})

# Whether to build for regular MPI apps or with PMPI bindings (for tools).
option(MUSTER_USE_PMPI "Build with PMPI bindings?" FALSE)

# Tests are off by default.  Set this to ON to build tests.
option(INSTALL_TESTS "Do you want to build the test programs?" FALSE)

# Muster requires boost for matrices and random number generators.
find_package(Boost REQUIRED)

# Find the MPI library and set some definitions
find_package(MPI)
if (MPI_FOUND)
  set(MUSTER_HAVE_MPI TRUE)
endif()

# Check for various timing functions, so we can support highest-resolution timers available.
include(CheckFunctionExists)

set(CMAKE_REQUIRED_INCLUDES  time.h sys/time.h)
set(CMAKE_REQUIRED_LIBRARIES rt)
check_function_exists(clock_gettime MUSTER_HAVE_CLOCK_GETTIME)

unset(CMAKE_REQUIRED_LIBRARIES)
check_function_exists(gettimeofday MUSTER_HAVE_GETTIMEOFDAY)

set(CMAKE_REQUIRED_INCLUDES mach/mach.h mach/mach_time.h)
unset(CMAKE_REQUIRED_LIBRARIES)
check_function_exists(mach_absolute_time MUSTER_HAVE_MACH_TIME)

# Configure a header file with all the variables we found.
configure_file(
  ${PROJECT_SOURCE_DIR}/muster-config.h.in
  ${PROJECT_BINARY_DIR}/muster-config.h)
install(FILES ${PROJECT_BINARY_DIR}/muster-config.h DESTINATION include)

# Configure subdirectories here, after everything else, so that they 
# inherit variables set in the global scope.
add_subdirectory(src)
add_subdirectory(tests)
