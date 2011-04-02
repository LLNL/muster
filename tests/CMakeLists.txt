# Support library for tests that don't need to be in the main muster library 
# (or they haven't been moved there yet)
add_library(test-support STATIC
  point.cpp
  point_set.cpp)

# Simple function to add an executable and link it to the test libraries.
function(add_test test_name src_name)
  add_executable(${test_name} ${src_name})
  target_link_libraries(${test_name} test-support muster)
  if (INSTALL_TESTS)
    install(TARGETS ${test_name} DESTINATION bin)
  endif()
endfunction()

# Function to add MPI tests, but only if MPI is enabled.
function(add_mpi_test test_name src_name)
  if (MUSTER_HAVE_MPI)
    add_test(${test_name} ${src_name})
    target_link_libraries(${test_name} ${MPI_LIBRARIES})
  endif()
endfunction()

add_test(point-cluster-test point_cluster_test.cpp)
add_test(mirkin-test mirkin_test.cpp)
add_test(pathological-cluster-test pathological_cluster_test.cpp)
add_test(nondistinct-test nondistinct_test.cpp)
add_test(convergence-test convergence_test.cpp)
add_test(write-members-test write_members_test.cpp)
add_test(bic-test bic_test.cpp)
add_test(reuse-test reuse_test.cpp)
add_test(random-test random_test.cpp)

add_mpi_test(par-cluster-test par_cluster_test.cpp)
add_mpi_test(par-cluster-speed-test par_cluster_speed_test.cpp)
add_mpi_test(par-bic-test par_bic_test.cpp)
add_mpi_test(multi-gather-test multi_gather_test.cpp)
add_mpi_test(gather-test gather_test.cpp)

include_directories(
  ${PROJECT_SOURCE_DIR}/external
  ${PROJECT_SOURCE_DIR}/src
  ${PROJECT_BINARY_DIR}
  ${Boost_INCLUDE_DIR}
  ${MPI_INCLUDE_PATH})
