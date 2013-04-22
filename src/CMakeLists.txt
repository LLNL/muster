set(MUSTER_SOURCES
	partition.cpp
	kmedoids.cpp
  binomial.cpp
  ../external/Timer.cpp
  ../external/timing.cpp)

set(MUSTER_HEADERS
 	kmedoids.h
 	random.h
 	counter.h
 	dissimilarity.h
 	partition.h
  binomial.h
  gather.h
  packable_vector.h
 	bic.h)

if (MUSTER_HAVE_MPI)
  list(APPEND MUSTER_SOURCES
    par_partition.cpp
    par_kmedoids.cpp
    trial.cpp)

  list(APPEND MUSTER_HEADERS
 	  par_partition.h
 	  par_kmedoids.h
 	  multi_gather.h
 	  trial.h
 	  id_pair.h
    mpi_bindings.h
    ../external/Timer.h
    ../external/timing.h
    ../external/stl_utils.h
    ../external/mpi_utils.h)
endif()

add_library(muster ${MUSTER_SOURCES})
if (MUSTER_HAVE_MPI)
  target_link_libraries(muster ${MPI_LIBRARIES})
endif()

if (MUSTER_HAVE_CLOCK_GETTIME)
  target_link_libraries(muster rt)
endif()

include_directories(
  ${PROJECT_SOURCE_DIR}/external
  ${Boost_INCLUDE_DIR}
  ${PROJECT_BINARY_DIR}
  ${MPI_INCLUDE_PATH})

install(FILES ${MUSTER_HEADERS} DESTINATION include)
install(TARGETS muster DESTINATION lib)
