set(MUSTER_SOURCES
	partition.cpp
	kmedoids.cpp
	density_based.cpp
  ../external/Timer.cpp
  ../external/timing.cpp)

if (MUSTER_HAVE_MPI)
  set(MUSTER_SOURCES ${MUSTER_SOURCES}
    par_partition.cpp 
    par_kmedoids.cpp 
    trial.cpp)
endif()

add_library(muster ${MUSTER_SOURCES})
if (MUSTER_HAVE_MPI)
  target_link_libraries(muster ${MPI_LIBRARIES})
endif()

include_directories(
  ${PROJECT_SOURCE_DIR}/external
  ${Boost_INCLUDE_DIR}
  ${PROJECT_BINARY_DIR})

if (MUSTER_HAVE_MPI)
	include_directories(${MPI_INCLUDE_PATH})
endif()

# include_HEADERS = \
# 	kmedoids.h \
# 	random.h \
# 	counter.h \
# 	dissimilarity.h \
# 	partition.h \
# 	bic.h \
# 	mpi_bindings.h


# if HAVE_MPI
# include_HEADERS += \
# 	par_partition.h \
# 	par_kmedoids.h \
# 	multi_gather.h \
# 	trial.h \
# 	id_pair.h
# endif

# INCLUDES = \
# 	$(BOOST_CPPFLAGS) \
# 	-I$(top_srcdir)/external

