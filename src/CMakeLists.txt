set(MUSTER_SOURCES
	partition.C
	kmedoids.C
  ../external/Timer.C
  ../external/timing.C)

if (MUSTER_HAVE_MPI)
  set(MUSTER_SOURCES ${MUSTER_SOURCES}
    par_partition.C 
    par_kmedoids.C 
    trial.C)
endif()

add_library(muster ${MUSTER_SOURCES})
if (MUSTER_HAVE_MPI)
  target_link_libraries(muster ${MPI_LIBRARIES})
endif()

include_directories(
  ${PROJECT_SOURCE_DIR}/external
  ${Boost_INCLUDE_DIR}
  ${PROJECT_BINARY_DIR}
  ${MPI_INCLUDE_DIR})

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

