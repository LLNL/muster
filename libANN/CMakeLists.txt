set(LIBANN_SOURCES
	ANN.cpp
	bd_fix_rad_search.cpp
	bd_pr_search.cpp
	bd_search.cpp
	bd_tree.cpp
	brute.cpp
	kd_dump.cpp
	kd_fix_rad_search.cpp
	kd_pr_search.cpp
	kd_search.cpp
	kd_split.cpp
	kd_tree.cpp
	kd_util.cpp
	perf.cpp)

add_library(ANN ${LIBANN_SOURCES})

include_directories(
  ${PROJECT_SOURCE_DIR}/libANN
  ${PROJECT_SOURCE_DIR}/external
  ${Boost_INCLUDE_DIR}
  ${PROJECT_BINARY_DIR})



