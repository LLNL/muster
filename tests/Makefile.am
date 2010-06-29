noinst_HEADERS = \
	point.h \
	gaussian.h \
	color.h \
	spherical_clustering_generator.h

EXTRA_DIST = $(noinst_HEADERS)

noinst_LTLIBRARIES = \
	libtestsupport.la

libtestsupport_la_SOURCES = \
	point.cpp

#test program for matrix
bin_PROGRAMS = \
	point-cluster-test \
	mirkin-test \
	pathological-cluster-test \
	nondistinct-test \
	convergence-test \
	write-members-test \
	bic-test \
	reuse-test

if HAVE_MPI
bin_PROGRAMS += \
	par-cluster-test \
	par-cluster-speed-test \
	par-bic-test \
	multi-gather-test
endif

mirkin_test_SOURCES = mirkin_test.cpp
point_cluster_test_SOURCES = point_cluster_test.cpp
bic_test_SOURCES = bic_test.cpp
pathological_cluster_test_SOURCES = pathological_cluster_test.cpp
nondistinct_test_SOURCES = nondistinct_test.cpp
convergence_test_SOURCES = convergence_test.cpp
reuse_test_SOURCES = reuse_test.cpp
write_members_test_SOURCES = write_members_test.cpp

par_cluster_speed_test_SOURCES = par_cluster_speed_test.cpp
par_cluster_test_SOURCES = par_cluster_test.cpp
par_bic_test_SOURCES = par_bic_test.cpp
multi_gather_test_SOURCES = multi_gather_test.cpp

LDADD = \
	libtestsupport.la \
	../src/libmuster.la

INCLUDES = \
	-I$(top_srcdir)/external \
	-I$(top_srcdir)/src \
	$(BOOST_CPPFLAGS)

