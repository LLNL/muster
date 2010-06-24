noinst_HEADERS = \
	point.h \
	gaussian.h \
	color.h \
	spherical_clustering_generator.h

EXTRA_DIST = $(noinst_HEADERS)

noinst_LTLIBRARIES = \
	libtestsupport.la

libtestsupport_la_SOURCES = \
	point.C

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

mirkin_test_SOURCES = mirkin_test.C
point_cluster_test_SOURCES = point_cluster_test.C
bic_test_SOURCES = bic_test.C
pathological_cluster_test_SOURCES = pathological_cluster_test.C
nondistinct_test_SOURCES = nondistinct_test.C
convergence_test_SOURCES = convergence_test.C
reuse_test_SOURCES = reuse_test.C
write_members_test_SOURCES = write_members_test.C

par_cluster_speed_test_SOURCES = par_cluster_speed_test.C
par_cluster_test_SOURCES = par_cluster_test.C
par_bic_test_SOURCES = par_bic_test.C
multi_gather_test_SOURCES = multi_gather_test.C

LDADD = \
	libtestsupport.la \
	../src/libmuster.la

INCLUDES = \
	-I$(top_srcdir)/external \
	-I$(top_srcdir)/src \
	$(BOOST_CPPFLAGS)

