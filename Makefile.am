lib_LTLIBRARIES = libcluster.la

noinst_HEADERS = \
	kmedoids.h \
	random.h \
	counter.h \
	color.h \
	dissimilarity.h \
	partition.h \
	MersenneTwister.h \
	Stat.h \
	point.h \
	bic.h \
	trial.h \
	id_pair.h \
	gaussian.h

libcluster_la_SOURCES = \
	partition.C \
	kmedoids.C \
	point.C \
	trial.C

if HAVE_MPI
noinst_HEADERS += \
	par_partition.h \
	par_kmedoids.h \
	multi_gather.h

libcluster_la_SOURCES += \
	par_partition.C \
	par_kmedoids.C
endif

#test program for matrix
noinst_PROGRAMS = \
	point-cluster-test \
	mirkin-test \
	pathological-cluster-test \
	nondistinct-test \
	convergence-test \
	reuse-test \
	write-members-test \
	bic-test

if HAVE_MPI
noinst_PROGRAMS += \
	par-cluster-test \
	par-cluster-speed-test \
	par-bic-test \
	multi-gather-test
endif

mirkin_test_SOURCES = mirkin_test.C
point_cluster_test_SOURCES = point_cluster_test.C
par_cluster_test_SOURCES = par_cluster_test.C
par_bic_test_SOURCES = par_bic_test.C
bic_test_SOURCES = bic_test.C
par_cluster_speed_test_SOURCES = par_cluster_speed_test.C
multi_gather_test_SOURCES = multi_gather_test.C
pathological_cluster_test_SOURCES = pathological_cluster_test.C
nondistinct_test_SOURCES = nondistinct_test.C
convergence_test_SOURCES = convergence_test.C
reuse_test_SOURCES = reuse_test.C
write_members_test_SOURCES = write_members_test.C

LDADD = libcluster.la ../libwavelet/libwavelet.la ../callpath/libcallpath.la

INCLUDES = \
	$(BOOST_CPPFLAGS) \
  -I$(top_srcdir)/libwavelet \
  -I$(top_srcdir)/callpath
