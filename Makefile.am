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
	point.h

libcluster_la_SOURCES = \
	partition.C \
	kmedoids.C \
	point.C

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
noinst_PROGRAMS = point-cluster-test mirkin-test 

if HAVE_MPI
noinst_PROGRAMS += \
	par-cluster-test \
	multi-gather-test
endif

mirkin_test_SOURCES = mirkin_test.C
point_cluster_test_SOURCES = point_cluster_test.C
par_cluster_test_SOURCES = par_cluster_test.C
multi_gather_test_SOURCES = multi_gather_test.C


LDADD = libcluster.la ../libwavelet/libwavelet.la

INCLUDES = \
	$(BOOST_CPPFLAGS) \
  -I$(top_srcdir)/libwavelet
