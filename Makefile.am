lib_LTLIBRARIES = libcluster.la

noinst_HEADERS = \
	kmedoids.h \
	ltqnorm.h \
	random.h \
	counter.h \
	dissimilarity.h \
	partition.h\
	MersenneTwister.h \
	RNGenerator.h \
	Stat.h

libcluster_la_SOURCES = \
	partition.C \
	kmedoids.C \
	RNGenerator.C \
	ltqnorm.C

#test program for matrix
noinst_PROGRAMS = point-cluster-test mirkin-test

mirkin_test_SOURCES = mirkin_test.C
point_cluster_test_SOURCES = point_cluster_test.C

LDADD = libcluster.la

INCLUDES=$(BOOST_CPPFLAGS)
