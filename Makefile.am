lib_LTLIBRARIES = libcluster.la

noinst_HEADERS = \
	KMedoids.h \
	MersenneTwister.h \
	RNGenerator.h \
	Stat.h \
	ltqnorm.h \
	random.h \
	counter.h \
	dissimilarity.h \
	partition.h

libcluster_la_SOURCES = \
	partition.C \
	KMedoids.C \
	RNGenerator.C \
	ltqnorm.C

#test program for matrix
noinst_PROGRAMS = pointClusterTest mirkin-test

mirkin_test_SOURCES = mirkin_test.C
pointClusterTest_SOURCES = PointClusterTest.C

LDADD = libcluster.la

INCLUDES=$(BOOST_CPPFLAGS)
