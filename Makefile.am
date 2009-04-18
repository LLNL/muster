lib_LTLIBRARIES = libcluster.la

noinst_HEADERS = \
	CLARA.h \
	ClusterDataSet.h \
	Dissimilarity.h \
	hash.h \
	KMedoids.h \
	matrix.h \
	MersenneTwister.h \
	Point.h \
	RNGenerator.h \
	Stat.h \
	ltqnorm.h \
	random.h \
	counter.h

libcluster_la_SOURCES = \
	CLARA.C \
	ClusterDataSet.C \
	KMedoids.C \
	Point.C \
	RNGenerator.C \
	ltqnorm.C

#test program for matrix
noinst_PROGRAMS = matrixTest pointClusterTest mirkin-test

matrixTest_SOURCES = MatrixTest.C
mirkin_test_SOURCES = mirkin_test.C
pointClusterTest_SOURCES = PointClusterTest.C

LDADD = libcluster.la

