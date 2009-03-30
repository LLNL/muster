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
	ltqnorm.h

libcluster_la_SOURCES = \
	CLARA.C \
	ClusterDataSet.C \
	KMedoids.C \
	Point.C \
	RNGenerator.C \
	ltqnorm.C

#test program for matrix
noinst_PROGRAMS = matrixTest pointClusterTest

matrixTest_SOURCES = MatrixTest.C
pointClusterTest_SOURCES = PointClusterTest.C

LDADD = libcluster.la

