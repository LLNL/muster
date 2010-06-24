lib_LTLIBRARIES = libmuster.la

include_HEADERS = \
	kmedoids.h \
	random.h \
	counter.h \
	dissimilarity.h \
	partition.h \
	bic.h \
	mpi_bindings.h

EXTERNAL_SOURCES = \
	../external/Timer.C \
	../external/timing.C

libmuster_la_SOURCES = \
	partition.C \
	kmedoids.C \
	$(EXTERNAL_SOURCES)

if HAVE_MPI
include_HEADERS += \
	par_partition.h \
	par_kmedoids.h \
	multi_gather.h \
	trial.h \
	id_pair.h

libmuster_la_SOURCES += \
	par_partition.C \
	par_kmedoids.C \
	trial.C
endif

INCLUDES = \
	$(BOOST_CPPFLAGS) \
	-I$(top_srcdir)/external

