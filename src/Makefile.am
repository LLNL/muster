lib_LTLIBRARIES = libmuster.la

include_HEADERS = \
	kmedoids.h \
	density_based.h \
	random.h \
	counter.h \
	dissimilarity.h \
	partition.h \
	bic.h \
	mpi_bindings.h

EXTERNAL_SOURCES = \
	../external/Timer.cpp \
	../external/timing.cpp

libmuster_la_SOURCES = \
	partition.cpp \
	kmedoids.cpp \
	density_based.cpp \
	$(EXTERNAL_SOURCES)

if HAVE_MPI
include_HEADERS += \
	par_partition.h \
	par_kmedoids.h \
	multi_gather.h \
	trial.h \
	id_pair.h

libmuster_la_SOURCES += \
	par_partition.cpp \
	par_kmedoids.cpp \
	trial.cpp
endif

INCLUDES = \
	$(BOOST_CPPFLAGS) \
	-I$(top_srcdir)/external

