AUTOMAKE_OPTIONS = foreign   # Don't require NEWS, README, INSTALL, etc.
ACLOCAL_AMFLAGS = -I m4

if INSTALL_TESTS
TEST_DIR = test
else
TEST_DIR = test
endif

SUBDIRS = src $(TEST_DIR)

# Build-time configuration of libraries.
include_HEADERS=cluster-config.h
DISTCLEANFILES = cluster-config.h _configs.sed

dist_noinst_SCRIPTS = \
	build/do_configure_darwin

EXTRA_DIST = LICENSE README
