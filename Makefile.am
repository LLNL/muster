AUTOMAKE_OPTIONS = foreign   # Don't require NEWS, README, INSTALL, etc.
ACLOCAL_AMFLAGS = -I m4

if INSTALL_TESTS
TEST_DIR = tests
else
TEST_DIR = tests
endif

SUBDIRS = src $(TEST_DIR) external

# Build-time configuration of libraries.
include_HEADERS=muster-config.h
DISTCLEANFILES = muster-config.h _configs.sed

dist_noinst_SCRIPTS = \
	build/do_configure_darwin

EXTRA_DIST = LICENSE README
