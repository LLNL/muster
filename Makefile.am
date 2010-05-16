AUTOMAKE_OPTIONS = foreign   # Don't require NEWS, README, INSTALL, etc.
ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src test

# Build-time configuration of libraries.
include_HEADERS=cluster-config.h
DISTCLEANFILES = cluster-config.h _configs.sed

dist_noinst_SCRIPTS = \
	build/do_configure_atlas \
	build/do_configure_bgl \
	build/do_configure_bgl_fen \
	build/do_configure_bgp \
	build/do_configure_bgp_fen \
	build/do_configure_darwin

EXTRA_DIST = LICENSE README
