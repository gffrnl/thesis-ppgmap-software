EXTRA_DIST = third-party vendor

AM_CPPFLAGS = -I$(top_srcdir)/include     \
              -I$(top_srcdir)/third-party \
              -I$(top_srcdir)/vendor

lib_LTLIBRARIES = libfrlap.la
libfrlap_la_SOURCES = src/frlap.c \
                      src/frlap1gd.c \
                      src/frlap1gd/dcgtor_spectral.c \
                      src/frlap1gd/dcgtor_huang_oberman.c \
                      src/frlap1gd/dcgtor_huang_oberman_linear.c \
                      src/frlap1gd/dcgtor_huang_oberman_quadratic.c \
                      src/frlap1gd/dcgtor_gorenflo_mainardi.c \
                      src/frlap1gd/dcgtor_centered_3_point_periodized.c \
                      src/frlap1gd/dcgtor_centered_5_point_periodized.c


pkginclude_frlap1gddir = ${pkgincludedir}/frlap1gd

pkginclude_HEADERS = include/frlap/frlap.h include/frlap/frlap1gd.h
pkginclude_frlap1gd_HEADERS = include/frlap/frlap1gd/dcgtors.h


uninstall-hook:
	-rmdir ${pkginclude_frlap1gddir}
	-rmdir ${pkgincludedir}

purge-emacs:
	rm -rf *~ include/*/*~ include/*/*/*~ src/*/*~

CLEANFILES = src/*.o
