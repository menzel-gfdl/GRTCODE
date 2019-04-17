if USE_SINGLE_PRECISION
PREC_FLAG = -DSINGLE_PRECISION
endif

lib_LTLIBRARIES = libshortwave.la
libshortwave_la_SOURCES = rayleigh.c shortwave.c solar_flux.c
libshortwave_la_LIBADD = ${top_builddir}/rs_utils/src/librs_utils.la
libshortwave_la_CPPFLAGS = $(PREC_FLAG) -I${top_builddir}/rs_utils/src
