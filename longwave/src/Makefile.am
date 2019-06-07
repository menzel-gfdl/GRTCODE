if USE_SINGLE_PRECISION
PREC_FLAG = -DSINGLE_PRECISION
endif

lib_LTLIBRARIES = liblongwave.la
liblongwave_la_SOURCES = longwave.c
liblongwave_la_LIBADD = ${top_builddir}/rs_utils/src/librs_utils.la
liblongwave_la_CPPFLAGS = $(PREC_FLAG) -I${top_builddir}/rs_utils/src
