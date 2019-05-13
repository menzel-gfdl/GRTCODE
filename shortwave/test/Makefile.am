if USE_SINGLE_PRECISION
PREC_FLAG = -DSINGLE_PRECISION
endif

check_PROGRAMS = test_shortwave
test_shortwave_SOURCES = test_shortwave.c
test_shortwave_LDADD = ${top_builddir}/shortwave/src/libshortwave.la \
                       ${top_builddir}/rs_utils/src/librs_utils.la \
                       -lm
test_shortwave_CPPFLAGS = $(PREC_FLAG) \
                          -I${top_builddir}/shortwave/src \
                          -I${top_builddir}/rs_utils/src
TESTS = test_shortwave
