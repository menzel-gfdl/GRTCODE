if USE_SINGLE_PRECISION
PREC_FLAG = -DSINGLE_PRECISION
endif

check_PROGRAMS = test_longwave
test_longwave_SOURCES = test_longwave.c
test_longwave_LDADD = ${top_builddir}/longwave/src/liblongwave.la \
                      ${top_builddir}/rs_utils/src/librs_utils.la \
                      -lm
test_longwave_CPPFLAGS = $(PREC_FLAG) \
                         -I${top_builddir}/longwave/src \
                         -I${top_builddir}/rs_utils/src
TESTS = test_longwave
