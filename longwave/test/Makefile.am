AM_CPPFLAGS = -I${top_builddir}/longwave/src \
              -I${top_builddir}/rs_utils/src
if USE_SINGLE_PRECISION
AM_CPPFLAGS += -DSINGLE_PRECISION
endif

check_PROGRAMS = test_longwave
test_longwave_SOURCES = test_longwave.c
test_longwave_LDADD = ${top_builddir}/longwave/src/liblongwave.a \
                      ${top_builddir}/rs_utils/src/librs_utils.a \
                      -lm
if ENABLE_CUDA
test_longwave_LDADD += -lcudart
endif
TESTS = test_longwave
