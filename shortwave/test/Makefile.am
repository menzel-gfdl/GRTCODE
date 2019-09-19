AM_CPPFLAGS = -I${top_builddir}/shortwave/src \
              -I${top_builddir}/rs_utils/src
if USE_SINGLE_PRECISION
AM_CPPFLAGS += -DSINGLE_PRECISION
endif

check_PROGRAMS = test_shortwave
test_shortwave_SOURCES = test_shortwave.c
test_shortwave_LDADD = ${top_builddir}/shortwave/src/libshortwave.a \
                       ${top_builddir}/rs_utils/src/librs_utils.a \
                       -lm
if ENABLE_CUDA
test_shortwave_LDADD += -lcudart
nodist_EXTRA_test_shortwave_SOURCES = dummy.cpp
endif
TESTS = test_shortwave
