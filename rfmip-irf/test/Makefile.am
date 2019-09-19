AM_CPPFLAGS = -I${top_builddir}/rs_utils/src
if USE_SINGLE_PRECISION
AM_CPPFLAGS += -DSINGLE_PRECISION
endif

check_PROGRAMS = check_results
check_results_SOURCES = check_results.c
check_results_LDADD = ${top_builddir}/rs_utils/src/librs_utils.a \
                      -lnetcdf -lhdf5 -lhdf5_hl -lm
TESTS = test_rfmip_irf
