if USE_SINGLE_PRECISION
PREC_FLAG = -DSINGLE_PRECISION
endif

check_PROGRAMS = check_results
check_results_SOURCES = check_results.c
check_results_LDADD = ${top_builddir}/rs_utils/src/librs_utils.la \
                      -lnetcdf -lhdf5 -lhdf5_hl -lm
check_results_CPPFLAGS = $(PREC_FLAG) \
                         -I${top_builddir}/rs_utils/src
TESTS = test_rfmip_irf
