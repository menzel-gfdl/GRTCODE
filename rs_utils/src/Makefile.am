if USE_SINGLE_PRECISION
PREC_FLAG = -DSINGLE_PRECISION
endif

lib_LTLIBRARIES = librs_utils.la
librs_utils_la_SOURCES = argparse.c curtis_godson.c device.c optics.c parse_csv.c spectral_grid.c utils.c verbosity.c
librs_utils_la_CPPFLAGS = $(PREC_FLAG)
