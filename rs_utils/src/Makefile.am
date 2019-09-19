if USE_SINGLE_PRECISION
AM_CPPFLAGS = -DSINGLE_PRECISION
endif
lib_LIBRARIES = librs_utils.a
librs_utils_a_SOURCES = argparse.c curtis_godson.c device.c optics.c parse_csv.c \
                        spectral_grid.c utils.c verbosity.c
if ENABLE_CUDA
gpu_object = rs_utils_gpu.$(OBJEXT)
objects = $(librs_utils_a_OBJECTS)
librs_utils_a_LIBADD = $(gpu_object)
include ${top_srcdir}/cuda.mk
endif
