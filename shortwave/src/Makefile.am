AM_CPPFLAGS = -I${top_builddir}/rs_utils/src
if USE_SINGLE_PRECISION
AM_CPPFLAGS += -DSINGLE_PRECISION
endif
lib_LIBRARIES = libshortwave.a
libshortwave_a_SOURCES = rayleigh.c shortwave.c solar_flux.c
if ENABLE_CUDA
gpu_object = shortwave_gpu.$(OBJEXT)
objects = $(libshortwave_a_OBJECTS)
libshortwave_a_LIBADD = $(gpu_object)
include ${top_srcdir}/cuda.mk
endif
