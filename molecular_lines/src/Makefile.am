AM_CPPFLAGS = -I${top_builddir}/rs_utils/src
if USE_SINGLE_PRECISION
AM_CPPFLAGS += -DSINGLE_PRECISION
endif
if ENABLE_CUDA
cuda_objects = cuda_kernels.c
endif
lib_LIBRARIES = libmolecular_lines.a
libmolecular_lines_a_SOURCES = cfcs.c collision_induced_absorption.c kernel_utils.c kernels.c \
                               launch.c molecular_lines.c molecules.c ozone_continuum.c \
                               parse_HITRAN_file.c RFM_voigt.c spectral_bin.c tips2017.c \
                               water_vapor_continuum.c $(cuda_objects)
if ENABLE_CUDA
gpu_object = molecular_lines_gpu.$(OBJEXT)
objects = $(libmolecular_lines_a_OBJECTS)
libmolecular_lines_a_LIBADD = $(gpu_object)
include ${top_srcdir}/cuda.mk
endif
