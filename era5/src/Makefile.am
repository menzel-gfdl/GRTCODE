AM_CPPFLAGS = -I${top_builddir}/shortwave/src \
              -I${top_builddir}/longwave/src \
              -I${top_builddir}/gas-optics/src \
              -I${top_builddir}/utilities/src
if USE_SINGLE_PRECISION
AM_CPPFLAGS += -DSINGLE_PRECISION
endif
bin_PROGRAMS = era5
era5_SOURCES = atmosphere.c era5.c
era5_LDADD = ${top_builddir}/shortwave/src/libshortwave.a \
             ${top_builddir}/longwave/src/liblongwave.a \
             ${top_builddir}/gas-optics/src/libgas_optics.a \
             ${top_builddir}/utilities/src/libgrtcode_utilities.a \
             ${top_builddir}/utilities/src/libargparse.a
if ENABLE_CUDA
era5_LDADD += -lcudart
nodist_EXTRA_era5_SOURCES = dummy.cpp
endif
