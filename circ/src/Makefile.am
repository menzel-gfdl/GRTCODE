AM_CPPFLAGS = -I${top_builddir}/shortwave/src \
              -I${top_builddir}/longwave/src \
              -I${top_builddir}/molecular_lines/src \
              -I${top_builddir}/rs_utils/src
if USE_SINGLE_PRECISION
AM_CPPFLAGS += -DSINGLE_PRECISION
endif
if USE_DISORT
AM_CPPFLAGS += -DDISORT
endif
bin_PROGRAMS = circ
circ_SOURCES = utilities.c cloud_optics.c disort_shortwave.c atmosphere.c circ.c
circ_LDADD = ${top_builddir}/shortwave/src/libshortwave.a \
             ${top_builddir}/longwave/src/liblongwave.a \
             ${top_builddir}/molecular_lines/src/libmolecular_lines.a \
             ${top_builddir}/rs_utils/src/librs_utils.a \
             -lnetcdf -lm
if ENABLE_CUDA
circ_LDADD += -lcudart
nodist_EXTRA_circ_SOURCES = dummy.cpp
endif
