AM_CPPFLAGS = -I${top_builddir}/shortwave/src \
              -I${top_builddir}/longwave/src \
              -I${top_builddir}/molecular_lines/src \
              -I${top_builddir}/rs_utils/src
if USE_SINGLE_PRECISION
AM_CPPFLAGS += -DSINGLE_PRECISION
endif
bin_PROGRAMS = rfmip_irf
rfmip_irf_SOURCES = atmosphere.c rfmip-irf.c
rfmip_irf_LDADD = ${top_builddir}/shortwave/src/libshortwave.a \
                  ${top_builddir}/longwave/src/liblongwave.a \
                  ${top_builddir}/molecular_lines/src/libmolecular_lines.a \
                  ${top_builddir}/rs_utils/src/librs_utils.a
if ENABLE_CUDA
rfmip_irf_LDADD += -lcudart
nodist_EXTRA_rfmip_irf_SOURCES = dummy.cpp
endif
