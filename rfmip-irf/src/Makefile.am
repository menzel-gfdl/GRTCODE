if USE_SINGLE_PRECISION
PREC_FLAG = -DSINGLE_PRECISION
endif

bin_PROGRAMS = rfmip_irf
rfmip_irf_SOURCES = atmosphere.c rfmip-irf.c
rfmip_irf_LDADD = ${top_builddir}/shortwave/src/libshortwave.la \
                  ${top_builddir}/longwave/src/liblongwave.la \
                  ${top_builddir}/molecular_lines/src/libmolecular_lines.la \
                  ${top_builddir}/rs_utils/src/librs_utils.la \
                  -lnetcdf -lhdf5 -lhdf5_hl -lm
rfmip_irf_CPPFLAGS = $(PREC_FLAG) \
                     -I${top_builddir}/shortwave/src \
                     -I${top_builddir}/longwave/src \
                     -I${top_builddir}/molecular_lines/src \
                     -I${top_builddir}/rs_utils/src
