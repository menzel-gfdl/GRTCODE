if USE_SINGLE_PRECISION
PREC_FLAG = -DSINGLE_PRECISION
endif

lib_LTLIBRARIES = libmolecular_lines.la
libmolecular_lines_la_SOURCES = cfcs.c kernel_utils.c kernels.c launch.c molecular_lines.c molecules.c ozone_continuum.c parse_HITRAN_file.c RFM_voigt.c spectral_bin.c tips2017.c water_vapor_continuum.c
libmolecular_lines_la_LIBADD = ${top_builddir}/rs_utils/src/librs_utils.la
libmolecular_lines_la_CPPFLAGS = $(PREC_FLAG) -I${top_builddir}/rs_utils/src
