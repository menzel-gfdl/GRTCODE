if USE_SINGLE_PRECISION
PREC_FLAG = -DSINGLE_PRECISION
endif

if FORTRAN_BINDINGS
lib_LTLIBRARIES = libmolecular_lines_fortran.la
libmolecular_lines_fortran_la_SOURCES = molecular_lines_fortran.F90
libmolecular_lines_fortran_la_LIBADD = ${top_builddir}/molecular_lines/src/libmolecular_lines.la \
                                       ${top_builddir}/rs_utils/fortran-bindings/librs_utils_fortran.la
libmolecular_lines_fortran_la_CPPFLAGS = $(PREC_FLAG) -I${top_builddir}/rs_utils/fortran-bindings
endif
