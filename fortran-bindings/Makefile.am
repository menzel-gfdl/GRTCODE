if USE_SINGLE_PRECISION
PREC_FLAG = -DSINGLE_PRECISION
endif

if FORTRAN_BINDINGS
lib_LTLIBRARIES = libmalloc_structs.la librs_utils_fortran.la libmolecular_lines_fortran.la
libmalloc_structs_la_SOURCES = malloc_structs.c
libmalloc_structs_la_CPPFLAGS = -I${top_builddir}/rs_utils/src \
                                -I${top_builddir}/molecular_lines/src \
                                -I${top_builddir}/shortwave/src
librs_utils_fortran_la_SOURCES = rs_utils_fortran.F90
librs_utils_fortran_la_LIBADD = ${top_builddir}/rs_utils/src/librs_utils.la \
                                ${top_builddir}/fortran-bindings/libmalloc_structs.la \
                                ${top_builddir}/shortwave/src/libshortwave.la
librs_utils_fortran_la_CPPFLAGS = $(PREC_FLAG)
libmolecular_lines_fortran_la_SOURCES = molecular_lines_fortran.F90
libmolecular_lines_fortran_la_LIBADD = ${top_builddir}/molecular_lines/src/libmolecular_lines.la \
                                       ${top_builddir}/fortran-bindings/librs_utils_fortran.la \
                                       ${top_builddir}/shortwave/src/libshortwave.la
libmolecular_lines_fortran_la_CPPFLAGS = $(PREC_FLAG) -I${top_builddir}/fortran-bindings
endif
