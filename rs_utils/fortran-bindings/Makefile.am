if USE_SINGLE_PRECISION
PREC_FLAG = -DSINGLE_PRECISION
endif

if FORTRAN_BINDINGS
lib_LTLIBRARIES = libmalloc_structs.la librs_utils_fortran.la
libmalloc_structs_la_SOURCES = malloc_structs.c
libmalloc_structs_la_CPPFLAGS = -I${top_builddir}/rs_utils/src \
                                -I${top_builddir}/molecular_lines/src
librs_utils_fortran_la_SOURCES = rs_utils_fortran.F90
librs_utils_fortran_la_LIBADD = ${top_builddir}/rs_utils/src/librs_utils.la \
                                ${top_builddir}/rs_utils/fortran-bindings/libmalloc_structs.la
librs_utils_fortran_la_CPPFLAGS = $(PREC_FLAG)
endif
