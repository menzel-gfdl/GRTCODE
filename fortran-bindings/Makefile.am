if USE_SINGLE_PRECISION
AM_CPPFLAGS = -DSINGLE_PRECISION
endif

if FORTRAN_BINDINGS
lib_LIBRARIES = libmalloc_structs.a librs_utils_fortran.a libmolecular_lines_fortran.a
libmalloc_structs_a_SOURCES = malloc_structs.c
libmalloc_structs_a_CPPFLAGS = -I${top_builddir}/rs_utils/src \
                               -I${top_builddir}/molecular_lines/src \
                               -I${top_builddir}/shortwave/src
librs_utils_fortran_a_SOURCES = rs_utils_fortran.F90
librs_utils_fortran_a_LIBADD = ${top_builddir}/rs_utils/src/librs_utils.a \
                               ${top_builddir}/fortran-bindings/libmalloc_structs.a \
                               ${top_builddir}/shortwave/src/libshortwave.a
libmolecular_lines_fortran_a_SOURCES = molecular_lines_fortran.F90
libmolecular_lines_fortran_a_LIBADD = ${top_builddir}/molecular_lines/src/libmolecular_lines.a \
                                      ${top_builddir}/fortran-bindings/librs_utils_fortran.a \
                                      ${top_builddir}/shortwave/src/libshortwave.a
libmolecular_lines_fortran_a_CPPFLAGS = -I${top_builddir}/fortran-bindings
endif
