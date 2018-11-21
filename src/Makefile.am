## Build the c library (and optional fortran bindings).
if USE_SINGLE_PRECISION
  PREC_FLAG = -DSINGLE_PRECISION
  PREC_TYPE = c_float
else
  PREC_TYPE = c_double
endif
if CREATE_FORTRAN_BINDINGS
  FOBJ = molecular_lines_f.F90 molecular_lines_fhl.F90
  FHEADERS = molecular_lines_f.mod molecular_lines_fhl.mod
endif
lib_LTLIBRARIES = libmolecular_lines.la
libmolecular_lines_la_SOURCES = molecular_lines.c host_launch.c calc_optical_depth.c spectral_bin.c ozone_continuum.c water_vapor_continuum.c RFM_voigt.c tips2017.c parse_HITRAN_file.c parse_csv.c molecules.c utils.c verbosity.c ${FOBJ}
libmolecular_lines_ladir = ${includedir}
libmolecular_lines_la_HEADERS = molecular_lines.h floating_point_type.h molecules.h parse_HITRAN_file.h ozone_continuum.h spectral_bin.h water_vapor_continuum.h $(FHEADERS)
libmolecular_lines_la_CPPFLAGS = ${PREC_FLAG}

## Build the python bindings.
if HAVE_PYTHON
pkgpython_PYTHON = molecular_lines.py __init__.py
endif
CLEANFILES = $(FHEADERS) molecular_lines.py
EXTRA_DIST = molecular_lines.py.in
molecular_lines.py: molecular_lines.py.in
	sed -e 's/TYPE/${PREC_TYPE}/g' -e 's,LIBDIR,"${libdir}",g' $< > $@
