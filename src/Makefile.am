## Enable single precision.
if USE_SINGLE_PRECISION
PREC_FLAG = -DSINGLE_PRECISION
PREC_TYPE = c_float
else
PREC_TYPE = c_double
endif


## Compile the c/cuda library.
if ENABLE_CUDA
CUOBJ = kernels.cu.c
endif
if CREATE_FORTRAN_BINDINGS
FLIB = libmolecular_linesf.la
endif
lib_LTLIBRARIES = libmolecular_lines.la $(FLIB)
libmolecular_lines_la_SOURCES = molecular_lines.c launch.c kernels.c kernel_utils.c spectral_bin.c ozone_continuum.c water_vapor_continuum.c RFM_voigt.c tips2017.c parse_HITRAN_file.c parse_csv.c molecules.c utils.c verbosity.c $(CUOBJ)
libmolecular_lines_ladir = $(includedir)
libmolecular_lines_la_HEADERS = molecular_lines.h floating_point_type.h molecules.h parse_HITRAN_file.h ozone_continuum.h spectral_bin.h water_vapor_continuum.h
libmolecular_lines_la_CPPFLAGS = $(PREC_FLAG)
if ENABLE_CUDA
libmolecular_lines_la_CFLAGS = -dc -x cu $(NVCCFLAGS) --compiler-options="$(CFLAGS)"
endif


## Compile the fortran bindings.
if CREATE_FORTRAN_BINDINGS
libmolecular_linesf_la_SOURCES = molecular_lines_f.F90 molecular_lines_fhl.F90
libmolecular_linesf_ladir = $(includedir)
libmolecular_linesf_la_HEADERS = molecular_lines_f.mod molecular_lines_fhl.mod
libmolecular_linesf_la_CPPFLAGS = $(PREC_FLAG)
libmolecular_linesf_la_LIBADD = libmolecular_lines.la
endif


## Build the python bindings.
if HAVE_PYTHON
pkgpython_PYTHON = molecular_lines.py __init__.py
CLEANFILES = $(FHEADERS) molecular_lines.py
EXTRA_DIST = molecular_lines.py.in
molecular_lines.py: molecular_lines.py.in
	sed -e 's/TYPE/${PREC_TYPE}/g' -e 's,LIBDIR,"${libdir}",g' $< > $@
endif
