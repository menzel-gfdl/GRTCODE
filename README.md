# GRTcode Molecular Lines


# Overview
This c library calculates optical depth values for atmospheric
columns on GPU and CPU architectures by explicitly treating the
rotational-vibrational lines of atmospheric molecules.  Molecular
line parameters are passed into the model as ascii HITRAN database
files, and optical depth values in each layer of the atmospheric
column are calculated on each point of a user-defined input spectral
grid.  Both FORTRAN-90 and Python (>= 3.5) bindings are optionally included.


# Requirements
This library is setup to use the GNU Build System.  The core library is
written in c (ISO/IEC 9899:1999 standard) and thus requires a c compiler,
such as the freely available gcc.  In addition, a set of FORTRAN
(ISO/IEC TS 29113:2012) and Python (3.5+) bindings are included, which
require a fortran compiler (such as gfortran) and the python numpy module
respectively.  In order to run on a NVIDIA GPU, the library also requires a
c++ compiler (such as g++), CUDA, and the NVCC compiler.


# Source Code
The source code currently resides in
[this Gitlab repository](https://gitlab.gfdl.noaa.gov/Raymond.Menzel/grtcodev2),
on branch use_autotools.  To obtain the code, run

```
$ git clone https://gitlab.gfdl.noaa.gov/Raymond.Menzel/grtcodev2.git molecular_lines
$ cd molecular_lines
$ git checkout use_autotools
```


# Building

### CPU-only
To build using default settings, run the normal GNU Build System commands:

```
$ autoreconf --install
$ ./configure [--prefix <where_to_install>] [--enable-fortran-bindings] \
              [--enable-single-precision]
$ make
$ make install
```

If you have Python (3.5+) installed on your system, python bindings will
also be installed.  These bindings require numpy, as well as setting the
PYTHONPATH environment variable.  For example, in bash:

```
export PYTHONPATH="<prefix>/lib/python3.5/site-packages/molecular_lines"
```

As usual, the default compilers and flags can be
overridden by specifying CC, CFLAGS, LDFLAGS, FC, and FFLAGS when running
configure in the usual fashion:

```
$ ./configure CC=icc CFLAGS='-O3 -qopenmp' LDFLAGS='-qopenmp' \
              --enable-fortran-bindings FC=ifort FFLAGS='-O3 -qopenmp'
```

This library can take advantage of thread-level parallelism through the
use of OpenMP (the default settings run single-threaded).  If your compiler
supports OpenMP, make sure you activate it by setting the corresponding
compiler option.

### With GPUs
If you have a c++ compiler, CUDA, and the NVCC compiler installed, the library
can built to run on NVIDIA GPUs by running

```
$ make -f Makefile.nvcc
$ make -f Makefile.nvcc install
```

Once again it is recommended that you also run

```
$ make -f Makefile.nvcc test
```

The make variable PREFIX can be used to specify where the
library will be installed (default is the current directory).

### Tip
One way to see how many CUDA-enabled GPUs your system has is to try running:

```
$ nvidia-smi --list-gpus
```
Each CUDA-enabled GPU that is found will be listed with with its device id,
followed by its model name.  Example output on a system with two Tesla K40c
GPUs will look like:

```
GPU 0: Tesla K40c (UUID: GPU-xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx)
GPU 1: Tesla K40c (UUID: GPU-xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx)
```


# APIs
A C API, Low-level Fortran API, and High-level Fortran API
are provided.  The Low-level Fortran Api simply provides direct
bindings to the C API, while the High-level Fortran API
provides a more simplified interface (at the expense of less control).
All APIs are documented with Doxygen.  To view Doxygen-generated HTML
describing each API, please open the file APIs.html in a browser and click
on the Modules tab.


# Extras
In addition to the source code and Makefiles, a few other directories
containing example programs and input data files are included in the
base of this repository.  As a warning, the provided example programs
make use of these provided input data directories, so any change to
their paths or file names will prevent the example programs from
running properly.


### Example Codes
A short, complete example program for each of the three APIs described
above is included in the examples directory.  These examples
demonstrate how to properly use each library function, and can be
built by simply running

```
$ make test
```

or

```
$ make -f Makefile.nvcc test
```


### Example HITRAN Database Files
Example ASCII HITRAN database file for water vapor, carbon dixoide,
ozone, nitrous oxide, methane, carbon monoxide, and oxygen are
included in the HITRAN_files directory.  As described on the API
pages, input ASCII HITRAN database files for each molecule that will be
included in the optical depth calculation are required.  These files
contain all molecular lines contained in the HITRAN 2012 database
for each species, and are recommended for use by new users or users
who are not interested in customizing the molecular spectra (i.e., by
removing certain lines or adding new lines).


### Example Continua Input Files
Example data files required when running with the water vapor and
ozone continua are provided in the water_vapor_continnum and
ozone_continumm directories respectively.  As described on the API
pages, the path to each of these directories is required in order
for the continua effects to be included in the optical depth
calculations.  These directories are recommended for use by users
who want to run with standard continua values.  Please note that if
the user would like to provide their own continuum files, the names
and format (csv with the same number of columns) of the new files must
match the names of the files in these provided directories (this
restriction will probably be removed in a future release).


### Example Namelist
An example namelist (for use with the High-level Fortran API) is
include in the namelist directory.  As described on the
High-level Fortran API page, when using this API a path to a
namelist file is required to override default values.  Users of the
High-level Fortran API are encouraged to use this file as template
when creating their own input namelist files.
