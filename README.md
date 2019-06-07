# GPU-able Radiative Transfer code (GRTcode)

# Overview
GRTcode consists of a column-based line-by-line radiative transfer
model capable of running on NVIDIA GPUs via Cuda-c and many-core CPUs
via OpenMP.  This package consists of several general purpose libraries
and a set of binaries that specifically target popular radiative transfer
benchmarks.

### Libraries

##### rs_utils
A utility library containing helper functions, error handlers, etc.

##### molecular_lines
A library that calculates the optical propeties (i.e. optical depths) of an
atmospheric column assumed to be made up of a series of "pressure levels".
This library uses the HITRAN database to calculate spectra for gases including
water vapor, ozone, carbon dioxide, etc. and several species of CFCs, HCFCs, and HFCs.

##### longwave
A four-stream longwave solver that uses a 2-term Pade approximation.

##### shortwave
A two-stream solver delta-Eddington and Adding methods.

### Binaries

##### rfmip-irf
Program that calculates longwave and shortwave radiative fluxes for the RFMIP-IRF
benchmark cases.

# Requirements
This library is setup to use the GNU Build System.  The libraries and binaries are
written in c (ISO/IEC 9899:1999 standard) and thus requires a c compiler,
such as the freely available gcc.  In order to run on NVIDIA GPUs,
a c++ compiler (such as g++), CUDA, and the NVCC compiler are also required.
Some of the binaries also require netCDF and HDF5.

# Building

### CPU-only
To build using default settings, run the normal GNU Build System commands:

```
$ autoreconf --install
$ ./configure [--prefix <where_to_install>] [--enable-single-precision]
$ make
$ make install
```

As usual, the default compilers and flags can be
overridden by specifying CC, CFLAGS, LDFLAGS, FC, and FFLAGS when running
configure in the usual fashion:

```
$ ./configure CC=icc CFLAGS='-O3 -qopenmp' LDFLAGS='-qopenmp'
```

This library can take advantage of thread-level parallelism through the
use of OpenMP (the default settings run single-threaded).  If your compiler
supports OpenMP, make sure you activate it by setting the corresponding
compiler option.
