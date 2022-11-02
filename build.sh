#!/bin/bash -e

source gpubox.env.sh
HERE=$PWD

# Install hdf5
#INSTALL_DEPENDENCIES="yes"
if [ -n "$INSTALL_DEPENDENCIES" ]; then
  cd $PREFIX
  if [ ! -d "hdf5" ]; then
    git clone https://github.com/HDFGroup/hdf5.git
  fi
  cd hdf5
  autoreconf -i
  ./autogen.sh
  ./configure --prefix=${PREFIX} CC=${CC}
  make
  make install
  
  # Install }}netcdf
  cd $PREFIX
  if [ ! -d "netcdf-c" ]; then
    git clone https://github.com/Unidata/netcdf-c.git
  fi
  cd netcdf-c
  git checkout v4.9.0
  ./bootstrap
  ./configure --prefix=${PREFIX} CC=${CC} CPPFLAGS="${CPPFLAGS}" LDFLAGS="${LDFLAGS}"
  make
  make install
else
  echo "Dependencies are assumed to be installed in ${PREFIX}/lib"
fi

# Build the code.
cd "$HERE"
mkdir -p bld-gpu
cd bld-gpu
CC="${CC}" CPPFLAGS="${CPPFLAGS}" CFLAGS="${CFLAGS}" CXX="${CXX}" CXXFLAGS="${CXXFLAGS}" LDFLAGS="${LDFLAGS}" NVCC="${NVCC}" NVCCFLAGS="${NVCCFLAGS}" NVCCLDFLAGS="${NVCCLDFLAGS}" MPICC="${MPICC}" MPICXX="${MPICXX}" RPATH="${RPATH}" make -f ../Makefile-gpu

cd "$HERE"
mkdir -p bld-host
cd bld-host
CC="${CC}" CPPFLAGS="${CPPFLAGS}" CFLAGS="${CFLAGS}" CXX="${CXX}" CXXFLAGS="${CXXFLAGS}" LDFLAGS="${LDFLAGS}" MPICC="${MPICC}" MPICXX="${MPICXX}" RPATH="${RPATH}" make -f ../Makefile-host
