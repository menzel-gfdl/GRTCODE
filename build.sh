#!/bin/bash -e

# Handle command line arguments.
arglist="$0 [-h|--help] env [-d]"
counter=0
while [[ $# -gt 0 ]]; do
  key=$1
  case $key in
    -h|--help)
      echo "$arglist"
      echo "\nPositional arguments:"
      echo "env:             Environment script to source."
      echo "\nOptional arguments:"
      echo "-d:              Install dependencies."
      echo "-h, --help:      Print this help message."
      exit 0
    ;;
    -d)
      install_dependencies="yes"
      shift
    ;;
    *)
      counter=$((counter+1))
      env="$key"
      shift
    ;;
  esac
done
if [ "$counter" -ne "1" ]; then
  echo "Error: incorrect number of arguments."
  echo "usage: $arglist"
  exit 1
fi

# Source environment.
source $env
here="$PWD"

if [ -n "$install_dependencies" ]; then
  # Install hdf5
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

# Build the GPU-capable version.
cd "$here"
mkdir -p bld-gpu
cd bld-gpu
CC="${CC}" CPPFLAGS="${CPPFLAGS}" CFLAGS="${CFLAGS}" CXX="${CXX}" CXXFLAGS="${CXXFLAGS}" LDFLAGS="${LDFLAGS}" NVCC="${NVCC}" NVCCFLAGS="${NVCCFLAGS}" NVCCLDFLAGS="${NVCCLDFLAGS}" MPICC="${MPICC}" MPICXX="${MPICXX}" RPATH="${RPATH}" make -f ../Makefile-gpu

# Build the host-only version.
cd "$here"
mkdir -p bld-host
cd bld-host
CC="${CC}" CPPFLAGS="${CPPFLAGS}" CFLAGS="${CFLAGS}" CXX="${CXX}" CXXFLAGS="${CXXFLAGS}" LDFLAGS="${LDFLAGS}" MPICC="${MPICC}" MPICXX="${MPICXX}" RPATH="${RPATH}" make -f ../Makefile-host
