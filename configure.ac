AC_PREREQ([2.59])
AC_INIT([gfdl-grtcode], [0.0], [gfdl.climate.model.info@noaa.gov])
AC_CANONICAL_HOST
AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE([-Wall -Werror foreign dist-zip subdir-objects])
AC_CONFIG_MACRO_DIR([m4])
AC_PROG_CC
AC_C_CONST
AC_PROG_CPP
AC_PROG_CXX
AM_PROG_AR
AC_PROG_RANLIB
AC_PROG_INSTALL
AC_SYS_LARGEFILE

#Provide an option to build in single precision.
AC_ARG_ENABLE([single_precision],
[--enable-single-precision Use single precision],
[case "${enableval}" in
  yes) single_precision=true ;;
  no)  single_precision=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-single-precision]) ;;
esac],[single_precision=false])
AM_CONDITIONAL([USE_SINGLE_PRECISION], [test x$single_precision = xtrue])

#Provide an option to build with CUDA.
AC_ARG_VAR([NVCC], [nvcc compiler command])
AC_ARG_VAR([NVCCFLAGS], [nvcc compiler flags])
AC_ARG_VAR([NVCCLDFLAGS], [nvcc linker flags])
AC_ARG_ENABLE([cuda],
[--enable-cuda Enable ability to run on NVIDIA GPUs.],
[case "${enableval}" in
  yes) cuda=true ;;
  no)  cuda=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-cuda]) ;;
esac],[cuda=false])
AM_CONDITIONAL([ENABLE_CUDA], [test x$cuda = xtrue])
if test x$cuda = xtrue
then
  AC_PATH_PROG([NVCC], [nvcc], [no])
  if test x$NVCC != xno
  then
    AC_MSG_CHECKING([whether nvcc works])
    cat > conftest.cu <<EOF
__global__ void test(void)
{
    int tid = blockIdx.x*blockDim.x + threadIdx.x;
    return;
}
int main(void)
{
    test<<<1, 32, 0, 0>>>();
    return 0;
}
EOF
    if $NVCC -x cu -w -o conftest.x conftest.cu && ./conftest.x
    then
      AC_MSG_RESULT([yes])
      rm -f conftest.cu conftest.x
    else
      AC_MSG_RESULT([no])
      rm -f conftest.cu conftest.x
      AC_MSG_ERROR([nvcc could not compile/run executable])
    fi
  else
    AC_MSG_ERROR([could not find nvcc])
  fi
fi

#Provide an option to build fortran bindings.
AC_ARG_ENABLE([fortran],
[--enable-fortran Build fortran bindings],
[case "${enableval}" in
  yes) fortran=true ;;
  no)  fortran=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-fortran]) ;;
esac],[fortran=false])
AM_CONDITIONAL([FORTRAN_BINDINGS], [test x$fortran = xtrue])
AC_PROG_FC

#Provide an option to use DISORT.
AC_ARG_ENABLE([use_disort],
[--enable-disort Use disort solver],
[case "${enableval}" in
  yes) use_disort=true ;;
  no)  use_disort=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-disort]) ;;
esac],[use_disort=false])
AM_CONDITIONAL([USE_DISORT], [test x$use_disort = xtrue])

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile
                 rs_utils/Makefile rs_utils/src/Makefile
                 molecular_lines/Makefile molecular_lines/src/Makefile
                 shortwave/Makefile shortwave/src/Makefile shortwave/test/Makefile
                 longwave/Makefile longwave/src/Makefile longwave/test/Makefile
                 rfmip-irf/Makefile rfmip-irf/src/Makefile rfmip-irf/test/Makefile
                 circ/Makefile circ/src/Makefile circ/test/Makefile
                 fortran-bindings/Makefile])
AC_OUTPUT
