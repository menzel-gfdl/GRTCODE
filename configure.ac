AC_INIT([molecular_lines], [0.0], [gfdl.climate.model.info@noaa.gov])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_PROG_CC
AC_PROG_FC
AM_PROG_AR
AC_PROG_LIBTOOL


AC_ARG_ENABLE([fortran_bindings],
[--enable-fortran-bindings Create Fortran bindings],
[case "${enableval}" in
  yes) fortran_bindings=true ;;
  no)  fortran_bindings=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-fortran-bindings]) ;;
esac],[fortran_bindings=false])
AM_CONDITIONAL([CREATE_FORTRAN_BINDINGS], [test x$fortran_bindings = xtrue])


AC_ARG_ENABLE([single_precision],
[--enable-single-precision Use single precision],
[case "${enableval}" in
  yes) single_precision=true ;;
  no)  single_precision=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-single-precision]) ;;
esac],[single_precision=false])
AM_CONDITIONAL([USE_SINGLE_PRECISION], [test x$single_precision = xtrue])


AM_PATH_PYTHON([3.5],,[:])
AM_CONDITIONAL([HAVE_PYTHON], [test "$PYTHON" != :])


AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
 Makefile
 src/Makefile
])
AC_OUTPUT
