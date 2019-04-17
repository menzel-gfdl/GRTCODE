AC_PREREQ([2.59])
AC_INIT([gfdl-grtcode], [0.0], [gfdl.climate.model.info@noaa.gov])
AC_CANONICAL_HOST
AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE([-Wall -Werror foreign dist-zip subdir-objects])
AC_CONFIG_MACRO_DIR([m4])
AC_PROG_CC
AC_C_CONST
AC_PROG_CPP
AM_PROG_AR
LT_PREREQ([2.4])
LT_INIT()
AC_PROG_INSTALL
AC_SYS_LARGEFILE
AC_ARG_ENABLE([single_precision],
[--enable-single-precision Use single precision],
[case "${enableval}" in
  yes) single_precision=true ;;
  no)  single_precision=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-single-precision]) ;;
esac],[single_precision=false])
AM_CONDITIONAL([USE_SINGLE_PRECISION], [test x$single_precision = xtrue])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile
                 molecular_lines/Makefile molecular_lines/src/Makefile
                 rs_utils/Makefile rs_utils/src/Makefile
                 shortwave/Makefile shortwave/src/Makefile])
AC_OUTPUT
