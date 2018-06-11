#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([skgen], [0.0.1], [BUG-REPORT-ADDRESS])
AC_CONFIG_SRCDIR([scr/main2.cpp])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_MACRO_DIR([m4])

# Checks for programs.
AC_PROG_CXX
AC_PROG_RANLIB
# Checks for libraries.

# Checks for header files.
#AC_CHECK_LIBS([xc openblas gsl m])

AC_CHECK_HEADERS([stdlib.h string.h])
AX_BOOST_BASE([1.53],,[AC_MSG_ERROR([skgen needs Boost, but it was not found in your system])])

export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:"/opt/etsf/lib/pkgconfig/"
PKG_CHECK_MODULES([XC], [libxc <= 3.0.0])
PKG_CHECK_MODULES([GSL], [gsl])
# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL

# Checks for library functions.
AC_CHECK_FUNCS([pow sqrt])

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([scr/Makefile])
AC_OUTPUT
