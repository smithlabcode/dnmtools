dnl This file is part of dnmtools
dnl
dnl Copyright (C) 2022: Andrew D. Smith
dnl
dnl Authors: Andrew D. Smith
dnl
dnl This is free software: you can redistribute it and/or modify it
dnl under the terms of the GNU General Public License as published by
dnl the Free Software Foundation, either version 3 of the License, or
dnl (at your option) any later version.
dnl
dnl This software is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
dnl General Public License for more details.

AC_INIT([dnmtools], [1.2.1], [andrewds@usc.edu],
        [dnmtools], [https://github.com/smithlabcode/dnmtools])
dnl the config.h is not currently #included in the source, and only
dnl used to keep command lines short.
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([subdir-objects foreign])

AC_CONFIG_MACRO_DIR([m4])
AC_LANG(C++)
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX_11([noext], [mandatory])
AC_PROG_RANLIB

dnl OpenMP happens here
AC_OPENMP([C++])
AS_VAR_IF(OPENMP_CXXFLAGS, [], [
dnl check for the OpenMP library; can't be later
AC_SEARCH_LIBS([omp_get_num_threads], [gomp omp], [],
               [AC_MSG_FAILURE([OpenMP library not found])])
dnl now we get setup for the right OpenMP flags
ADS_OPENMP([], [AC_MSG_FAILURE([OpenMP must be installed to build dnmtools])])
])dnl end of OpenMP stuff

dnl recursively configure abismal and smithlab_cpp
AC_CONFIG_SUBDIRS([src/abismal])
AX_SUBDIRS_CONFIGURE([src/smithlab_cpp], [--enable-hts])

dnl check for HTSLib if requested
hts_fail_msg="
Failed to locate HTSLib on your system. Please use the LDFLAGS and
CPPFLAGS variables to specify the directories where the HTSLib library
and headers can be found.
"

gsl_fail_msg="
Failed to locate the GNU Scientific Library (GSL) on your
system. Please use the LDFLAGS and CPPFLAGS variables to specify the
directories where the GSL library and headers can be found.
"

zlib_fail_msg="
Failed to locate the ZLib on your system. Please use the LDFLAGS and
CPPFLAGS variables to specify the directories where the ZLib library
and headers can be found.
"

dnl check for required libraries
AC_SEARCH_LIBS([hts_version], [hts], [], [AC_MSG_FAILURE([$hts_fail_msg])])
AC_SEARCH_LIBS([zlibVersion], [z], [], [AC_MSG_FAILURE([$zlib_fail_msg])])
AC_SEARCH_LIBS([cblas_dgemm], [gslcblas], [], [AC_MSG_FAILURE([$gsl_fail_msg])])
AC_SEARCH_LIBS([gsl_blas_dgemm], [gsl], [], [AC_MSG_FAILURE([$gsl_fail_msg])])

AC_CONFIG_FILES([
  Makefile
])
AC_OUTPUT
