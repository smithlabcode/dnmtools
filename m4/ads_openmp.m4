#
# SYNOPSIS
#
#   ADS_OPENMP([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   This macro is copied from the one in autoconf-archive with
#   modifications to work with the clang of some macOS versions as of
#   2022. The documentation is mostly copied from the original, with
#   the initials ADS where the new stuff appears.
#
#   This macro tries to find out how to compile programs that use OpenMP a
#   standard API and set of compiler directives for parallel programming
#   (see http://www-unix.mcs/)
#
#   On success, it sets the OPENMP_CFLAGS/OPENMP_CXXFLAGS/OPENMP_F77FLAGS
#   output variable to the flag (e.g. -omp) used both to compile *and* link
#   OpenMP programs in the current language.
#
#   NOTE: You are assumed to not only compile your program with these flags,
#   but also link it with them as well.
#
#   If you want to compile everything with OpenMP, you should set:
#
#     CFLAGS="$CFLAGS $OPENMP_CFLAGS"
#     #OR#  CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"
#     #OR#  FFLAGS="$FFLAGS $OPENMP_FFLAGS"
#
#   (depending on the selected language).
#
#   The user can override the default choice by setting the corresponding
#   environment variable (e.g. OPENMP_CFLAGS).
#
#   ACTION-IF-FOUND is a list of shell commands to run if an OpenMP flag is
#   found, and ACTION-IF-NOT-FOUND is a list of commands to run it if it is
#   not found. If ACTION-IF-FOUND is not specified, the default action will
#   define HAVE_OPENMP.
#
#   (ADS) This macro attempts to take care of both compiling and
#   linking, which means you will need an appropriate library in your
#   LIBS already. This can be accomplished by using AC_SEARCH_LIBS
#   with the function `omp_get_num_threads` and libraries like gomp
#   and omp.
#
#   (ADS) OpenMP is not enabled unless the OPENMP_FLAGS are set as
#   outputs, which can be done with AC_SUBST or by calling AC_OPENMP
#   first. The latter might allow for more elegant handling of code
#   that should run either as single or multithreaded.
#
# LICENSE
#
#   Copyright (c) 2008 Steven G. Johnson <stevenj@alum.mit.edu>
#   Copyright (c) 2015 John W. Peterson <jwpeterson@gmail.com>
#   Copyright (c) 2016 Nick R. Papior <nickpapior@gmail.com>
#   Copyright (c) 2022 Andrew D. Smith <andrewds@usc.edu>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <https://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

#serial 14 garbage

AC_DEFUN([ADS_OPENMP], [
AC_PREREQ([2.69]) dnl for _AC_LANG_PREFIX

AC_CACHE_CHECK([for OpenMP flag of _AC_LANG compiler], ax_cv_[]_AC_LANG_ABBREV[]_openmp,
[
AC_REQUIRE([_ADS_OPENMP_SAFE_WD])
save[]_AC_LANG_PREFIX[]FLAGS=$[]_AC_LANG_PREFIX[]FLAGS
ax_cv_[]_AC_LANG_ABBREV[]_openmp=unknown
# Flags to try:  -fopenmp (gcc), -mp (SGI & PGI),
#                -qopenmp (icc>=15), -openmp (icc),
#                -xopenmp (Sun), -omp (Tru64),
#                -qsmp=omp (AIX),
#                -Wp,-fopenmp (macOS)
#                none
ax_openmp_flags="-fopenmp -openmp -qopenmp -mp -xopenmp -omp -qsmp=omp -Wp,-fopenmp none"
if test "x$OPENMP_[]_AC_LANG_PREFIX[]FLAGS" != x; then
  ax_openmp_flags="$OPENMP_[]_AC_LANG_PREFIX[]FLAGS $ax_openmp_flags"
fi
for ax_openmp_flag in $ax_openmp_flags; do
  case $ax_openmp_flag in
    none) []_AC_LANG_PREFIX[]FLAGS=$save[]_AC_LANG_PREFIX[] ;;
    *) []_AC_LANG_PREFIX[]FLAGS="$save[]_AC_LANG_PREFIX[]FLAGS $ax_openmp_flag" ;;
  esac
  AC_LINK_IFELSE([AC_LANG_SOURCE([[
@%:@include <omp.h>

static void
parallel_fill(int * data, int n)
{
  int i;
@%:@pragma omp parallel for
  for (i = 0; i < n; ++i)
    data[i] = i;
}

int
main()
{
  int arr[100000];
  omp_set_num_threads(2);
  parallel_fill(arr, 100000);
  return 0;
}
]])],[ax_cv_[]_AC_LANG_ABBREV[]_openmp=$ax_openmp_flag; break],[])
done
[]_AC_LANG_PREFIX[]FLAGS=$save[]_AC_LANG_PREFIX[]FLAGS
rm -rf mp penmp mp.dSYM penmp.dSYM
])
if test "x$ax_cv_[]_AC_LANG_ABBREV[]_openmp" = "xunknown"; then
  m4_default([$2],:)
else
  if test "x$ax_cv_[]_AC_LANG_ABBREV[]_openmp" != "xnone"; then
    OPENMP_[]_AC_LANG_PREFIX[]FLAGS=$ax_cv_[]_AC_LANG_ABBREV[]_openmp
  fi
  m4_default([$1], [AC_DEFINE(HAVE_OPENMP,1,[Define if OpenMP is enabled])])
fi
])dnl ADS_OPENMP

# _ADS_OPENMP_SAFE_WD
# ------------------
# AC_REQUIREd by ADS_OPENMP.  Checks both at autoconf time and at
# configure time for files that ADS_OPENMP clobbers.
AC_DEFUN([_ADS_OPENMP_SAFE_WD],
[m4_syscmd([test ! -e penmp && test ! -e mp && test ! -e penmp.dSYM && test ! -e mp.dSYM])]dnl
[m4_if(sysval, [0], [], [m4_fatal(m4_normalize(
  [ADS_OPENMP clobbers files named 'mp' and 'penmp'.
   To use ADS_OPENMP you must not have either of these files
   at the top level of your *build* tree.]))])]dnl
[if test -e penmp || test -e mp || test -e penmp.dSYM || test -e mp.dSYM ; then
  AC_MSG_ERROR(m4_normalize(
    [ADS@&t@_OPENMP clobbers files named 'mp' and 'penmp', and
     directoreis named 'mp.dSYM' and 'penmp.dSYM'.
     Aborting configure because one of these files already exists.]))
fi])
