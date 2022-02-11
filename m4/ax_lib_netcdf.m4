_ldflags_# ===========================================================================
#      https://www.gnu.org/software/autoconf-archive/ax_lib_netcdf4.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_LIB_NETCDF()
#   AX_LIB_NETCDFF()
#
# DESCRIPTION
#
#   This macro provides tests of the availability of the NetCDF v4 c and Fortran
#   libraries.
#
#   If NetCDF c is successfully found, this macro calls
#
#     AC_SUBST(NETCDF_VERSION)
#     AC_SUBST(NETCDF_CC)
#     AC_SUBST(NETCDF_CFLAGS)
#     AC_SUBST(NETCDF_CPPFLAGS)
#     AC_SUBST(NETCDF_LDFLAGS)
#     AC_SUBST(NETCDF_LIBS)
#
#   If NetCDF Fortran library is successfully found, this macro calls
#
#     AC_SUBST(NETCDFF_VERSION)
#     AC_SUBST(NETCDF_FC)
#     AC_SUBST(NETCDF_FCFLAGS)
#     AC_SUBST(NETCDF_FCCPPFLAGS)
#     AC_SUBST(NETCDF_FCLIBS)
#     AC_DEFINE(HAVE_NETCDF)
#     AC_DEFINE(HAVE_NETCDFF)
#
#   It also sets
#
#     with_netcdf="yes"
#     with_netcdf_fortran="yes"    (if NetCDF has Fortran support)
#
#   If NetCDF is disabled or not found, this macros sets
#
#     with_netcdf="no"
#     with_netcdf_fortran="no"
#
#   Your configuration script can test $with_netcdf to take any further
#   actions. NETCDF_{C,CPP,LD}FLAGS may be used when building with C or
#   C++. NETCDF_FC{FLAGS,CPPFLAGS,LIBS} and NETCDF_LDFLAGS should be used when
#   building Fortran applications.
#
#   To use the macro, one would code one of the following in "configure.ac"
#   before AC_OUTPUT:
#
#     1) # Check for NetCDF support
#        AX_LIB_NETCDF()
#
#     2) # Check for NetCDF Fortran support
#        AX_LIB_NETCDFF()
#
#   One could test $with_netcdf for the outcome or display it as follows
#
#     echo "NetCDF v4 support:  $with_netcdf"
#
#   One could also for example, override the default CC in "configure.ac" to
#   enforce compilation with the compiler that NetCDF v4 was built with:
#
#     AX_LIB_NETCDF()
#     if test "$with_netcdf" = "yes"; then
#             CC="$NETCDF_CC"
#     else
#             AC_MSG_ERROR([Unable to find NetCDF.])
#     fi
#
# LICENSE
#
#   Copyright (c) 2016 Timothy Brown <tbrown@freeshell.org>
#                 2018 Seth Underwood <underwoo@underwoo.io>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

AC_DEFUN([AX_LIB_NETCDF], [
AC_LANG_ASSERT(C)dnl
AC_REQUIRE([AC_PROG_SED])dnl
AC_REQUIRE([AC_PROG_AWK])dnl
AC_REQUIRE([AC_PROG_GREP])dnl

AC_ARG_WITH([netcdf],
[AS_HELP_STRING(
  [--with-netcdf=[yes/no/PATH]],
  [base directory of NetCDF installation])],
[],
[with_netcdf="yes"])
AS_IF([test $with_netcdf != no],
  [_AX_C_LIB_NETCDF()])dnl
])

AC_DEFUN([_AX_C_LIB_NETCDF], [
# Was a path to NetCDF given?
AS_IF([test $with_netcdf != yes], [_ax_c_lib_netcdf_prefix=$with_netcdf])
test -z ${_ax_c_lib_netcdf_prefix+x} && _ax_c_lib_netcdf_bin=$_ax_c_lib_netcdf_prefix/bin
# Check to see if NC_CONFIG is in the path, or NETCDF_PREFIX/bin.
AC_PATH_PROGS([NC_CONFIG], [nc-config], [], [$_ax_c_lib_netcdf_bin$PATH_SEPARATOR$PATH])
# If nc-config found, use it to get information on the NetCDF library
AS_IF([test ! -z "${NC_CONFIG+x}"],[
  _ax_c_lib_netcdf_cc=$(eval $NC_CONFIG --cc 2> /dev/null | $AWK '{print $[]1}')
  AS_IF([test "x$_ax_c_lib_netcdf_cc" = "xccache"],
        [_ax_c_lib_netcdf_cc=$(eval $NC_CONFIG --cc 2> /dev/null | $AWK '{print $[]2}')])
  _ax_c_lib_netcdf_version=$(eval $NC_CONFIG --version 2> /dev/null | $AWK '{print $[]2}')
  test -z $_ax_c_lib_netcdf_prefix && _ax_c_lib_netcdf_prefix=$(eval $NC_CONFIG --prefix 2> /dev/null)
  _ax_c_lib_netcdf_cflags=$(eval $NC_CONFIG --cflags 2> /dev/null)
  # Some systems (e.g. Cray) use a compiler wrapper that automatically add the required flag, which
  # has nc-config --cflags not return the -I flag.  This causes issues if trying to not use the
  # wrapper.  Thus, we add the -I flag if the cflags variable doesn't have it.  This should be
  # harmless.
  _ax_c_lib_netcdf_includedir=$(eval $NC_CONFIG --includedir 2> /dev/null)
  echo "$_ax_c_lib_netcdf_cflags" | $GREP -e "-I$_ac_c_lib_netcdf_includedir" 2>&1 > /dev/null || _ax_c_lib_netcdf_cflags="$_ax_c_lib_netcdf_cflags -I$_ax_c_lib_netcdf_includedir"
  # LIBS and LDFLAGS sorted based on prefix (e.g. -L and -l)
  for arg in $(eval $NC_CONFIG --libs 2> /dev/null) ; do :
    AS_CASE([$arg],
            [-L*], [echo $_ax_c_lib_netcdf_ldflags | $GREP -e "$arg" 2>&1 >/dev/null || _ax_c_lib_netcdf_ldflags="$arg $_ax_c_lib_netcdf_ldflags"],
            [-l*], [echo $_ax_c_lib_netcdf_libs | $GREP -e "$arg" 2>&1 >/dev/null || _ax_c_lib_netcdf_libs="$arg $_ax_c_lib_netcdf_libs"])
  done])
AS_IF([test ! -z "${_ax_c_lib_netcdf_prefix}"], [
  # Check if nc-config was able to give us specific information.  If not,
  # Guess the locations for libs and headers, and libs flags
  test -z "$_ax_c_lib_netcdf_cflags" && _ax_c_lib_netcdf_cflags=-I$_ax_c_lib_netcdf_prefix/include
  test -z "$_ax_c_lib_netcdf_ldflags" && _ax_c_lib_netcdf_ldflags=-L$_ax_c_lib_netcdf_prefix/lib
  test -z "$_ax_c_lib_netcdf_libs" && _ax_c_lib_netcdf_libs=-lnetcdf
])
AS_VAR_SET_IF([_ax_c_lib_netcdf_cc],
              [AS_VAR_SET([NETCDF_CC],[$_ax_c_lib_netcdf_cc])],
              [AS_VAR_SET([NETCDF_CC],[unknown])])
# Check for the netCDF library
# This is essentially the AC_SEARCH_LIBS, but as we need to loop through
# different libs and ldflag combinations to verify if -lnetcdf works.
AS_VAR_PUSHDEF([ax_netcdf_lib_search], [ax_cv_lib_netcdf])dnl
AS_VAR_PUSHDEF([ax_netcdf_ldflags_search], [ax_cv_ldflags_netcdf])
AC_CACHE_CHECK([for C netCDF library], [ax_netcdf_lib_search], [
_ax_c_lib_netcdf_save_LIBS=$LIBS
_ax_c_lib_netcdf_save_LDFLAGS=$LDFLAGS
# Did the user set the NETCDF_LDFLAGS and/or NETCDF_LIBS?
test -z "$NETCDF_LDFLAGS" && NETCDF_LDFLAGS=$_ax_c_lib_netcdf_ldflags
test -z "$NETCDF_LIBS" && NETCDF_LIBS=$_ax_c_lib_netcdf_libs
AC_LANG_CONFTEST([AC_LANG_CALL([], [nc_inq_libvers])])
for ax_netcdf_ldflag in '' "$NETCDF_LDFLAGS"; do
  AS_IF([test -z "$ax_netcdf_ldflag"],
        [ax_netcdf_ldflags_res="none required" && LDFLAGS="$_ax_c_lib_netcdf_save_LDFLAGS"],
        [ax_netcdf_ldflags_res="$ax_netcdf_ldflag" && LDFLAGS="$ax_netcdf_ldflag $_ax_c_lib_netcdf_save_LDFLAGS"])
  for ax_netcdf_lib in '' "$NETCDF_LIBS" '-lnetcdf' "-lnetcdf $NETCDF_LIBS"; do
    AS_IF([test -z "$ax_netcdf_lib"],
          [ax_netcdf_lib_res="none required" && LIBS="$_ax_c_lib_netcdf_save_LIBS"],
          [ax_netcdf_lib_res=$ax_netcdf_lib && LIBS="$ax_netcdf_lib $_ax_c_lib_netcdf_save_LIBS"])
    AC_LINK_IFELSE([], [
      AS_VAR_SET([ax_netcdf_lib_search], [$ax_netcdf_lib_res])
      AS_VAR_SET([ax_netcdf_ldflags_search], [$ax_netcdf_ldflags_res])])
    AS_VAR_SET_IF([ax_netcdf_lib_search], [break 2])
  done
done
AS_VAR_SET_IF([ax_netcdf_lib_search], [], [AS_VAR_SET([ax_netcdf_lib_search], [no])])
AS_VAR_SET_IF([ax_netcdf_ldflags_search], [], [AS_VAR_SET([ax_netcdf_ldflags_search], [no])])
rm conftest.$ac_ext
LIBS=$_ax_c_lib_netcdf_save_LIBS
LDFLAGS=$_ax_c_lib_netcdf_save_LDFLAGS
AS_VAR_COPY([ax_netcdf_lib_res], [ax_netcdf_lib_search])
AS_VAR_COPY([ax_netcdf_ldflags_res], [ax_netcdf_ldflags_search])
AS_IF([test "$ax_netcdf_lib_res" != no && test "$ax_netcdf_ldflags_res" != no],
      [test "$ax_netcdf_lib_res" = "none required" || NETCDF_LIBS=$ax_netcdf_lib_res
       test "$ax_netcdf_ldflags_res" = "none required" || NETCDF_LDFLAGS=$ax_netcdf_ldflags_res
       ax_netcdf_lib_search="yes"])
])

# Check for the netcdf.h file
# This test is essentially a copy of AC_CHECK_HEADER, but looping over different
# CFLAGS combinations.  We are not concerning ourselves with CPPFLAGS, as
# NETCDF_CPPFLAGS will be the same as NETCDF_CFLAGS.
AS_VAR_PUSHDEF([ax_netcdf_cflags_search], [ax_cv_cflags_netcdf])
AC_CACHE_CHECK([netcdf.h usability], [ax_netcdf_cflags_search], [
_ax_c_lib_netcdf_save_CFLAGS=$CFLAGS
# Check if the user already set NETCDF_CFLAGS
test -z "$NETCDF_CFLAGS" && NETCDF_CFLAGS=$_ax_c_lib_netcdf_cflags
AC_LANG_CONFTEST([AC_LANG_SOURCE([@%:@include <netcdf.h>])])
for ax_netcdf_cflags in '' "$NETCDF_CFLAGS"; do
  AS_IF([test -z "$ax_netcdf_cflags"],
        [ax_netcdf_cflags_res="none required"],
        [ax_netcdf_cflags_res="$ax_netcdf_cflags" && CFLAGS="$ax_netcdf_cflags $_ax_c_lib_netcdf_save_CFLAGS"])
  AC_COMPILE_IFELSE([], [AS_VAR_SET([ax_netcdf_cflags_search], [$ax_netcdf_cflags_res])])
  AS_VAR_SET_IF([ax_netcdf_cflags_search], [break])
done
AS_VAR_SET_IF([ax_netcdf_cflags_search], [], [AS_VAR_SET([ax_netcdf_cflags_search], [no])])
rm conftest.$ac_ext
CFLAGS=$_ax_c_lib_netcdf_save_CFLAGS
AS_VAR_COPY([ax_netcdf_cflags_res], [ax_netcdf_cflags_search])
AS_IF([test "$ax_netcdf_cflags_res" != no],
      [test "$ax_netcdf_cflags_res" = "none required" || NETCDF_CFLAGS=$ax_netcdf_cflags_res
       ax_netcdf_cflags_search="yes"])
])

AS_IF([test "$ax_netcdf_cflags_search" = "no" || test "$ax_netcdf_lib_search" = "no" || test "$ax_netcdf_ldflags_search" = "no"],
[AC_MSG_WARN([Unable to locate a working C netCDF library.  Please specify --with-netcdf=<LOCATION> where <LOCATION> is the full path to where netCDF is installed.])
with_netcdf="no"],[
# Print the netCDF version, either from what we got with nc-config, or directly
# from the library
AC_MSG_CHECKING([the netCDF library version])
AS_IF([test -z "$_ax_c_lib_netcdf_version"], [
CFLAGS="$NETCDF_CFLAGS $_ax_c_lib_netcdf_save_CFLAGS"
LDFLAGS="$NETCDF_LDFLAGS $_ax_c_lib_netcdf_save_LDFLAGS"
LIBS="$NETCDF_LIBS $_ax_c_lib_netcdf_save_LIBS"
AC_RUN_IFELSE([AC_LANG_PROGRAM([@%:@include <netcdf.h>], [printf("%s",cn_inq_libvers());])],
              [NETCDF_VERSION=$(./conftest$EXEEXT)],
              [NETCDF_VESION=unknown])
CFLAGS=$_ax_c_lib_netcdf_save_CFLAGS
LDFLAGS=$_ax_c_lib_netcdf_save_LDFLAGS
LIBS=$_ax_c_lib_netcdf_save_LDFLAGS],
[NETCDF_VERSION=$_ax_c_lib_netcdf_version])
AC_MSG_RESULT([$NETCDF_VERSION])
AC_SUBST([NETCDF_VERSION])
AC_SUBST([NETCDF_CC])
AC_SUBST([NETCDF_CFLAGS])
AC_SUBST([NETCDF_LDFLAGS])
AC_SUBST([NETCDF_LIBS])
])
AS_VAR_POPDEF([ax_netcdf_lib_search])
AS_VAR_POPDEF([ax_netcdf_ldflags_search])
AS_VAR_POPDEF([ax_netcdf_cflags_search])dnl
])

#######################################################################
# NetCDF Fortran
AC_DEFUN([AX_LIB_NETCDFF], [
AC_LANG_ASSERT(Fortran)dnl
AC_REQUIRE([AC_PROG_SED])dnl
AC_REQUIRE([AC_PROG_AWK])dnl
AC_REQUIRE([AC_PROG_GREP])dnl

AC_ARG_WITH([netcdf-fortran],
[AS_HELP_STRING(
  [--with-netcdf-fortran=[yes/no/PATH]],
  [base directory of NetCDF Fortran installation])],
[],
[with_netcdf="yes"])
AS_IF([test $with_netcdf != no],
  [_AX_FC_LIB_NETCDF()])dnl
])

AC_DEFUN([_AX_FC_LIB_NETCDF],[
# Was a path to NetCDF given?
AS_IF([test $with_netcdf != yes], [_ax_fc_lib_netcdf_prefix=$with_netcdf])
test -z ${_ax_fc_lib_netcdf_prefix+x} && _ax_fc_lib_netcdf_bin=$_ax_fc_lib_netcdf_prefix/bin
# Check to see if NF_CONFIG is in the path, or NETCDF_PREFIX/bin.
AC_PATH_PROGS([NF_CONFIG], [nf-config], [], [$_ax_fc_lib_netcdf_bin$PATH_SEPARATOR$PATH])
# If nf-config found, use it to get information on the NetCDF library
AS_IF([test ! -z "${NF_CONFIG+x}"],[
  _ax_fc_lib_netcdf_fc=$(eval $NF_CONFIG --fc 2> /dev/null)
  _ax_fc_lib_netcdf_version=$(eval $NF_CONFIG --version 2> /dev/null)
  test -z $_ax_fc_lib_netcdf_prefix && _ax_fc_lib_netcdf_prefix=$(eval $NF_CONFIG --prefix 2> /dev/null)
  _ax_fc_lib_netcdf_fcflags=$(eval $NF_CONFIG --fflags 2> /dev/null)
  # Some systems (e.g. Cray) use a compiler wrapper that automatically add the required flag, which
  # has nf-config --fflags not return the -I flag.  This causes issues if trying to not use the
  # wrapper.  Thus, we add the -I flag if the cflags variable doesn't have it.  This should be
  # harmless.
  _ax_fc_lib_netcdf_includedir=$(eval $NF_CONFIG --includedir 2> /dev/null)
  echo "$_ax_fc_lib_netcdf_cflags" | $GREP -e "-I$_ax_fc_lib_netcdf_includedir" 2>&1 > /dev/null || _ax_fc_lib_netcdf_cflags="$_ax_fc_lib_netcdf_cflags -I$_ax_fc_lib_netcdf_includedir"
  # LIBS and LDFLAGS sorted based on prefix (e.g. -L and -l)
  for arg in $(eval $NF_CONFIG --flibs 2> /dev/null) ; do
    AS_CASE([$arg],
            [-L*], [echo $_ax_fc_lib_netcdf_ldflags | $GREP -e "$arg" 2>&1 >/dev/null || _ax_fc_lib_netcdf_ldflags="$arg $_ax_c_lib_netcdf_ldflags"],
            [-l*], [echo $_ax_fc_lib_netcdf_libs | $GREP -e "$arg" 2>&1 >/dev/null || _ax_fc_lib_netcdf_libs="$arg $_ax_c_lib_netcdf_libs"])
  done])
AS_IF([test ! -z ${_ax_fc_lib_netcdf_prefix+x}], [
  # Check if nc-config was able to give us specific information.  If not,
  # Guess the locations for libs and headers, and libs flags
  test -z "$_ax_fc_lib_netcdf_fcflags" && _ax_fc_lib_netcdf_fcflags=-I$_ax_fc_lib_netcdf_prefix/include
  test -z "$_ax_fc_lib_netcdf_ldflags" && _ax_fc_lib_netcdf_ldflags=-L$_ax_fc_lib_netcdf_prefix/lib
  test -z "$_ax_fc_lib_netcdf_libs" && _ax_fc_lib_netcdf_libs=-lnetcdff
])
AS_VAR_SET_IF([_ax_fc_lib_netcdf_fc],
              [AS_VAR_SET([NETCDF_FC],[$_ax_fc_lib_netcdf_fc])],
              [AS_VAR_SET([NETCDF_FC],[unknown])])
# Check for the netCDF Fortran library
# This is essentially the AC_SEARCH_LIBS, but as we need to loop through
# different libs and ldflags combinations to verify if -lnetcdff works.
AS_VAR_PUSHDEF([ax_netcdff_lib_search], [ax_cv_lib_netcdff])dnl
AS_VAR_PUSHDEF([ax_netcdff_ldflags_search], [ax_cv_ldflags_netcdff])
AC_CACHE_CHECK([for Fortran netCDF library], [ax_netcdff_lib_search], [
_ax_fc_lib_netcdf_save_LIBS=$LIBS
_ax_fc_lib_netcdf_save_LDFLAGS=$LDFLAGS
# Did the user set the NETCDF_LDFLAGS and/or NETCDF_LIBS?
test -z "$NETCDF_FCLDFLAGS" && NETCDF_FCLDFLAGS=$_ax_fc_lib_netcdf_ldflags
test -z "$NETCDF_FCLIBS" && NETCDF_FCLIBS=$_ax_fc_lib_netcdf_libs
AC_LANG_CONFTEST([AC_LANG_PROGRAM([], [      print*, nf_inq_libvers()])])
for ax_netcdff_ldflag in '' "$NETCDF_FCLDFLAGS"; do
  AS_IF([test -z "$ax_netcdff_ldflag"],
        [ax_netcdff_ldflags_res="none required" && LDFLAGS="$_ax_fc_lib_netcdf_save_LDFLAGS"],
        [ax_netcdff_ldflags_res="$ax_netcdff_ldflag" && LDFLAGS="$ax_netcdff_ldflag $_ax_fc_lib_netcdf_save_LDFLAGS"])
  for ax_netcdff_lib in '' "$NETCDF_FCLIBS" '-lnetcdff' "-lnetcdff $NETCDF_FCLIBS"; do
    AS_IF([test -z "$ax_netcdff_lib"],
          [ax_netcdff_lib_res="none required" && LIBS="$_ax_fc_lib_netcdf_save_LIBS"],
          [ax_netcdff_lib_res=$ax_netcdff_lib && LIBS="$ax_netcdff_lib $_ax_fc_lib_netcdf_save_LIBS"])
    AC_LINK_IFELSE([], [
      AS_VAR_SET([ax_netcdff_lib_search], [$ax_netcdff_lib_res])
      AS_VAR_SET([ax_netcdff_ldflags_search], [$ax_netcdff_ldflags_res])])
    AS_VAR_SET_IF([ax_netcdff_lib_search], [break 2])
  done
done
AS_VAR_SET_IF([ax_netcdff_lib_search], [], [AS_VAR_SET([ax_netcdff_lib_search], [no])])
AS_VAR_SET_IF([ax_netcdff_ldflags_search], [], [AS_VAR_SET([ax_netcdff_ldflags_search], [no])])
rm conftest.$ac_ext
LIBS=$_ax_fc_lib_netcdf_save_LIBS
LDFLAGS=$_ax_fc_lib_netcdf_save_LDFLAGS
AS_VAR_COPY([ax_netcdff_lib_res], [ax_netcdff_lib_search])
AS_VAR_COPY([ax_netcdff_ldflags_res], [ax_netcdff_ldflags_search])
AS_IF([test "$ax_netcdff_lib_res" != no && test "$ax_netcdff_ldflags_res" != no],
      [test "$ax_netcdff_lib_res" = "none required" || NETCDF_FCLIBS=$ax_netcdff_lib_res
       test "$ax_netcdff_ldflags_res" = "none required" || NETCDF_FCLDFLAGS=$ax_netcdff_ldflags_res
       ax_netcdff_lib_search="yes"])
])
AS_VAR_PUSHDEF([ax_netcdf_fcflags_search], [ax_cv_fcflags_netcdf])
AC_CACHE_CHECK([netcdf.inc usability], [ax_netcdf_fcflags_search], [
_ax_fc_lib_netcdf_save_FCFLAGS=$FCFLAGS
# Check if the user already set NETCDF_FCFLAGS
test -z "$NETCDF_FCFLAGS" && NETCDF_FCFLAGS=$_ax_fc_lib_netcdf_fcflags
AC_LANG_CONFTEST([AC_LANG_SOURCE([include 'netcdf.inc'])])
for ax_netcdf_fcflags in '' "$NETCDF_FCFLAGS"; do
  AS_IF([test -z "$ax_netcdf_fcflags"],
        [ax_netcdf_fcflags_res="none required"],
        [ax_netcdf_fcflags_res="$ax_netcdf_fcflags" && FCFLAGS="$ax_netcdf_fcflags $_ax_fc_lib_netcdf_save_FCFLAGS"])
  AC_COMPILE_IFELSE([], [AS_VAR_SET([ax_netcdf_fcflags_search], [$ax_netcdf_fcflags_res])])
  AS_VAR_SET_IF([ax_netcdf_fcflags_search], [break])
done
AS_VAR_SET_IF([ax_netcdf_fcflags_search], [], [AS_VAR_SET([ax_netcdf_fcflags_search], [no])])
rm conftest.$ac_ext
FCFLAGS=$_ax_fc_lib_netcdf_save_FCFLAGS
AS_VAR_COPY([ax_netcdf_cflags_res], [ax_netcdf_cflags_search])
AS_IF([test "$ax_netcdf_fcflags_res" != no],
      [test "$ax_netcdf_fcflags_res" = "none required" || NETCDF_FCFLAGS=$ax_netcdf_fcflags_res
       ax_netcdf_fcflags_search="yes"])
])
# Check for the netcdf module, and if it is usable with the current FC compiler
# Here we assume the netcdf module will be with the netcdf.inc file, which is
# usually a good assumption.
AS_VAR_PUSHDEF([ax_netcdf_mod_check], [ax_cv_netcdf_mod_available])
AC_CACHE_CHECK([netCDF Fortran module usability], [ax_netcdf_mod_check], [
_ax_fc_lib_netcdf_save_FCFLAGS=$FCFLAGS
FCFLAGS="$NETCDF_FCFLAGS $_ax_fc_lib_netcdf_save_FCFLAGS"
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [use netcdf])], [AS_VAR_SET([ax_netcdf_mod_check], [yes])], [AS_VAR_SET([ax_netcdf_mod_check], [no])])
FCFLAGS=$_ax_fc_lib_netcdf_save_FCFLAGS
])
# Check to see if the F90 bindings are usable in the netcdf library.  As we
# already know if the netCDF Fortran library is available, we can shortcircuit
# this test using the results above.
# Some compilers will not correctly mangle the module routine names unless the
# file extension is .f90.  We will switch the file extension for this test to
# .f90, and switch it back prior to returning.
_ax_fc_lib_netcdf_save_srcext=$ac_ext
AC_FC_SRCEXT([f90])
AS_VAR_PUSHDEF([ax_netcdf_lib_f90], [ax_cv_netcdf_lib_f90_available])
AC_CACHE_CHECK([netCDF Fortran 90 binding available], [ax_netcdf_lib_f90], [
AS_IF([test $ax_netcdff_lib_search = "no"],
      [ax_netcdf_lib_f90=no],[
  _ax_fc_lib_netcdf_save_FCFLAGS=$FCFLAGS
  _ax_fc_lib_netcdf_save_LDFLAGS=$LDFLAGS
  _ax_fc_lib_netcdf_save_LIBS=$LIBS
  FCFLAGS="$NETCDF_FCFLAGS $_ax_fc_lib_netcdf_save_FCFLAGS"
  LDFLAGS="$NETCDF_FCLDFLAGS $_ax_fc_lib_netcdf_save_LDFLAGS"
  LIBS="$NETCDF_FCLIBS $_ax_fc_lib_netcdf_save_LIBS"
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [      use netcdf
      print*, nf90_inq_libvers()])], [AS_VAR_SET([ax_netcdf_lib_f90], [yes])], [AS_VAR_SET([ax_netcdf_lib_f90], [no])])
  FCFLAGS=$_ax_fc_lib_netcdf_save_FCFLAGS
  LDFLAGS=$_ax_fc_lib_netcdf_save_LDFLAGS
  LIBS=$_ax_fc_lib_netcdf_save_LIBS])
])
ac_ext=$_ax_fc_lib_netcdf_save_srcext
AS_IF([test "$ax_netcdf_fcflags_search" = "no" || test "$ax_netcdf_mod_check" = "no" || test "$ax_netcdf_lib_f90" = "no" || \
  test "$ax_netcdff_lib_search" = "no" || test "$ax_netcdff_ldflags_search" = "no"],
[AC_MSG_WARN([Unable to locate a working Fortran netCDF library.  Please specify --with-netcdf-fortran=<LOCATION> where <LOCATION> is the full path to where netCDF is installed.])
with_netcdf_fortran="no"],[
# Print the netCDF version, either from what we got with nf-config, or directly
# from the library
AC_MSG_CHECKING([the Fortran netCDF library version])
AS_IF([test -z "$_ax_c_lib_netcdf_version"], [
FCFLAGS="$NETCDF_FCFLAGS $_ax_fc_lib_netcdf_save_FCFLAGS"
LDFLAGS="$NETCDF_FCLDFLAGS $_ax_fc_lib_netcdf_save_LDFLAGS"
LIBS="$NETCDF_FCLIBS $_ax_fc_lib_netcdf_save_LIBS"
AC_RUN_IFELSE([AC_LANG_PROGRAM([], [      include 'netcdf.inc'
      print*, nf_inq_libvers()])],
              [NETCDFF_VERSION=$(./conftest$EXEEXT)],
              [NETCDFF_VESION=unknown])
FCFLAGS=$_ax_fc_lib_netcdf_save_FCFLAGS
LDFLAGS=$_ax_fc_lib_netcdf_save_LDFLAGS
LIBS=$_ax_fc_lib_netcdf_save_LDFLAGS],
[NETCDFF_VERSION=$_ax_fc_lib_netcdf_version])
AC_MSG_RESULT([$NETCDFF_VERSION])
AC_SUBST([NETCDFF_VERSION])
AC_SUBST([NETCDF_FC])
AC_SUBST([NETCDF_FCFLAGS])
AC_SUBST([NETCDF_FCLDFLAGS])
AC_SUBST([NETCDF_FCLIBS])
AC_DEFINE([HAVE_NETCDFF], [1], [Defined if you have netCDF Fortran support])
])
AS_VAR_POPDEF([ax_netcdf_fcflags_search])
AS_VAR_POPDEF([ax_netcdf_mod_check])
AS_VAR_POPDEF([ax_netcdf_lib_f90])
AS_VAR_POPDEF([ax_netcdff_lib_search])
AS_VAR_POPDEF([ax_netcdff_ldflags_search])
])
