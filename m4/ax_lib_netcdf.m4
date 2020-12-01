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
AS_IF([test ! -z ${NC_CONFIG+x}],[
  _ax_c_lib_netcdf_cc=$(eval $NC_CONFIG --cc | $AWK '{print $[]1}')
  AS_IF([test "x$_ax_c_lib_netcdf_cc" = "xccache"],
        [_ax_c_lib_netcdf_cc=$(eval $NC_CONFIG --cc | $ASK '{print $[]1}')])
  _ax_c_lib_netcdf_version=$(eval $NC_CONFIG --version | $AWK '{print $[]2}')
  test -z $_ax_c_lib_netcdf_prefix && _ax_c_lib_netcdf_prefix=$(eval $NC_CONFIG --prefix)
  _ax_c_lib_netcdf_cflags=$(eval $NC_CONFIG --cflags)
  # LIBS and LDFLAGS sorted based on prefix (e.g. -L and -l)
  for arg in $(eval $NC_CONFIG --libs) ; do :
    AS_CASE([$arg],
            [-L*], [echo $_ax_c_lib_netcdf_ldflags | $GREP -e "$arg" 2>&1 >/dev/null \
                    || _ax_c_lib_netcdf_ldflags="$arg $_ax_c_lib_netcdf_ldflags"],
            [-l*], [echo $_ax_c_lib_netcdf_libs | $GREP -e "$arg" 2>&1 >/dev/null \
                    || _ax_c_lib_netcdf_libs="$arg $_ax_c_lib_netcdf_libs"])
  done],
  [test ! -z ${_ax_c_lib_netcdf_prefix+x}], [
  # nc-config not found, and we were given a location of NetCDF.
  # Guess the locations for libs and headers locations
  _ax_c_lib_netcdf_cflags=-I$_ax_c_lib_netcdf_prefix/include
  _ax_c_lib_netcdf_ldflags=-L$_ax_c_lib_netcdf_prefix/lib
  _ax_c_lib_netcdf_libs=-lnetcdf
])

# Check for nc_create in netcdf
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
AC_LANG_CONFTEST([AC_LANG_CALL([], [nc_inq_libver])])
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
       test "$ax_netcdf_ldflags_res" = "none required" || NETCDF_LDFLAGS=$ax_netcdf_ldflags_res])
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
        [ax_netcdf_cflags_res="non required"],
        [ax_netcdf_cflags_res="$ax_netcdf_cflags" && CFLAGS="$ax_netcdf_cflags $_ax_c_lib_netcdf_save_CFLAGS"])
  AC_COMPILE_IFELSE([], [AS_VAR_SET([ax_netcdf_cflags_search], [$ax_netcdf_cflags_res])])
  AS_VAR_SET_IF([ax_netcdf_cflags_search], [break])
done
AS_VAR_SET_IF([ax_netcdf_cflags_search], [], [AS_VAR_SET([ax_netcdf_cflags_search], [no])])
rm conftest.$ac_ext
CFLAGS=$_ax_c_lib_netcdf_save_CFLAGS
AS_VAR_COPY([ax_netcdf_cflags_res], [ax_netcdf_cflags_search])
AS_IF([test "$ax_netcdf_cflags_res" != no],
      [test "$ax_netcdf_lib_res" = "none required" || NETCDF_CFLAGS=$ax_netcdf_cflags_res
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

AC_REQUIRE([AC_PROG_SED])
AC_REQUIRE([AC_PROG_AWK])
AC_REQUIRE([AC_PROG_GREP])

# Add a default --with-netcdf-fortran configuration option.
AC_ARG_WITH([netcdf-fortran],
  [AS_HELP_STRING(
    [--with-netcdf-fortran=[yes/no/PATH]],
    [base directory of NetCDF Fortran installation])],
  [if test "$withval" = "no"; then
     with_netcdf_fortran="no"
   elif test "$withval" = "yes"; then
     with_netcdf_fortran="yes"
   else
     with_netcdf_fortran="yes"
     NETCDFF_PREFIX="${withval}"
     NF_CONFIG="${withval}/bin/nf-config"
   fi],
   [with_netcdf_fortran="yes"]
)

# Try and find NetCDF tools and options.
if test "$with_netcdf_fortran" = "yes"; then
  if test -z "$NF_CONFIG"; then
    # Check to see if NF_CONFIG is in the path.
    AC_PATH_PROGS([NF_CONFIG], [nf-config], [])
    NETCDFF_PREFIX=$(AS_DIRNAME([$(AS_DIRNAME(["$NF_CONFIG"]))]))
  else
    AC_MSG_CHECKING([Using provided NetCDF Fortran prefix])
    AC_MSG_RESULT([$NF_CONFIG])
  fi

  AC_MSG_CHECKING([for NetCDF Fortran libraries])

  if test ! -f "$NF_CONFIG" || test ! -x "$NF_CONFIG"; then
    AC_MSG_RESULT([no])
    AC_MSG_WARN([

Unable to locate NetCDF Fortran compilation helper script 'nf-config'.
Please specify --with-netcdf-fortran=<LOCATION> as the full path prefix
where NetCDF Fortran has been installed.
])
    with_netcdf_fortran="no"
  else
    # Get the actual compiler used
    NETCDF_FC=$(eval $NF_CONFIG --fc | $AWK '{print $[]1}')

    # Look for version
    NETCDFF_VERSION=$(eval $NF_CONFIG --version | $AWK '{print $[]2}')

    # Look for the CFLAGS
    test -z "$NETCDF_FCCPPFLAGS" && NETCDF_FCCPPFLAGS=$(eval $NF_CONFIG --cflags)

    # Look for the FCFLAGS
    test -z "$NETCDF_FCFLAGS" && NETCDF_FCFLAGS=$(eval $NF_CONFIG --fflags)

    # Look for the LIBS and LDFLAGS
    NETCDF_tmp_fclibs=$(eval $NF_CONFIG --flibs)

    # Sort out the tmp libs based on their prefixes
    for arg in $NETCDF_tmp_fclibs ; do
      case "$arg" in
        -L*) echo $NETCDF_FCLDFLAGS_tmp | $GREP -e "$arg" 2>&1 >/dev/null \
                || NETCDF_FCLDFLAGS_tmp="$arg $NETCDF_FCLDFLAGS_tmp"
              ;;
        -l*) echo $NETCDF_FCLIBS_tmp | $GREP -e "$arg" 2>&1 >/dev/null \
                || NETCDF_FCLIBS_tmp="$arg $NETCDF_FCLIBS_tmp"
              ;;
      esac
    done
    test -z "$NETCDF_FCLDFLAGS" && NETCDF_FCLDFLAGS="$NETCDF_FCLDFLAGS_tmp"
    test -z "$NETCDF_FCLIBS" && NETCDF_FCLIBS="$NETCDF_FCLIBS_tmp"

    AC_MSG_RESULT([yes (version $[NETCDFF_VERSION])])

    # Save previous FC settings
    #ax_lib_netcdf_save_FC=$FC
    ax_lib_netcdf_save_FCFLAGS=$FCFLAGS
    ax_lib_netcdf_save_LIBS=$LIBS
    ax_lib_netcdf_save_LDFLAGS=$LDFLAGS

    #FC=$NETCDF_FC
    FCFLAGS=$NETCDF_FCFLAGS
    LIBS=$NETCDF_FCLIBS
    LDFLAGS=$NETCDF_FCLDFLAGS

    AC_LANG_PUSH([Fortran])
    AC_FC_SRCEXT([F90])

    ## Checking presence of netcdf.inc
    AC_MSG_CHECKING(for netcdf.inc is available)
    AC_COMPILE_IFELSE(
      AC_LANG_PROGRAM([],
        [#include <netcdf.inc>]),
        [ac_cv_netcdff_inc=yes],
        [ac_cv_netcdff_inc=no])
    AC_MSG_RESULT($ac_cv_netcdff_inc)

    ## Checking presence of netcdf Fortran 90 bindings
    AC_MSG_CHECKING(if F90 netCDF bindings are available)
    AC_COMPILE_IFELSE(
      AC_LANG_PROGRAM([],
        [use netcdf]),
        [ac_cv_netcdff_mod=yes],
        [ac_cv_netcdff_mod=no])
    AC_MSG_RESULT($ac_cv_netcdff_mod)

    ## Check fo use of library
    AC_MSG_CHECKING(if F90 netCDF library usable)
    AC_LINK_IFELSE(
      AC_LANG_PROGRAM([],
        [[use netcdf
        implicit none
        integer :: ncid, status
        status = nf90_open(path = "foo.nc", mode = nf90_nowrite, ncid = ncid)
        ]]),
        [ac_cv_libnetcdff=yes],
        [ac_cv_libnetcdff=no])
    AC_MSG_RESULT($ac_cv_libnetcdff)

    AC_LANG_POP([Fortran])

    if test "$ac_cv_netcdf_inc" = "no" || \
       test "$ac_cv_netcdff_mod" = "no" || \
       test "$ac_cv_libnetcdff" = "no" ; then
         AC_MSG_WARN([Unable to compile NetCDF Fortran test program])
         with_netcdf_fortran="no"
    fi

    # Restore previous FC settings
    #FC=$ax_lib_netcdf_save_FC
    FCFLAGS=$ax_lib_netcdf_save_FCFLAGS
    LIBS=$ax_lib_netcdf_save_LIBS
    LDFLAGS=$ax_lib_netcdf_save_LDFLAGS

    AC_SUBST([NETCDFF_VERSION])
    AC_SUBST([NETCDF_FC])
    AC_SUBST([NETCDF_FCFLAGS])
    AC_SUBST([NETCDF_FCLDFLAGS])
    AC_SUBST([NETCDF_FCLIBS])
    AC_DEFINE([HAVE_NETCDFF], [1], [Defined if you have netCDF Fortran support])
  fi
fi
])
