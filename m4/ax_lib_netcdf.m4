# ===========================================================================
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
#     AC_SUBST(NETCDF_FFLAGS)
#     AC_SUBST(NETCDF_FCCPPFLAGS)
#     AC_SUBST(NETCDF_FLIBS)
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

AC_REQUIRE([AC_PROG_SED])
AC_REQUIRE([AC_PROG_AWK])
AC_REQUIRE([AC_PROG_GREP])

# Add a default --with-netcdf configuration option.
AC_ARG_WITH([netcdf],
  [AS_HELP_STRING(
    [--with-netcdf=[yes/no/PATH]],
    [base directory of NetCDF installation])],
  [if test "$withval" = "no"; then
     with_netcdf="no"
   elif test "$withval" = "yes"; then
     with_netcdf="yes"
   else
     with_netcdf="yes"
     NETCDF_PREFIX="${withval}"
     NC_CONFIG="${withval}/bin/nc-config"
   fi],
   [with_netcdf="yes"]
)

# Try and find NetCDF tools and options.
if test "$with_netcdf" = "yes"; then
  if test -z "$NC_CONFIG"; then
    # Check to see if NC_CONFIG is in the path.
    AC_PATH_PROGS([NC_CONFIG], [nc-config], [])
    NETCDF_PREFIX=$(AS_DIRNAME([$(AS_DIRNAME(["$NC_CONFIG"]))]))
  else
    AC_MSG_CHECKING([Using provided NetCDF prefix])
    AC_MSG_RESULT([$NC_CONFIG])
  fi

  AC_MSG_CHECKING([for NetCDF libraries])

  if test ! -f "$NC_CONFIG" || test ! -x "$NC_CONFIG"; then
    AC_MSG_RESULT([no])
    AC_MSG_WARN([

Unable to locate NetCDF compilation helper script 'nc-config'.
Please specify --with-netcdf=<LOCATION> as the full path prefix
where NetCDF has been installed.
])
    with_netcdf="no"
  else
    # Get the actual compiler used
    NETCDF_CC=$(eval $NC_CONFIG --cc | $AWK '{print $[]1}')
    if test "$NETCDF_CC" = "ccache"; then
      NETCDF_CC=$(eval $NC_CONFIG --cc | $AWK '{print $[]2}')
    fi

    # Look for version
    NETCDF_VERSION=$(eval $NC_CONFIG --version | $AWK '{print $[]2}')

    # Look for the CFLAGS
    test -z "$NETCDF_CFLAGS" && NETCDF_CFLAGS=$(eval $NC_CONFIG --cflags)

    # Look for the CPPFLAGS
    test -z "$NETCDF_CPPFLAGS" && NETCDF_CPPFLAGS=$(eval $NC_CONFIG --cflags)

    # Look for the LIBS and LDFLAGS
    NETCDF_tmp_clibs=$(eval $NC_CONFIG --libs)

    # Sort out the tmp libs based on their prefixes
    for arg in $NETCDF_tmp_clibs ; do
      case "$arg" in
        -L*) echo $NETCDF_LDFLAGS_tmp | $GREP -e "$arg" 2>&1 >/dev/null \
                || NETCDF_LDFLAGS_tmp="$arg $NETCDF_LDFLAGS_tmp"
              ;;
        -l*) echo $NETCDF_LIBS_tmp | $GREP -e "$arg" 2>&1 >/dev/null \
                || NETCDF_LIBS_tmp="$arg $NETCDF_LIBS_tmp"
              ;;
      esac
    done

    test -z "$NETCDF_LDFLAGS" && NETCDF_LDFLAGS=$NETCDF_LDFLAGS_tmp
    test -z "$NETCDF_LIBS" && NETCDF_LIBS=$NETCDF_LIBS_tmp

    AC_MSG_RESULT([yes (version $[NETCDF_VERSION])])

    # Save current CC settings
    #ax_lib_netcdf_save_CC=$CC
    ax_lib_netcdf_save_CFLAGS=$CFLAGS
    ax_lib_netcdf_save_CPPFLAGS=$CPPFLAGS
    ax_lib_netcdf_save_LIBS=$LIBS
    ax_lib_netcdf_save_LDFLAGS=$LDFLAGS

    #CC=$NETCDF_CC
    CFLAGS=$NETCDF_CFLAGS
    CPPFLAGS=$NETCDF_CPPFLAGS
    LIBS=$NETCDF_LIBS
    LDFLAGS=$NETCDF_LDFLAGS

    AC_CHECK_HEADER([netcdf.h], [ac_cv_netcdf_h=yes], [ac_cv_netcdf_h=no])
    AC_CHECK_LIB([netcdf], [nc_create], [ac_cv_libnetcdf=yes],[ac_cv_libnetcdf=no])

    if test "$ac_cv_netcdf_h" = "no" || \
       test "$ac_cv_libnetcdf" = "no" ; then
         AC_MSG_WARN([Unable to compile NetCDF c test program])
         with_netcdf="no"
    fi

    #CC=$ax_lib_netcdf_save_CC
    CFLAGS=$ax_lib_netcdf_save_CFLAGS
    CPPFLAGS=$ax_lib_netcdf_save_CPPFLAGS
    LIBS=$ax_lib_netcdf_save_LIBS
    LDFLAGS=$ax_lib_netcdf_save_LDFLAGS

    AC_SUBST([NETCDF_VERSION])
    AC_SUBST([NETCDF_CC])
    AC_SUBST([NETCDF_CFLAGS])
    AC_SUBST([NETCDF_LDFLAGS])
    AC_SUBST([NETCDF_LIBS])
    AC_DEFINE([HAVE_NETCDF], [1], [Defined if you have netCDF support])
  fi
fi
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
