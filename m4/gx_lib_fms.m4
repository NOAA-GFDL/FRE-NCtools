#cspell:ignore LIBFMS xyes
# ===========================================================================
#
# SYNOPSIS
#
#   GX_LIB_FMS()
#
# DESCRIPTION
#
#   This macro checks for the presence of the FMS library.  FMS is a Fortran
#   library.  Thus, this macro runs AC_LANG_ASSERT([Fortran]) to ensure that
#   the Fortran compiler has been set.
#
#   If the FMS library is found, this macro calls
#
#      AC_SUBST([LIBFMS_FCFLAGS])
#      AC_SUBST([LIBFMS_LIBS])
#
#   It also sets
#
#      with_fms="yes"
#
#   If FMS is disabled or not found, this macro sets
#
#      with_fms="no"
#
#   Your configuration script can test $with_fms to take any further
#   actions.
#
#   To use the macro, do the following in "configure.ac" before
#   AC_OUTPUT:
#
#      GX_LIB_FMS()
#
# LICENSE
#
#   Copyright (c) 2025 Seth Underwood <Seth.Underwood@noaa.gov>
#
#   Copying and distribution of this file, with or without
#   modification, are permitted in any medium without royalty provided
#   the copyright notice and this notice are preserved. This file is
#   offered as-is, without any warranty.

AC_DEFUN([GX_LIB_FMS], [
AC_LANG_ASSERT([Fortran])dnl
AS_IF([test "x$1" = "required"],
[AS_VAR_SET([default_with_fms], [yes])],
[AS_VAR_SET([default_with_fms], [no])])

AC_ARG_WITH([fms],
[AS_HELP_STRING([--with-fms=[yes/no/PATH]],
  [Use an external libFMS library, at optional PATH location (default no)])],
  [],
  [with_fms="$default_with_fms"])
AS_IF([test "x$with_fms" != "xno"],
[_GX_TEST_LIB_FMS()])dnl
])

AC_DEFUN([_GX_TEST_LIB_FMS], [
# Check if the user supplied a path to the FMS library
AS_IF([test "x$with_fms" != "xyes"], [_gx_lib_fms_prefix="$with_fms"])
# Set the lib and include paths from the prefix
AS_IF([test ! -z ${_gx_lib_fms_prefix+x}], [
AS_VAR_SET([_gx_lib_fms_lib], [${_gx_lib_fms_prefix}/lib])
AS_VAR_SET([_gx_lib_fms_include], [${_gx_lib_fms_prefix}/include])])

# Some compilers do not correctly mangle the module routine names
# unless the Fortran file extension is .f90.  This is a workaround for
# that issue.  The Fortran file extension will be reverted prior to
# returning.
_gx_fc_lib_fms_save_srcext=$ac_ext
AC_FC_SRCEXT([f90])

# Check if the FMS Fortran module, and if it usable with the current FC compiler.
AC_CACHE_CHECK([FMS module usability], [gx_cv_fms_mod_check],[
_gx_fc_lib_fms_save_FCFLAGS="$FCFLAGS"
FCFLAGS="$FCFLAGS -I$_gx_lib_fms_include"
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [use fms])],
[AS_VAR_SET([gx_cv_fms_mod_check], [yes])],
[AS_VAR_SET([gx_cv_fms_mod_check], [no])])
FCFLAGS="$_gx_fc_lib_fms_save_FCFLAGS"
])
# Check if the FMS library is usable with the current FC compiler.
AC_CACHE_CHECK([FMS library is available], [gx_cv_fms_lib_check],[
_gx_fc_lib_fms_save_FCFLAGS="$FCFLAGS"
_gx_fc_lib_fms_save_LDFLAGS="$LDFLAGS"
_gx_fc_lib_fms_save_LIBS="$LIBS"
FCFLAGS="$FCFLAGS -I$_gx_lib_fms_include"
LDFLAGS="$LDFLAGS -L$_gx_lib_fms_lib"
LIBS="$LIBS -lFMS"
AC_LINK_IFELSE([AC_LANG_PROGRAM([], [use fms])],
[AS_VAR_SET([gx_cv_fms_lib_check], [yes])],
[AS_VAR_SET([gx_cv_fms_lib_check], [no])])
FCFLAGS="$_gx_fc_lib_fms_save_FCFLAGS"
LDFLAGS="$_gx_fc_lib_fms_save_LDFLAGS"
LIBS="$_gx_fc_lib_fms_save_LIBS"
])
# Reset the Fortran file extension
ac_ext=$_gx_fc_lib_fms_save_srcext
# Warn if the user requested FMS, but we were unable to find a working
# library
AS_IF([test "$gx_cv_fms_mod_check" = "no" -o "$gx_cv_fms_lib_check" = "no"],
[AC_MSG_WARN([Unable to find a working FMS library.  Please specify --with-fms=<LOCATION> where location is the path to the FMS library.])
with_fms="no"],[
# Return the variables
with_fms="yes"
AC_SUBST([LIBFMS_FCFLAGS], ["-I$_gx_lib_fms_include"])
AC_SUBST([LIBFMS_LIBS], ["-L$_gx_lib_fms_lib -lFMS"])
])
])
