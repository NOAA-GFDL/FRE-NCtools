# ===========================================================================
#
# SYNOPSIS
#
#   GX_OPENACC_FLAGS()
#
# DESCRIPTION
#
#  Checks C compiler flags for openacc support.
#
# LICENSE
#
#   Copyright (c) 2024 Ryan Mulhall
#
#   This program is free software; you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation; either version 3 of the License, or (at your
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

# ----------------------------------------------------------------------
#
#  Will set OPENACC_CFLAGS to openacc flags for a given compiler if accepted.
#
#  Mainly for nvhpc, offloading with gcc is not currently supported,
#  although the build will still work.
#
#
#  Will set OPENACC_CFLAGS to OpenACC flags for a given compiler if accepted.
#
AC_DEFUN([GX_OPENACC_FLAGS],[
AC_CACHE_CHECK([for C compiler OpenACC flags], [gx_cv_openacc_flags], [

  AC_LANG_PUSH([C])
  gx_cv_openacc_flags=unknown
  gx_openacc_flags_CFLAGS_save=$CFLAGS

  # Try common OpenACC flags (-acc for NVHPC/PGI, -fopenacc for GCC/Clang)
  for ac_flag in '-acc' '-fopenacc'; do
    CFLAGS="$gx_openacc_flags_CFLAGS_save $ac_flag"
    AC_LINK_IFELSE([AC_LANG_PROGRAM(
      [[
#ifdef HAVE_OPENACC_H
# include <openacc.h>
#else
  extern int acc_get_device_type(void);
#endif
      ]],
      [[
        int dev = acc_get_device_type();
        (void)dev;
      ]])],
      [gx_cv_openacc_flags="$ac_flag"; break])
  done

  CFLAGS="$gx_openacc_flags_CFLAGS_save"
  AC_LANG_POP([C])
])

if test "x$gx_cv_openacc_flags" = "xunknown"; then
  OPENACC_CFLAGS=""
  m4_default([$2], [AC_MSG_WARN([OpenACC flags not found or not supported by compiler])])
else
  OPENACC_CFLAGS="$gx_cv_openacc_flags"
fi

AC_SUBST([OPENACC_CFLAGS])
])
