# ===========================================================================
#
# SYNOPSIS
#
#   GX_CC_OPENACC_FLAGS()
#
# DESCRIPTION
#
#  Checks C compiler flags for openacc support.
#
# LICENSE
#
#   Copyright (c) 2023 Ryan Mulhall
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
#  Will set CC_OPENACC_FLAGS to openacc flags for a given compiler if accepted.
#  Will also add an accepted openacc flag to LDFLAGS
#
#  First checks for general openacc support flag, then optimizations and targets.
#
#  Mainly for nvhpc, gcc support for openacc is limited so some (seq) directives
#  are ifdef'd out when using gcc.
#
AC_DEFUN([GX_CC_OPENACC_FLAGS],[
AC_CACHE_CHECK([whether C compiler accepts OpenACC flags], [gx_cv_cc_openacc_flags],[

AC_LANG_ASSERT(C)
gx_cv_cc_openacc_flags=unknown
gx_cc_openacc_flags_CFLAGS_save=$CFLAGS

dnl check for base openacc flag
for ac_flag in '-acc' \
               '-fopenacc'; do
  AC_LINK_IFELSE([AC_LANG_SOURCE(
          extern int acc_get_device_type();
          int main(int argc, char** argv){
              acc_get_device_type();
              return 0;
          })],
     [gx_cv_cc_openacc_flags="$gx_cc_openacc_flags_CFLAGS_save ${ac_flag}"]; break)
done
rm -f conftest.err conftest.$ac_objext conftest.$ac_ext

LDFLAGS="${ac_flag} ${LDFLAGS}"

dnl check for optimization, target and any other compiler-specific flags to use for gpu offloading
dnl TODO gcc flags for the offload, -foffload=nvptx-none could work if configured properly
for extra_flags in '-O2 -tp native -gpu=ccnative -Minfo=accel -Mnoinline' \
                   '-O2'; do
  test "x$extra_flags" != xnone && CFLAGS="${ac_flag} ${extra_flags}"
  AC_LINK_IFELSE([AC_LANG_SOURCE(
          extern int acc_get_device_type();
          int main(int argc, char** argv){
              acc_get_device_type();
              return 0;
          })],
     [gx_cv_cc_openacc_flags="$ac_flag $extra_flags"]; break)
done
rm -f conftest.err conftest.$ac_objext conftest.$ac_ext

CFLAGS="$gx_cc_openacc_flags_CFLAGS_save"

if test "x$gx_cv_cc_openacc_flags" = xunknown; then
  m4_default([$2],
              [AC_MSG_ERROR([no])])
else
  CC_OPENACC_FLAGS="${gx_cv_cc_openacc_flags}"
  AC_MSG_RESULT([yes])
fi
],
AC_SUBST([CC_OPENACC_FLAGS])
)])
