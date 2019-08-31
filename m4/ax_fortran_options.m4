# AX_REAL_KIND8_FLAG([ACTION-IF-SUCCESS], [ACTION-IF-FAILURE = FAILURE])
# ----------------------------------------------------------------------
# Look for the compiler flag that sets the default REAL kind to KIND=8.
# Call ACTION-IF-SUCCESS (defaults to nothing) if successful (i.e. can
# compile with default REAL(KIND=8)) and ACTION-IF-FAILURE (defaults
# to failing with an error message) if not.
#
# Sets the variable FC_REAL_KIND8_FLAG to hold the flag.
#
# The known flags are:
# -fdefault-real-8: gfortran
#    -real_size 64: Intel compiler
#        -s real64: Cray
#              -r8: Portland Group compiler
#     -qrealsize=8: IBM compiler
AC_DEFUN([AX_REAL_KIND8_FLAG],[
AC_LANG_PUSH([Fortran])
AC_CACHE_CHECK([for Fortran default REAL KIND 8 flag], [ax_cv_fc_real_kind8_flag],[
ax_cv_fc_real_kind8_flag=unknown
ax_real_kind8_flag_FCFLAGS_save=$FCFLAGS
for ac_flag in none \
               '-fdefault-real-8' \
               '-real_size 64' \
               '-s real64' \
               '-r8' \
               '-qrealsize=8'; do
  test "x$ac_flag" != xnone && FCFLAGS="$ax_real_kind8_flag_FCFLAGS_save ${ac_flag}"
  AC_COMPILE_IFELSE([[
     program test
     interface
     subroutine test_sub(a)
     real(kind=8) :: a
     end subroutine test_sub
     end interface
     real :: b=1.0
     call test_sub(b)
     end program test]],
     [ax_cv_fc_real_kind8_flag=$ac_flag; break])
done
rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
FCFLAGS=$ax_real_kind8_flag_FCFLAGS_save
])
if test "x$ax_cv_fc_real_kind8_flag" = xunknown; then
  m4_default([$3],
              [AC_MSG_ERROR([Fortran cannot set default real kind to 8])])
else
  FC_REAL_KIND8_FLAG=$ax_cv_fc_real_kind8_flag
  if test "x$FC_REAL_KIND8_FLAG" = xnone; then
    FC_REAL_KIND8_FLAG=
  fi
  $2
fi
AC_LANG_POP([Fortran])
AC_SUBST([FC_REAL_KIND8_FLAG])
])

# AX_REAL_KIND4_FLAG([ACTION-IF-SUCCESS], [ACTION-IF-FAILURE = FAILURE])
# ----------------------------------------------------------------------
# Look for the compiler flag that sets the default REAL kind to KIND=4.
# Call ACTION-IF-SUCCESS (defaults to nothing) if successful (i.e. can
# compile with default REAL(KIND=4)) and ACTION-IF-FAILURE (defaults
# to failing with an error message) if not.
#
# Sets the variable FC_REAL_KIND4_FLAG to hold the flag.
#
# The known flags are:
#             none: gfortran (gfortran does not have an option to set the
#                   default REAL kind to KIND=4)
#    -real_size 32: Intel compiler
#        -s real32: Cray
#              -r4: Portland Group compiler
#     -qrealsize=4: IBM compiler
AC_DEFUN([AX_REAL_KIND4_FLAG],[
AC_LANG_PUSH([Fortran])
AC_CACHE_CHECK([for Fortran default REAL KIND 4 flag], [ax_cv_fc_real_kind4_flag],[
ax_cv_fc_real_kind4_flag=unknown
ax_real_kind4_flag_FCFLAGS_save=$FCFLAGS
for ac_flag in none \
               '-fdefault-real-4' \
               '-real_size 32' \
               '-s real32' \
               '-r4' \
               '-qrealsize=4'; do
  test "x$ac_flag" != xnone && FCFLAGS="$ax_real_kind4_flag_FCFLAGS_save ${ac_flag}"
  AC_COMPILE_IFELSE([[
     program test
     interface
     subroutine test_sub(a)
     real(kind=4) :: a
     end subroutine test_sub
     end interface
     real :: b=1.0
     call test_sub(b)
     end program test]],
     [ax_cv_fc_real_kind4_flag=$ac_flag; break])
done
rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
FCFLAGS=$ax_real_kind4_flag_FCFLAGS_save
])
if test "x$ax_cv_fc_real_kind4_flag" = xunknown; then
  m4_default([$3],
              [AC_MSG_ERROR([Fortran cannot set default real kind to 4])])
else
  FC_REAL_KIND4_FLAG=$ax_cv_fc_real_kind4_flag
  if test "x$FC_REAL_KIND4_FLAG" = xnone; then
    FC_REAL_KIND4_FLAG=
  fi
  $2
fi
AC_LANG_POP([Fortran])
AC_SUBST([FC_REAL_KIND4_FLAG])
])

# AX_CRAY_POINTER_FLAG([ACTION-IF-SUCCESS], [ACTION-IF-FAILURE = FAILURE])
# -----------------------------------------------------------------------------
# Look for the compiler flag that allows Fortran Cray Pointers.  Cray
# pointers are an are part of a non-standard extension that provides a
# C-like pointer in Fortran.  Call ACTION-IF-SUCCESS (defaults to
# nothing) if successful (i.e. can use Cray pointers) and
# ACTION-IF-FAILURE (defaults to failing with an error message) if not.
#
# Sets the variable FC_CRAY_POINTER_FLAG to hold the flag.
#
# The known flags are:
# -fcray-pointer: gfortran
#           none: Intel compiler (No option required for Cray Pointers)
#        unknown: Cray
# -Mcray=pointer: Portland Group compiler
#           none: IBM compiler (No option required for Cray Pointers)
AC_DEFUN([AX_CRAY_POINTER_FLAG],[
AC_LANG_PUSH([Fortran])
AC_CACHE_CHECK([for Fortran Cray Pointer flag], [ax_cv_fc_cray_ptr_flag],[
ax_cv_fc_cray_ptr_flag=unknown
ax_cray_ptr_flag_FCFLAGS_save=$FCFLAGS
for ac_flag in none \
               '-fcray-pointer' \
               '-Mcray=pointer'; do
  test "x$ac_flag" != xnone && FCFLAGS="$ax_cray_ptr_flag_FCFLAGS_save ${ac_flag}"
  AC_COMPILE_IFELSE([[
     program test
     integer(kind=8) :: ipt
     integer iarri(10)
     pointer (ipt, iarr)
     end program test]],
     [ax_cv_fc_cray_ptr_flag=$ac_flag; break])
done
rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
FCFLAGS=$ax_cray_ptr_flag_FCFLAGS_save
])
if test "x$ax_cv_fc_cray_ptr_flag" = xunknown; then
  m4_default([$3],
              [AC_MSG_ERROR([Fortran cannot use Cray pointers])])
else
  FC_CRAY_POINTER_FLAG=$ax_cv_fc_cray_ptr_flag
  if test "x$FC_CRAY_POINTER_FLAG" = xnone; then
    FC_CRAY_POINTER_FLAG=
  fi
  $2
fi
AC_LANG_POP([Fortran])
AC_SUBST([FC_CRAY_POINTER_FLAG])
])
