# Prelude
AC_INIT([FRE NC Tools], [1.0], [oar.gfdl.help@noaa.gov], [fre-nctools], [https://github.com/NOAA-GFDL/FRE-NCtools])
AC_PREREQ([2.63])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])

AC_CONFIG_SRCDIR([postprocessing/land_utils/combine-ncc.F90])
AC_CONFIG_MACRO_DIR([m4])

# Check for programs
AC_PROG_CC([icc gcc])
AM_PROG_CC_C_O

AC_PROG_FC([ifort gfortran])
AC_PROG_FC_C_O
AC_FC_SRCEXT([f90])
AC_FC_SRCEXT([F90])

# Check for library
AX_LIB_NETCDF()
if test "$with_netcdf" = "no"; then
  AC_MSG_ERROR([Unable to find NetCDF c library.])
fi

AX_LIB_NETCDFF()
if test "$with_netcdf_fortran" = "no"; then
  AC_MSG_ERROR([Unable to find NetCDF Fortran library.])
fi

# Check for header files

# Check for typedefs, structures, and compiler characteristics

# Check for library functions
AC_FC_FREEFORM()

# Set compiler family specific compiler options
if test $CC = gcc; then
  BASE_CFLAGS=
fi
if test $CC = icc; then
  BASE_CFLAGS="-traceback"
fi
if test $FC = gfortran; then
  BASE_FCFLAGS="-ffree-line-length-none"
  R8_FCFLAGS="-fdefault-real-8"
  R4_FCFLAGS=""
  MODDIR_FCFLAGS="-J"
fi
if test $FC = ifort; then
  BASE_FCFLAGS="-fltconsistency -fno-alias -stack_temps -safe_cray_ptr -ftz -assume byterecl -i4 -real_size 64 -traceback"
  R8_FCFLAGS="-real_size 64"
  R4_FCFLAGS="-real_size 32"
  MODDIR_FCFLAGS="-module"
fi

AC_SUBST([BASE_CFLAGS])
AC_SUBST([BASE_FCFLAGS])
AC_SUBST([R8_FCFLAGS])
AC_SUBST([R4_FCFLAGS])

# Output files
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile
                 postprocessing/combine_blobs/Makefile
                 postprocessing/combine_restarts/Makefile
                 postprocessing/iceberg_comb/Makefile
                 postprocessing/land_utils/Makefile
                 postprocessing/list_ncvars/Makefile
                 postprocessing/mppnccombine/Makefile
                 postprocessing/plevel/Makefile
                 postprocessing/split_ncvars/Makefile
                 postprocessing/timavg/Makefile
               ])
AC_OUTPUT
