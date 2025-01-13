#***********************************************************************
#                   GNU Lesser General Public License
#
# This file is part of the GFDL FRE NetCDF tools package (FRE-NCTools).
#
# FRE-NCTools is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.
#
# FRE-NCTools is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
# for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with FRE-NCTools (LICENSE.md).  If not, see
# <http://www.gnu.org/licenses/>.
#***********************************************************************

# Prelude
AC_PREREQ([2.63])
AC_INIT(
  [FRE NCTools],
  [2024.05],
  [https://github.com/NOAA-GFDL/FRE-NCtools/issues],
  [fre-nctools],
  [https://github.com/NOAA-GFDL/FRE-NCtools])

# Place the git description information in the configure script
m4_define([git_revision],
  [m4_esyscmd_s([git describe --tags --always --dirty 2> /dev/null || echo "unknown version"])])
AC_REVISION([git_revision])
m4_define([git_hashval],
  [m4_esyscmd_s([git rev-parse HEAD 2> /dev/null || echo "unknown git hash"])])
AC_REVISION([git_hashval])
m4_define([git_date],
  [m4_esyscmd_s([git show --no-patch --format=%as HEAD 2> /dev/null || date +%Y-%m-%d])])
AC_REVISION([git_date])
m4_define([copyright_year],
  [m4_substr(git_date, 0, 4)])

AC_CONFIG_SRCDIR([src/land_utils/combine-ncc.F90])
AC_CONFIG_MACRO_DIR([m4])
AC_REQUIRE_AUX_FILE([tap-driver.sh])

# Remove -Werror from automake options until OS vendors brings in a newer
# version of automake.
AM_INIT_AUTOMAKE([foreign -Wall -Woverride tar-pax subdir-objects])
AM_EXTRA_RECURSIVE_TARGETS([docs install-docs])

# Check for programs
AC_PROG_RANLIB
AM_PROG_AR

AC_ARG_ENABLE([quad-precision],
  [AS_HELP_STRING([--enable-quad-precision],
    [build using higher (quad) precision])],
    [],
    [enable_quad_precision=no])

AC_PROG_CC([icc gcc nvc])
AM_PROG_CC_C_O
## When autoconf v2.70 is more available, this can be replaced with:
##AC_C__GENERIC
GX_C__GENERIC
AS_IF([test $gx_cv_c__Generic = no],
  AC_MSG_ERROR([The C compiler does not support the generic selection C-11 standard.  Please use a C-11 compliant compiler.])
)

AC_PROG_FC([ifort gfortran nvfortran])
AC_PROG_FC_C_O
AC_FC_SRCEXT([f90])
AC_FC_SRCEXT([F90])
AC_FC_LINE_LENGTH([unlimited])
AC_FC_FREEFORM()
AC_FC_MODULE_EXTENSION()
AC_FC_MODULE_FLAG()
AC_FC_MODULE_OUTPUT_FLAG()
GX_FC_DEFAULT_REAL_KIND8_FLAG()
GX_FC_DEFAULT_REAL_KIND4_FLAG()

AC_SEARCH_LIBS([ceil],[m])

# Check for the MPI compilers
# Default MPI build is off
build_mpi=no
LX_FIND_MPI()
AS_IF([test x"$with_mpi" = x"yes" && test x"$have_C_mpi" = x"no"], [
  AC_MSG_ERROR([MPI build requested but no MPI compiler found.])
])
AM_CONDITIONAL([WITH_MPI], [test x"$with_mpi" = x"yes"])
AM_CONDITIONAL([WITH_MPI_TESTS], [test x"$with_mpi" = x"yes" && test -z "$SKIP_MPI_TEST"])

# Build with OpenACC support.  Default is 'no'
AC_ARG_ENABLE([gpu],
  [AS_HELP_STRING([--enable-gpu],
    [Builds with OpenACC. This will result in a second executable for fregrid, fregrid_gpu. (default no)])])
AS_IF([test ${enable_gpu:-no} = yes],
  [enable_gpu=yes],
  [enable_gpu=no])
# check compile flags
AS_IF([test ${enable_gpu} = yes],
      [GX_OPENACC_FLAGS()])
AM_CONDITIONAL([ENABLE_GPU], [test "$enable_gpu" = "yes"])

AC_CHECK_FUNCS(gettid, [], [])

# Check if Fortran compiler has Cray pointer support
GX_FC_CRAY_POINTER_FLAG()

 # Add Cray Pointer support flag
if test ! -z "$FC_CRAY_POINTER_FLAG"; then
   FCFLAGS="$FCFLAGS $FC_CRAY_POINTER_FLAG"
fi

# Check for programs
AC_PROG_GREP
AC_PROG_EGREP
AC_PROG_SED

AC_CHECK_PROGS(ASCIIDOC, [asciidoc])
if test -n "$ASCIIDOC"; then
  AC_MSG_CHECKING([for asciidoc version])
  asciidoc_version=`$ASCIIDOC --version 2>/dev/null`
  case "${asciidoc_version}" in
  asciidoc' '8*)
    AC_MSG_RESULT([${asciidoc_version}])
    ;;
  *)
    AC_MSG_RESULT([${asciidoc_version} (unknown)])
    ;;
  esac
fi
AM_CONDITIONAL([WITH_ASCIIDOC], [test -n "$ASCIIDOC"])

# Although these are supplied with the NetCDF C library, we'll double
# check that ncgen and ncdump are available
AC_ARG_VAR([NCGEN], [Path to ncgen])
AC_PATH_PROG(NCGEN, [ncgen])
AS_IF([test -z "$NCGEN"],
  [AC_MSG_ERROR([Cannot find ncgen.  Please ensure ncgen is in PATH or set the environment variable NCGEN to the full path.])])

AC_ARG_VAR([NCDUMP], [Path to ncdump])
AC_PATH_PROG(NCDUMP, [ncdump])
AS_IF([test -z "$NCDUMP"],
  [AC_MSG_ERROR([Cannot find ncdump.  Please ensure ncdump is in PATH or set the environment variable NCDUMP to the full path.])])

AC_ARG_VAR([NCRCAT], [Path to ncrcat])
AC_PATH_PROG(NCRCAT, [ncrcat])
AS_IF([test -z "$NCRCAT"],
  [AC_MSG_ERROR([Cannot find ncrcat.  Please ensure ncrcat is in PATH or set the environment variable NCRCAT to the full path.])])

AC_ARG_VAR([NCATTED], [Path to ncatted])
AC_PATH_PROG(NCATTED, [ncatted])
AS_IF([test -z "$NCATTED"],
  [AC_MSG_ERROR([Cannot find ncatted.  Please ensure ncatted is in PATH or set the environment variable NCATTED to the full path.])])

AC_ARG_VAR([NCKS], [Path to ncks])
AC_PATH_PROG(NCKS, [ncks])
AS_IF([test -z "$NCKS"],
  [AC_MSG_ERROR([Cannot find ncks.  Please ensure ncks is in PATH or set the environment variable NCKS to the full path.])])

AC_ARG_VAR([PERL], [Path to perl])
AC_PATH_PROG(PERL, [perl])
AX_PROG_PERL_MODULES([List::MoreUtils], dnl
  eval AS_TR_CPP(HAVE_PLMOD_LIST_MOREUTILS)=yes,
  AC_MSG_WARN([Perl module List::MoreUtils not found.  Some applications will not be installed.]))
AM_CONDITIONAL([HAVE_PLMOD_LIST_MOREUTILS], [test x$HAVE_PLMOD_LIST_MOREUTILS = xyes])

AX_PYTHON_MODULE([xarray],[],[python3])
if test ${HAVE_PYMOD_XARRAY} = yes; then
  AM_CONDITIONAL([SKIP_XARRAY_TESTS], false )
else
  AM_CONDITIONAL([SKIP_XARRAY_TESTS], true )
fi

# Look for mpirun, mpiexec, srun application
AC_ARG_VAR([MPIRUN], [Path to mpirun])
AS_IF([test x"$with_mpi" = x"yes"], [
  AC_CHECK_PROGS(MPIRUN, [srun aprun mpirun mpiexec])
  AS_IF([test -z "$MPIRUN"],
    [AC_MSG_WARN([Cannot find a suitable MPI launch application.  MPI tests will be skipped.])])
])

# Check if openacc.h exists
if test "$enable_gpu" = yes ; then
   AC_CHECK_HEADERS([openacc.h], [], [AC_MSG_ERROR(["Cannot find OpenACC header file"])] )
   AC_MSG_WARN(Currently only NVIDIA compilers are supported to compile with OpenACC in FRE-NCTOOLS)
fi

# Check for getrusage.
AC_CHECK_FUNCS([getrusage])

# Check for library
AX_LIB_NETCDF()
if test "$with_netcdf" = "no"; then
  AC_MSG_ERROR([Unable to find NetCDF c library.])
fi

AC_LANG_PUSH([Fortran])
AX_LIB_NETCDFF()
if test "$with_netcdf_fortran" = "no"; then
  AC_MSG_ERROR([Unable to find NetCDF Fortran library.])
fi
AC_LANG_POP([Fortran])

# Check if Linux sched_getaffinity is available
AC_CHECK_FUNCS([sched_getaffinity], [], [])

# Check for typedefs, structures, and compiler characteristics
# Check if the C compiler supports a working `long double` with more range
# or precision than `double`.
AS_IF([test ${enable_quad_precision:-no} = yes], [
AC_TYPE_LONG_DOUBLE_WIDER
AS_IF([test $ac_cv_type_long_double_wider = no],
  [AC_MSG_ERROR([system does not support quad precision.  Try setting --enable-mixed-precision=no])
  ])
])

# Check for library functions

# Check if Test25-input link points to an existing directory
AM_CONDITIONAL([HAVE_TEST25_INPUT], [test -e ${srcdir}/t/Test25-input/grid_spec.nc])

# Set macros
AS_VAR_SET([GIT_REVISION], ["git_revision"])dnl
AC_DEFINE([GIT_REVISION], ["git_revision"],
          [Holds the 'git describe' information if configure ran within a git working directory])dnl
AC_SUBST(GIT_REVISION)dnl
AS_VAR_SET([GIT_HEADHASH], ["git_hashval"])dnl
AC_DEFINE([GIT_HEADHASH], ["git_hashval"],
          [Holds the 'git rev-parse HEAD' information if configure ran within a git working directory])dnl
AC_SUBST(GIT_HEADHASH)dnl
AS_VAR_SET([COPYRIGHT_YEAR], ["copyright_year"])dnl
AC_DEFINE([COPYRIGHT_YEAR], ["copyright_year"],
          [Year to place in Copyright notices printed with --version information])dnl
AC_SUBST(COPYRIGHT_YEAR)dnl

# Code for setting rpath based ncview's configure.in code.
AS_IF([test x$ac_compiler_gnu = xyes], [RPATH_FLAGS=""
  for word in $NETCDF_LDFLAGS $NETCDF_FCLDFLAGS
  do
    if test `expr $word : -L/` -eq 3
    then
      RPDIR=${word:2}
      regex=".*$RPDIR.*"
      if test `expr "$RPATH_FLAGS" : $regex` -eq 0
      then
        # RPDIR not in RPATH_FLAGS
        RPATH_FLAGS="$RPATHFLAGS -Wl,-rpath,$RPDIR"
      fi
    fi
  done])
AC_SUBST(RPATH_FLAGS) dnl

# Output files
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile
                 lib/libnfu/Makefile
                 man/Makefile
                 src/Makefile
                 tests/Makefile
                 tools/libfrencutils/Makefile
                 tools/libfrencutils_gpu/Makefile
                 tools/check_mask/Makefile
                 tools/cubic_utils/Makefile
                 tools/fregrid/Makefile
                 tools/fregrid_gpu/Makefile
                 tools/make_coupler_mosaic/Makefile
                 tools/make_hgrid/Makefile
                 tools/make_land_domain/Makefile
                 tools/make_quick_mosaic/Makefile
                 tools/make_regional_mosaic/Makefile
                 tools/make_solo_mosaic/Makefile
                 tools/make_topog/Makefile
                 tools/make_vgrid/Makefile
                 tools/mppncscatter/Makefile
                 tools/ncexists/Makefile
                 tools/nc_null_check/Makefile
                 tools/remap_land/Makefile
                 tools/river_regrid/Makefile
                 tools/runoff_regrid/Makefile
                 tools/transfer_to_mosaic_grid/Makefile
                 tools/simple_hydrog/Makefile
                 tools/simple_hydrog/postp/Makefile
                 tools/simple_hydrog/lakes/Makefile
                 tools/simple_hydrog/rmvpr/Makefile
                 tools/simple_hydrog/libfmslite/Makefile
                 t_gpu/Makefile
                 t_gpu/test_read_remap_file/Makefile
                 t_gpu/test_get_grid_cell_struct/Makefile
                 t_gpu/test_get_upbound_nxcells_2dx2d/Makefile
                 t_gpu/test_get_interp_order1/Makefile
])

AC_OUTPUT
