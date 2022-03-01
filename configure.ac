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
  [FRE NC Tools],
  [2.20.0],
  [oar.gfdl.help@noaa.gov],
  [fre-nctools],
  [https://github.com/NOAA-GFDL/FRE-NCtools])

# Place the git description informatin in the configure script
m4_define([git_revision],
  [m4_esyscmd_s([git describe --always --dirty 2> /dev/null || echo "unknown version"])])
AC_REVISION([$git_revision])

m4_define([git_hashval],
  [m4_esyscmd_s([git rev-parse HEAD 2> /dev/null || echo "unknown githash"])])
AC_REVISION([$git_hashval])

AC_CONFIG_SRCDIR([postprocessing/land_utils/combine-ncc.F90])
AC_CONFIG_MACRO_DIR([m4])
AC_REQUIRE_AUX_FILE([tap-driver.sh])

AM_INIT_AUTOMAKE([foreign -Wall -Werror -Woverride tar-ustar])
AM_EXTRA_RECURSIVE_TARGETS([docs install-docs])

# Check for programs
AC_PROG_RANLIB
AM_PROG_AR

AC_ARG_ENABLE([quad-precision],
  [AS_HELP_STRING([--enable-quad-precision],
    [build using higher (quad) precision])],
    [],
    [enable_quad_precision=no])

AC_ARG_ENABLE([ocean-model-grid-generator],
  [AS_HELP_STRING([--enable-ocean-model-grid-generator],
    [build and install the ocean model grid generator tool])])
AS_IF([test ${enable_ocean_model_grid_generator:-yes} = yes],
  [enable_ocean_model_grid_generator=yes],
  [enable_ocean_model_grid_generator=no])
AM_CONDITIONAL([ENABLE_OCEAN_MODEL_GRID_GENERATOR],
  [test $enable_ocean_model_grid_generator = yes])

AC_PROG_CC([icc gcc])
AM_PROG_CC_C_O

AC_PROG_FC([ifort gfortran])
AC_PROG_FC_C_O
AC_FC_SRCEXT([f90])
AC_FC_SRCEXT([F90])
AC_FC_LINE_LENGTH([unlimited])
AC_FC_FREEFORM()
AC_FC_MODULE_EXTENSION()
GX_FC_DEFAULT_REAL_KIND8_FLAG()
GX_FC_DEFAULT_REAL_KIND4_FLAG()

AC_SEARCH_LIBS([ceil],[m])

# Check for the MPI compilers
# Default MPI build is off
build_mpi=no
LX_FIND_MPI()
AS_IF([test x"$with_mpi" = x"yes"],
   [if test x"$have_C_mpi" = x"yes"; then
      build_mpi="yes"
    fi])
AM_CONDITIONAL([WITH_MPI], [test x"$build_mpi" = x"yes"])
AM_CONDITIONAL([WITH_MPI_TESTS], [test x"$build_mpi" = x"yes" -a -z "$SKIP_MPI_TEST"])

AC_CHECK_FUNCS(gettid, [], [])

# Check if Fortran compiler has Cray pointer support
GX_FC_CRAY_POINTER_FLAG()


 # Add Cray Pointer support flag
if test ! -z "$FC_CRAY_POINTER_FLAG"; then
   FCFLAGS="$FCFLAGS $FC_CRAY_POINTER_FLAG"
fi


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

AC_CHECK_PROGS(BATS, [bats])
AC_CHECK_PROGS(PROVE, [prove])

AM_CONDITIONAL([WITH_CHECK_PROGS], [test -n "$PROVE" -a -n "$BATS"])

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

# Other defs
AC_DEFINE([GIT_REVISION], "git_revision",
          [Holds the 'git describe' information if configure ran within a git working directory])
AC_DEFINE([GIT_HEADHASH], "git_hashval",
          [Holds the 'git rev-parse HEAD' information if configure ran within a git working directory])

AS_IF([test $enable_ocean_model_grid_generator = yes],
  AC_CONFIG_SUBDIRS([tools/ocean_model_grid_generator])
)

# Output files
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile
                 postprocessing/combine_restarts/Makefile
                 postprocessing/iceberg_comb/Makefile
                 postprocessing/land_utils/Makefile
                 postprocessing/list_ncvars/Makefile
                 postprocessing/mppnccombine/Makefile
                 postprocessing/plevel/Makefile
                 postprocessing/split_ncvars/Makefile
                 postprocessing/timavg/Makefile
                 tools/libfrencutils/Makefile
                 tools/check_mask/Makefile
                 tools/cubic_utils/Makefile
                 tools/fregrid/Makefile
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
                 t/Makefile
                 tools/simple_hydrog/Makefile
                 tools/simple_hydrog/postp/Makefile
                 tools/simple_hydrog/lakes/Makefile
                 tools/simple_hydrog/rmvpr/Makefile
                 tools/simple_hydrog/libfmslite/Makefile
])
AC_OUTPUT
