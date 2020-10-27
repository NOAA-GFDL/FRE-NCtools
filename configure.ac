# Prelude
AC_PREREQ([2.63])
AC_INIT([FRE NC Tools], [18.0.0], [oar.gfdl.help@noaa.gov], [fre-nctools], [https://github.com/NOAA-GFDL/FRE-NCtools])
# Place the git description informatin in the configure script
m4_define([git_revision],
  [m4_esyscmd_s([git describe --always --dirty 2> /dev/null || echo "unknown version"])])
AC_REVISION([$git_revision])

AC_CONFIG_SRCDIR([postprocessing/land_utils/combine-ncc.F90])
AC_CONFIG_MACRO_DIR([m4])
AC_REQUIRE_AUX_FILE([tap-driver.sh])

AM_INIT_AUTOMAKE([foreign -Wall -Werror -Woverride])
AM_EXTRA_RECURSIVE_TARGETS([docs install-docs])

# Check for programs
AC_PROG_RANLIB
AM_PROG_AR

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

AC_CHECK_FUNCS(gettid, [], [])


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

AX_LIB_NETCDFF()
if test "$with_netcdf_fortran" = "no"; then
  AC_MSG_ERROR([Unable to find NetCDF Fortran library.])
fi

# Check for typedefs, structures, and compiler characteristics

# Check for library functions

# Other defs
AC_DEFINE([GIT_REVISION], "git_revision",
          [Holds the 'git describe' information if configure ran within a git working directory])
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
 

])
AC_OUTPUT
