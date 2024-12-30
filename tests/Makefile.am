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
# License along with FRE-NCTools.  If not, see
# <http://www.gnu.org/licenses/>.
#***********************************************************************

if SKIP_XARRAY_TESTS
skipflag="skip"
else
skipflag=""
endif

if HAVE_PLMOD_LIST_MOREUTILS
have_plmod_list_moreutils="yes"
else
have_plmod_list_moreutils="no"
endif

AM_CFLAGS = $(NETCDF_CFLAGS) -I$(top_builddir)/src -I$(top_srcdir)/tools/libfrencutils
AM_FCFLAGS = -I$(top_builddir)/src
LDADD = -L$(top_builddir)/src -lver $(NETCDF_LDFLAGS) $(NETCDF_LIBS)

check_PROGRAMS = \
  create_5degree_topog \
  create_daily_tile_files \
  create_hydro_glcc_file \
  create_land_remap_files \
  create_mask_topog \
  create_ocean_mask_grid \
  create_plevel_test_ncfile \
  create_regional_file \
  create_river_regrid_file \
  create_runoff_file \
  create_temp_salt_file \
  create_timeavg_test_ncfiles \
  c_test_version \
  f_test_version

create_plevel_test_ncfile_SOURCES=create_plevel_test_ncfile.c
create_timeavg_test_ncfiles_SOURCES=create_timeavg_test_ncfiles.c
c_test_version_SOURCES=libtests/c_test_version.c
f_test_version_SOURCES=libtests/f_test_version.f90

TESTS_ENVIRONMENT = \
    top_srcdir=$(abs_top_srcdir); export top_srcdir; \
    builddir=$(abs_builddir); export builddir; \
    MPIRUN=$(MPIRUN); export MPIRUN; \
    skipflag=$(skipflag); export skipflag; \
    have_plmod_list_moreutils=$(have_plmod_list_moreutils); export have_plmod_list_moreutils; \
    PKGLIBEXECDIR="$(abs_top_builddir)/src"; export PKGLIBEXECDIR; \
    for d in $$( find $(abs_top_builddir)/tools -mindepth 1 -not -path '*/\.*' -type d ); do \
      PATH="$$d"'$(PATH_SEPARATOR)'"$$PATH"; \
    done; \
    PATH="$(abs_top_builddir)/src"'$(PATH_SEPARATOR)'"$$PATH"; export PATH;

tools_tests = \
    check_mask/ocean-mosaic-mask \
    check_mask/ocean-grid-mask \
    make_hgrid/c48-hgrid \
    make_hgrid/c48-nest \
    make_hgrid/c48-nest-telescoping \
    make_hgrid/from-file \
    make_hgrid/ocean-hgrid \
    make_vgrid/ocean-vgrid \
    make_solo_mosaic/ocean-mosaic \
    make_solo_mosaic/c48-mosaic \
    make_topog/ocean-topog \
    make_topog/ocean-topog-mpi \
    make_coupler_mosaic/coupled-model \
    make_coupler_mosaic/coupled-model-mpi \
    make_coupler_mosaic/coupled-model-nest \
    make_regional_mosaic/make_regional_mosaic \
    make_quick_mosaic/coupled-model-quick \
    mppncscatter/xyt-decomp \
    fregrid/cubedsphere \
    fregrid/extrapolate \
    fregrid/latlon \
    fregrid/remap-file \
    fregrid/remap-file-mpi \
    fregrid/remap-file-stretchgrid \
    remap_land/remap_land \
    river_regrid/river_regrid \
    runoff_regrid/runoff_regrid

combine_tests = \
    combine-ncc/combine \
    combine_restarts/combine \
    combine_restarts/dryrun \
    combine_restarts/extra_options \
    iceberg_comb/check-valid-file \
    iceberg_comb/combine \
    iceberg_comb/combine-no-icebergs \
    iceberg_comb/missing-file \
    iceberg_comb/outfile-exists \
    iceberg_comb/setnum-mismatch \
    is-compressed/compressed \
    is-compressed/not-compressed \
    mppnccombine/dry-run \
    mppnccombine/list-output-files \
    mppnccombine/input-output-kinds \
    mppnccombine/missing-file-fail \
    mppnccombine/missing-file-force \
    mppnccombine/non-zero-start-number \
    mppnccombine/output-kinds \
    mppnccombine/print-memory-usage \
    mppnccombine/remove-decomp-files \
    mppnccombine/six-digit-extension \
    mppnccombine/with-end-number \
    mppnccombine/xt-decomp \
    mppnccombine/xyt-decomp \
    mppnccombine/yt-decomp

land_util_tests = \
    decompress-ncc/decompress \
    scatter-ncc/scatter_1x1 \
    scatter-ncc/scatter_1x2 \
    scatter-ncc/scatter_2x1

TESTS = \
    libtests/print_version \
    list_ncvars/list_ncvars_ftn \
    list_ncvars/list_ncvars_sh \
    plevel/plevel \
    plevel/plevel-0 \
    plevel/plevel-3 \
    plevel/plevel-a \
    plevel/plevel-o \
    plevel/plevel-nctypes \
    plevel/plevel-t \
    split_ncvars/split_ncvars \
    split_ncvars/split_ncvars-i \
    split_ncvars/split_ncvars-l \
    split_ncvars/split_ncvars-o \
    split_ncvars/split_ncvars-p \
    split_ncvars/split_ncvars-s \
    split_ncvars/split_ncvars-v \
    timeavg/timavg \
    timeavg/timavg-m \
    $(combine_tests) \
    $(land_util_tests) \
    $(tools_tests)

EXTRA_DIST = \
  init.sh \
  init.cfg \
  nc_cmp \
  make_hgrid/create_fv3_grid_spec_files \
  make_hgrid/grid_spec.tile1.ncl \
  make_hgrid/grid_spec.tile2.ncl \
  make_hgrid/grid_spec.tile3.ncl \
  make_hgrid/grid_spec.tile4.ncl \
  make_hgrid/grid_spec.tile5.ncl \
  make_hgrid/grid_spec.tile6.ncl \
  $(TESTS)

# Shortcut targets to make it easier to run (very) expensive tests.
check-expensive:
	$(MAKE) check RUN_EXPENSIVE_TESTS=yes
check-very-expensive:
	$(MAKE) check-expensive RUN_VERY_EXPENSIVE_TESTS=yes
