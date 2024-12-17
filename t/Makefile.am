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
if WITH_MPI_TESTS
skip_MPI=""
else
skip_MPI="skip"
endif

if SKIP_XARRAY_TESTS
skipflag="skip"
else
skipflag=""
endif
TEST_EXTENSIONS = .sh

TESTS_ENVIRONMENT = top_srcdir=$(abs_top_srcdir); export top_srcdir; \
                    skip_mpi=$(skip_MPI); export skip_mpi; \
                    skipflag=$(skipflag); export skipflag; \
    for d in $$( find $(abs_top_builddir)/tools -mindepth 1 -not -path '*/\.*' -type d ); do \
      PATH="$$d"'$(PATH_SEPARATOR)'"$$PATH"; \
    done; \
                    PATH="$(abs_top_builddir)/src"'$(PATH_SEPARATOR)'"$$PATH"; export PATH;


if HAVE_TEST25_INPUT
# Data for test 25 is only available on some systems as the directory
# is a symlink
test25_data = Test25-input
endif

SH_LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
                $(top_srcdir)/tap-driver.sh

TESTS = Test01-check_programs_exist.sh \
        Test05-remap_c48_regular.sh \
        Test06-regrid_extrap.sh \
        Test07-runoff_regrid.sh \
        Test10-remap_land_res.sh \
        Test11-make_regional_mosaic.sh \
        Test12-mppnscatter.sh \
        Test13-make_quick_mosaic.sh \
        Test14-remap_ocean.sh \
        Test15-regrid_land.sh \
        Test16-river_regrid.sh \
        Test25-hydrology.sh \
        Test31-fregrid_stretched.sh \
        Test32-fregrid_no_stretched.sh

EXTRA_DIST = $(TESTS) \
             test_utils.bash \
             Test03-input \
             Test06-input \
             Test07-input \
             Test10-input \
             Test12-input \
             Test15-input \
             Test16-input \
             $(test25_data) \
             Test31-input \
             Test32-input \
             Test33-reference

clean-local:
	$(RM) *.nc Test*/*
