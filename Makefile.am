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
ACLOCAL_AMFLAGS = -I m4

SUBDIRS = lib \
          src \
          tools/make_topog \
          tools/make_vgrid \
          tools/mppncscatter \
          tools/ncexists \
          tools/remap_land \
          tools/river_regrid \
          tools/runoff_regrid \
          tools/transfer_to_mosaic_grid \
          tools/simple_hydrog \
          tools/simple_hydrog/libfmslite \
          tools/simple_hydrog/postp \
          tools/simple_hydrog/lakes \
          tools/simple_hydrog/rmvpr \
          tools/nc_null_check

SUBDIRS += man \
           tests

if ENABLE_GPU
SUBDIRS += t_gpu
endif

# Shortcut targets to make it easier to run (very) expensive tests.
check-expensive:
	$(MAKE) check RUN_EXPENSIVE_TESTS=yes
check-very-expensive:
	$(MAKE) check-expensive RUN_VERY_EXPENSIVE_TESTS=yes
