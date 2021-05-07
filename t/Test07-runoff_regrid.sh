#!/usr/bin/env bats

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

# Test remap runoff data from regular lat-lon grid onto cm2m grid
load test_utils
SETUP_FNCT="runoff_regrid_input"

@test "Test remap runoff data from regular lat-lon grid onto cm2m grid" {

  run command runoff_regrid \
		--input_file runoff.daitren.iaf.nc \
		--input_fld_name runoff \
		--output_mosaic ocean_mosaic.nc \
		--output_topog topog.nc \
		--output_file runoff.cm2m.nc
  [ "$status" -eq 0 ]

}

runoff_regrid_input() {
    generate_all_from_ncl 06
    generate_all_from_ncl 03
    generate_all_from_ncl 07
}
