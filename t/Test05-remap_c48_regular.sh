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

export SETUP_FNCT="generate_all_from_ncl 03"
load test_utils

@test "remap data from C48 to regular lat-lon grid" {

  mkdir output

  run command fregrid \
		--input_mosaic C48_mosaic.nc \
		--input_file 19800101.atmos_daily \
		--scalar_field zsurf,temp,t_surf \
		--nlon 144 \
		--nlat 90 \
		--interp_method conserve_order2 \
		--output_dir output \
		--output_file 19800101.atmos_daily.nc \
		--check_conserve \
		--remap_file C48_to_N45_remap.nc 
  [ "$status" -eq 0 ]

}
