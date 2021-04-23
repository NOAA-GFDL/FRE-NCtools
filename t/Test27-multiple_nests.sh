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

# Test grid for multiple same level and telescoping nests 
load test_utils

@test "Test grid for multiple same level and telescoping nests" {

#Make_hgrid: create three same level -level1- nests in tiles 2,5,6"
  run command make_hgrid \
		--grid_type gnomonic_ed \
		--nlon 96 \
		--grid_name C48_grid \
		--do_schmidt \
		--stretch_factor 1.0 \
		--target_lon -97.5 \
		--target_lat 36.5 \
		--nest_grids 3 \
		--parent_tile 2,5,6 \
        --refine_ratio 2,2,2 \
        --istart_nest 7,13,7 \
        --jstart_nest 7,7,23 \
        --iend_nest 58,68,40 \
        --jend_nest 58,68,48 \
        --halo 3 \
        --great_circle_algorithm \
        --verbose 1
  [ "$status" -eq 0 ]

#Make_hgrid: create two same level -level1- nests in tiles 2,5 
#            and one -level2- telescoping nest in tile7"
#          ( tile7 refers to the first nest on the first level)
  run command make_hgrid \
		--grid_type gnomonic_ed \
		--nlon 96 \
		--grid_name C48_grid \
		--do_schmidt \
		--stretch_factor 1.0 \
		--target_lon -97.5 \
		--target_lat 36.5 \
		--nest_grids 3 \
		--parent_tile 2,5,7 \
        --refine_ratio 2,2,2 \
        --istart_nest 7,13,7 \
        --jstart_nest 7,7,23 \
        --iend_nest 58,68,40 \
        --jend_nest 58,68,48 \
        --halo 3 \
        --great_circle_algorithm \
        --verbose 1
  [ "$status" -eq 0 ]

}
