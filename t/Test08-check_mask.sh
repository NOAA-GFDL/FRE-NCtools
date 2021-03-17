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

# Test check_mask: create mask_table to mask out some all-land domain
# to save processor usage for a sea-ice model, baltic1 experiment

@test "Test check_ mask for baltic1 experiment" {

  if [ ! -d "Test08" ] 
  then
  		mkdir Test08
  fi

  cd Test08
  ncgen -o baltic1_grid_spec.nc $top_srcdir/t/Test08-input/baltic1_grid_spec.ncl

  run command check_mask \
		--grid_file baltic1_grid_spec.nc \
		--min_pe 60 \
		--max_pe 80 
  [ "$status" -eq 0 ]

  cd ..
  rm -rf Test08

}
