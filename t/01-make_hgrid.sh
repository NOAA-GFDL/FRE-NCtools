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
# License along with FRE-NCTools (LICENSE.md).  If not, see
# <http://www.gnu.org/licenses/>.
#***********************************************************************

teardown () {
  echo "$output"
  rm -f ocean_hgrid.nc
}

@test "make_hgrid exists and is executable" {
  run command -v make_hgrid
  [ "$status" -eq 0 ]
  run make_hgrid -h
  [ "$status" -eq 2 ]
}

@test "make_hgrid creates ocean_hgrid" {
  run command make_hgrid \
      --grid_type tripolar_grid \
      --nxbnd 2 \
      --nybnd 7 \
      --xbnd -280,80 \
      --ybnd -82,-30,-10,0,10,30,90 \
      --dlon 1.0,1.0 \
      --dlat 1.0,1.0,0.6666667,0.3333333,0.6666667,1.0,1.0 \
      --grid_name ocean_hgrid \
      --center c_cell
  [ "$status" -eq 0 ]
}
