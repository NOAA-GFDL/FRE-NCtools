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

teardown () {
  echo "$output"
  rm -f ocean_vgrid.nc 
}

@test "make_vgrid exists and is executable" {
  run command -v make_vgrid
  [ "$status" -eq 0 ]
  run make_vgrid -h
  [ "$status" -eq 1 ]
}

@test "make_vgrid creates ocean_vgrid" {
  run command make_vgrid \
      --nbnds 3 \
      --bnds 0.,220.,5500. \
      --dbnds 10.,10.,367.14286 \
      --center c_cell \
      --grid_name ocean_vgrid
  [ "$status" -eq 0 ]
}
