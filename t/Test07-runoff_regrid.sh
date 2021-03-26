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

@test "Test  remap runoff data from regular lat-lon grid onto cm2m grid" {

  if [ ! -d "Test07" ] 
  then
  		mkdir Test07
  fi

  cd Test07
  ncgen -o ocean_hgrid.nc $top_srcdir/t/Test06-input/ocean_hgrid.ncl
  ncgen -o ocean_mosaic.nc $top_srcdir/t/Test06-input/ocean_mosaic.ncl
  ncgen -o ocean_vgrid.nc $top_srcdir/t/Test06-input/ocean_vgrid.ncl
  ncgen -o C48_grid.tile1.nc $top_srcdir/t/Test03-input/C48_grid.tile1.ncl
  ncgen -o C48_grid.tile2.nc $top_srcdir/t/Test03-input/C48_grid.tile2.ncl
  ncgen -o C48_grid.tile3.nc $top_srcdir/t/Test03-input/C48_grid.tile3.ncl
  ncgen -o C48_grid.tile4.nc $top_srcdir/t/Test03-input/C48_grid.tile4.ncl
  ncgen -o C48_grid.tile5.nc $top_srcdir/t/Test03-input/C48_grid.tile5.ncl
  ncgen -o C48_grid.tile6.nc $top_srcdir/t/Test03-input/C48_grid.tile6.ncl
  ncgen -o C48_mosaic.nc $top_srcdir/t/Test03-input/C48_mosaic.ncl
  ncgen -o runoff.daitren.iaf.nc $top_srcdir/t/Test07-input/runoff.daitren.iaf.ncl
  ncgen -o topog.nc $top_srcdir/t/Test07-input/topog.ncl

  run command runoff_regrid \
		--input_file runoff.daitren.iaf.nc \
		--input_fld_name runoff \
		--output_mosaic ocean_mosaic.nc \
		--output_topog topog.nc \
		--output_file runoff.cm2m.nc
  [ "$status" -eq 0 ]

  cd ..
  rm -rf Test07
}
