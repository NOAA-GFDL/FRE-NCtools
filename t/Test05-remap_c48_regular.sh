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

@test "remap data from C48 to regular lat-lon grid" {

  if [ ! -d "Test05" ] 
  then
  		mkdir Test05
  fi

  cd Test05
  ncgen -o OCCAM_p5degree.nc $top_srcdir/t/Test03-input/OCCAM_p5degree.ncl
  ncgen -o C48_mosaic.nc $top_srcdir/t/Test03-input/C48_mosaic.ncl
  ncgen -o C48_grid.tile1.nc $top_srcdir/t/Test03-input/C48_grid.tile1.ncl
  ncgen -o C48_grid.tile2.nc $top_srcdir/t/Test03-input/C48_grid.tile2.ncl
  ncgen -o C48_grid.tile3.nc $top_srcdir/t/Test03-input/C48_grid.tile3.ncl
  ncgen -o C48_grid.tile4.nc $top_srcdir/t/Test03-input/C48_grid.tile4.ncl
  ncgen -o C48_grid.tile5.nc $top_srcdir/t/Test03-input/C48_grid.tile5.ncl
  ncgen -o C48_grid.tile6.nc $top_srcdir/t/Test03-input/C48_grid.tile6.ncl
  ncgen -o 19800101.atmos_daily.tile1.nc $top_srcdir/t/Test03-input/19800101.atmos_daily.tile1.ncl
  ncgen -o 19800101.atmos_daily.tile2.nc $top_srcdir/t/Test03-input/19800101.atmos_daily.tile2.ncl
  ncgen -o 19800101.atmos_daily.tile3.nc $top_srcdir/t/Test03-input/19800101.atmos_daily.tile3.ncl
  ncgen -o 19800101.atmos_daily.tile4.nc $top_srcdir/t/Test03-input/19800101.atmos_daily.tile4.ncl
  ncgen -o 19800101.atmos_daily.tile5.nc $top_srcdir/t/Test03-input/19800101.atmos_daily.tile5.ncl
  ncgen -o 19800101.atmos_daily.tile6.nc $top_srcdir/t/Test03-input/19800101.atmos_daily.tile6.ncl

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

  cd ..
  rm -rf Test05
}

