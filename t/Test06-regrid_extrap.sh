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

# Test remap data onto cm2m ocean grid with extrapolation and vertical interpolation

@test "Test remap data onto cm2m ocean grid with extrapolation and vertical interpolation" {

  if [ ! -d "Test06" ] 
  then
  		mkdir Test06
  fi

  cd Test06
  ncgen -o ocean_hgrid.nc $top_srcdir/t/Test06-input/ocean_hgrid.ncl
  ncgen -o ocean_mosaic.nc $top_srcdir/t/Test06-input/ocean_mosaic.ncl
  ncgen -o ocean_vgrid.nc $top_srcdir/t/Test06-input/ocean_vgrid.ncl
  ncgen -o WOA09_ann_theta.nc $top_srcdir/t/Test06-input/WOA09_ann_theta.ncl

  run command make_hgrid \
		--grid_type regular_lonlat_grid \
		--nxbnd 2 \
		--nybnd 2 \
		--xbnd 0,360 \
		--ybnd -90,90 \
		--nlon 720 \
		--nlat 360 \
		--grid_name levitus_grid
  [ "$status" -eq 0 ]

  run command make_solo_mosaic \
		--num_tiles 1 \
		--dir . \
		--mosaic_name levitus_mosaic \
		--tile_file levitus_grid.nc \
		--periodx 360
  [ "$status" -eq 0 ]

  run command fregrid \
		--input_mosaic levitus_mosaic.nc \
		--input_file WOA09_ann_theta.nc \
		--scalar_field POTENTIAL_TEMP \
		--output_file WOA09_ann_theta_cm2g_extrap.nc \
		--output_mosaic ocean_mosaic.nc \
		--extrapolate \
		--dst_vgrid ocean_vgrid.nc \
		--check_conserve
  [ "$status" -eq 0 ]

  cd ..
  rm -rf Test06

}
