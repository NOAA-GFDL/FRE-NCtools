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

# Test make_regional_mosaic: create mosaic file for regional output and do remapping
#             for regional output

@test "Test make_regional_mosaic" {
  skip "TO DO: Files are too large" 

  if [ ! -d "Test11" ] 
  then
  		mkdir Test11
  fi

  cd Test11
  cp $top_srcdir/t/Test11-input/* . 

  ncgen -o rregionatmos_4xdaily_eq_avg.tile2.nc rregionatmos_4xdaily_eq_avg.tile2.cdl
  ncgen -o rregionatmos_4xdaily_eq_avg.tile3.nc rregionatmos_4xdaily_eq_avg.tile3.cdl
  ncgen -o rregionatmos_4xdaily_eq_avg.tile5.nc rregionatmos_4xdaily_eq_avg.tile5.cdl
  ncgen -o rregionatmos_4xdaily_eq_avg.tile6.nc rregionatmos_4xdaily_eq_avg.tile6.cdl

  run command make_regional_mosaic \
		--global_mosaic C192-rot-120.-0.-r1_mosaic.nc \
		--regional_file rregionatmos_4xdaily_eq_avg.tile2.nc
  [ "$status" -eq 0 ] 

  run command make_regional_mosaic \
		--global_mosaic C192-rot-120.-0.-r1_mosaic.nc \
		--regional_file rregionatmos_4xdaily_eq_avg.tile3.nc
  [ "$status" -eq 0 ] 

  run command make_regional_mosaic  \
		--global_mosaic C192-rot-120.-0.-r1_mosaic.nc  \
		--regional_file rregionatmos_4xdaily_eq_avg.tile5.nc
  [ "$status" -eq 0 ] 

  run command make_regional_mosaic  \
		--global_mosaic C192-rot-120.-0.-r1_mosaic.nc  \
		--regional_file rregionatmos_4xdaily_eq_avg.tile6.nc
  [ "$status" -eq 0 ] 

#Create regional mosaic
  run command make_solo_mosaic  \
		--num_tiles 4  \
		--dir ./  \
		--mosaic regional_mosaic  \
		--tile_file regional_grid.tile2.nc,regional_grid.tile3.nc,regional_grid.tile5.nc,regional_grid.tile6.nc
  [ "$status" -eq 0 ] 

#Create regular lat-lon grid (100:160, -15:15, size is 360x180)
  run command make_hgrid  \
		--grid_type regular_lonlat_grid  \
		--nxbnd 2  \
		--nybnd 2  \
		--xbnd 0,360  \
		--ybnd -2,2  \
		--nlon 720  \
		--nlat 10  \
		--grid_name latlon_grid
  [ "$status" -eq 0 ] 

#Create lat-lon mosaic
  run command make_solo_mosaic  \
		--num_tiles 1  \
		--dir ./  \
		--mosaic_name latlon_mosaic  \
		--tile_file latlon_grid.nc
  [ "$status" -eq 0 ] 

#Remap data
  run command fregrid  \
		--input_mosaic regional_mosaic.nc  \
		--output_mosaic latlon_mosaic.nc  \
		--input_file rregionatmos_4xdaily_eq_avg  \
		--output_file atmos_4xdaily_latlon.nc  \
		--scalar_field ps  \
		--check_conserve
  [ "$status" -eq 0 ] 

  cd ..
  rm -rf Test11

}


