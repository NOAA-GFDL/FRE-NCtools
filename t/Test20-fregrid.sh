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

# test regrid ocean restart file 
load input_util
SETUP_FNCT="cp $top_srcdir/t/Test20-input/*.nc ."

@test "Test fregrid ocean data" {

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

#Remap data from CM2.1 ocean grid onto regular lat-lon grid.
  run command fregrid   \
		--input_mosaic CM2.1_mosaic.nc   \
		--input_file ocean_temp_salt.res.nc   \
		--scalar_field temp,salt  \
		--output_file ocean_temp_salt.res.latlon.nc   \
		--output_mosaic latlon_mosaic.nc   \
		--check_conserve
  [ "$status" -eq 0 ]

}
