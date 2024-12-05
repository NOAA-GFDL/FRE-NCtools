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
# same as Test24 except checks result with float tolerance, added in for compatibility with gcc

load test_utils

@test "Test fregrid ocean data" {

  cp $top_srcdir/t/Test20-input/*.nc .

#Create regular lat-lon grid (100:160, -15:15, size is 360x180)
   make_hgrid  \
                --grid_type regular_lonlat_grid  \
                --nxbnd 2  \
                --nybnd 2  \
                --xbnd 0,360  \
                --ybnd -2,2  \
                --nlon 720  \
                --nlat 10  \
                --grid_name latlon_grid

#Create lat-lon mosaic
   make_solo_mosaic  \
                --num_tiles 1  \
                --dir ./  \
                --mosaic_name latlon_mosaic  \
                --tile_file latlon_grid.nc

#Remap data from CM2.1 ocean grid onto regular lat-lon grid.
   fregrid   \
		--input_mosaic CM2.1_mosaic.nc   \
		--input_file ocean_temp_salt.res.nc   \
		--scalar_field temp,salt  \
		--output_file ocean_temp_salt.res.latlon.nc   \
		--output_mosaic latlon_mosaic.nc   \
		--check_conserve

   run_and_check [ -e ocean_temp_salt.res.latlon.nc ]
   run_and_check [ -e $top_srcdir/t/Test20-reference/ocean_temp_salt.res.latlon.nc ]

   # checks values with float tolerance to bypass any result differences from compilers
   run_and_check $top_srcdir/tests/nc_cmp -d -t 0.000001  ocean_temp_salt.res.latlon.nc  $top_srcdir/t/Test20-reference/ocean_temp_salt.res.latlon.nc

   run_and_check [ -e latlon_mosaic.nc ]
   run_and_check [ -e $top_srcdir/t/Test20-reference/latlon_mosaic.nc ]

   run_and_check $top_srcdir/tests/nc_cmp -d  latlon_mosaic.nc  $top_srcdir/t/Test20-reference/latlon_mosaic.nc
}
