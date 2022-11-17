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

load test_utils

@test "Test make_hgrid do_transform with referece data" {

##if [ "$CC" != "icc" ]; then skip "Test fails reference check without icc"; fi

 #Create a equal distance gnomonic cubic grid using the do_cube_transform option
  make_hgrid \
      --grid_type gnomonic_ed \
      --nlon 256 \
      --do_cube_transform \
      --stretch_factor 3 \
      --target_lat 17.5 \
      --target_lon 172.5 \
      --grid_name C128_ctr_grid

   [ -e  C128_ctr_grid.tile1.nc ]
   [ -e $top_srcdir/t/Test33-reference/C128_ctr_grid.tile1.nc ]
   nccmp -V
   nccmp -d  C128_ctr_grid.tile1.nc $top_srcdir/t/Test33-reference/C128_ctr_grid.tile1.nc

   [ -e  C128_ctr_grid.tile2.nc ]
   [ -e $top_srcdir/t/Test33-reference/C128_ctr_grid.tile2.nc ]
   nccmp -V
   nccmp -d  C128_ctr_grid.tile2.nc $top_srcdir/t/Test33-reference/C128_ctr_grid.tile2.nc

   [ -e  C128_ctr_grid.tile3.nc ]
   [ -e $top_srcdir/t/Test33-reference/C128_ctr_grid.tile3.nc ]
   nccmp -V
   nccmp -d  C128_ctr_grid.tile3.nc $top_srcdir/t/Test33-reference/C128_ctr_grid.tile3.nc

   [ -e  C128_ctr_grid.tile4.nc ]
   [ -e $top_srcdir/t/Test33-reference/C128_ctr_grid.tile4.nc ]
   nccmp -V
   nccmp -d  C128_ctr_grid.tile4.nc $top_srcdir/t/Test33-reference/C128_ctr_grid.tile4.nc

   [ -e  C128_ctr_grid.tile5.nc ]
   [ -e $top_srcdir/t/Test33-reference/C128_ctr_grid.tile5.nc ]
   nccmp -V
   nccmp -d  C128_ctr_grid.tile5.nc $top_srcdir/t/Test33-reference/C128_ctr_grid.tile5.nc

   [ -e  C128_ctr_grid.tile6.nc ]
   [ -e $top_srcdir/t/Test33-reference/C128_ctr_grid.tile6.nc ]
   nccmp -V
   nccmp -d  C128_ctr_grid.tile6.nc $top_srcdir/t/Test33-reference/C128_ctr_grid.tile6.nc

}
