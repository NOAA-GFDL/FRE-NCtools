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

# test stretched grid

@test "Test stretched grid data lats 32.0 34.0 35.4" {

  if [ ! -d "Test31" ]
  then
  		mkdir Test31
  fi

  cd Test31
   cp $top_srcdir/t/Test31-input/ocean_hgrid.nc .
   cp $top_srcdir/t/Test31-input/ocean_mosaic.nc .
   cp $top_srcdir/t/Test31-input/ocean_topog.nc .

#Make streetched grid
  run_and_check make_hgrid \
               --grid_type gnomonic_ed --do_schmidt \
               --stretch_factor 2.5 \
               --target_lon 262.4 \
               --target_lat 32.0 \
               --nlon 512 \
               --grid_name C256_grid_32.0

  run_and_check make_hgrid \
               --grid_type gnomonic_ed --do_schmidt \
               --stretch_factor 2.5 \
               --target_lon 262.4 \
               --target_lat 32.0 \
               --nlon 512 \
               --grid_name C256_grid_34.0


 run_and_check make_hgrid \
               --grid_type gnomonic_ed --do_schmidt \
               --stretch_factor 2.5 \
               --target_lon 262.4 \
               --target_lat 35.4 \
               --nlon 512 \
               --grid_name C256_grid_35.4

#Create stretched grid mosaic
  run_and_check make_solo_mosaic \
                --num_tiles 6 \
                --dir ./ \
                --mosaic_name C256_mosaic_32.0 --tile_file C256_grid_32.0.tile1.nc,C256_grid_32.0.tile2.nc,C256_grid_32.0.tile3.nc,C256_grid_32.0.tile4.nc,C256_grid_32.0.tile5.nc,C256_grid_32.0.tile6.nc


  run_and_check make_solo_mosaic \
                --num_tiles 6 \
                --dir ./ \
                --mosaic_name C256_mosaic_34.0 --tile_file C256_grid_34.0.tile1.nc,C256_grid_34.0.tile2.nc,C256_grid_34.0.tile3.nc,C256_grid_34.0.tile4.nc,C256_grid_34.0.tile5.nc,C256_grid_34.0.tile6.nc


  run_and_check make_solo_mosaic \
                --num_tiles 6 \
                --dir ./ \
                --mosaic_name C256_mosaic_35.4 --tile_file C256_grid_35.4.tile1.nc,C256_grid_35.4.tile2.nc,C256_grid_35.4.tile3.nc,C256_grid_35.4.tile4.nc,C256_grid_35.4.tile5.nc,C256_grid_35.4.tile6.nc

# stretched grid lats 32.0 34.0 35.4
  run bash -c 'fregrid \
                --input_mosaic C256_mosaic_32.0.nc \
                --nlon 640 \
                --nlat 400 \
                --latBegin 15.0 \
                --latEnd 65.0 \
                --lonBegin 230.0 \
                --lonEnd 310.0 \
                --remap_file fregrid_remap_file_640_by_400_32.0.nc \
                --output_file out_32.0.nc --check_conserve \
                | awk 'NR==4' | rev | cut -c39-49 | rev'
  var_32_0=$(echo $output | awk '{ print sprintf("%.9f", $1); }')
  echo $output
  run_and_check expr ${var_32_0} \< 0.00001

  run bash -c 'fregrid \
                --input_mosaic C256_mosaic_34.0.nc \
                --nlon 640 \
                --nlat 400 \
                --latBegin 15.0 \
                --latEnd 65.0 \
                --lonBegin 230.0 \
                --lonEnd 310.0 \
                --remap_file fregrid_remap_file_640_by_400_34.0.nc \
                --output_file out_34.0.nc --check_conserve \
                | awk 'NR==4' | rev | cut -c39-49 | rev'
  var_34_0=$(echo $output | awk '{ print sprintf("%.9f", $1); }')
  echo $output
  run_and_check expr ${var_34_0} \< 0.00001


  run bash -c 'fregrid \
                --input_mosaic C256_mosaic_35.4.nc \
                --nlon 640 \
                --nlat 400 \
                --latBegin 15.0 \
                --latEnd 65.0 \
                --lonBegin 230.0 \
                --lonEnd 310.0 \
                --remap_file fregrid_remap_file_640_by_400_35.4.nc \
                --output_file out_35.4.nc --check_conserve \
                | awk 'NR==4' | rev | cut -c39-49 | rev'
  var_35_4=$(echo $output | awk '{ print sprintf("%.9f", $1); }')
  echo $output
  run_and_check expr ${var_35_4} \< 0.00001

  cd ..
#  rm -rf Test31
}
