#!/usr/bin/env bats
# test stretched grid 

@test "Test stretched grid data lats 32.0 34.0 35.4" {

  if [ ! -d "Test28" ] 
  then
  		mkdir Test28
  fi

  cd Test28
   cp $top_srcdir/t/Test28-input/ocean_hgrid.nc . 
   cp $top_srcdir/t/Test28-input/ocean_mosaic.nc .
   cp $top_srcdir/t/Test28-input/topog.nc .

#Make streetched grid 
  run command make_hgrid \
               --grid_type gnomonic_ed --do_schmidt \
               --stretch_factor 2.5 \
               --target_lon 262.4 \
               --target_lat 32.0 \
               --nlon 512 \
               --grid_name C256_grid_32.0
  [ "$status" -eq 0 ]

  run command make_hgrid \
               --grid_type gnomonic_ed --do_schmidt \
               --stretch_factor 2.5 \
               --target_lon 262.4 \
               --target_lat 32.0 \
               --nlon 512 \
               --grid_name C256_grid_34.0
  [ "$status" -eq 0 ]

  run command make_hgrid \
               --grid_type gnomonic_ed --do_schmidt \
               --stretch_factor 2.5 \
               --target_lon 262.4 \
               --target_lat 35.4 \
               --nlon 512 \
               --grid_name C256_grid_35.4
  [ "$status" -eq 0 ]


#Create stretched grid mosaic
  run command make_solo_mosaic \
                --num_tiles 6 \
                --dir ./ \
                --mosaic_name C256_mosaic_32.0 --tile_file C256_grid_32.0.tile1.nc,C256_grid_32.0.tile2.nc,C256_grid_32.0.tile3.nc,C256_grid_32.0.tile4.nc,C256_grid_32.0.tile5.nc,C256_grid_32.0.tile6.nc
  [ "$status" -eq 0 ]

  run command make_solo_mosaic \
                --num_tiles 6 \
                --dir ./ \
                --mosaic_name C256_mosaic_34.0 --tile_file C256_grid_34.0.tile1.nc,C256_grid_34.0.tile2.nc,C256_grid_34.0.tile3.nc,C256_grid_34.0.tile4.nc,C256_grid_34.0.tile5.nc,C256_grid_34.0.tile6.nc
  [ "$status" -eq 0 ]

  run command make_solo_mosaic \
                --num_tiles 6 \
                --dir ./ \
                --mosaic_name C256_mosaic_35.4 --tile_file C256_grid_35.4.tile1.nc,C256_grid_35.4.tile2.nc,C256_grid_35.4.tile3.nc,C256_grid_35.4.tile4.nc,C256_grid_35.4.tile5.nc,C256_grid_35.4.tile6.nc
  [ "$status" -eq 0 ]

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
  var_32_0=$output | awk '{ print sprintf("%.9f", $1); }'
  echo $output
  [[ ${var_32_0} < 0.0 ]]

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
  var_34_0=$output | awk '{ print sprintf("%.9f", $1); }'
  echo $output
  [[ ${var_34_0} < 0.0 ]]

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
  var_35_4=$output | awk '{ print sprintf("%.9f", $1); }'
  echo $output
  [[ ${var_35_4} < 0.0 ]]
 
  cd ..
#  rm -rf Test28
}
