#!/usr/bin/env bats
# test land mask stretched grid check area ratio 

@test "Test land mask stretched grid data lat 32.0 34.0 35.4" {

  if [ ! -d "Test30" ] 
  then
  		mkdir Test30
  fi

  cd Test30
  
  cp $top_srcdir/t/Test30-input/ocean_hgrid.nc . 
  cp $top_srcdir/t/Test30-input/ocean_mosaic.nc .
  cp $top_srcdir/t/Test30-input/ocean_topog.nc .

#Make no streetched grid 
  run command make_hgrid --grid_type gnomonic_ed \
                         --do_schmidt --stretch_factor 2.5 \
                         --target_lon 262.4 \
                         --target_lat 32.0 \
                         --nlon 512 \
                         --grid_name C256_grid_32.0
  [ "$status" -eq 0 ]

  run command make_hgrid --grid_type gnomonic_ed \
                         --do_schmidt --stretch_factor 2.5 \
                         --target_lon 262.4 \
                         --target_lat 34.0 \
                         --nlon 512 \
                         --grid_name C256_grid_34.0
  [ "$status" -eq 0 ]

  run command make_hgrid --grid_type gnomonic_ed \
                         --do_schmidt --stretch_factor 2.5 \
                         --target_lon 262.4 \
                         --target_lat 35.4 \
                         --nlon 512 \
                         --grid_name C256_grid_35.4
  [ "$status" -eq 0 ]

#Create no stretched grid mosaic
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

# land mask stretched grid lat 32.0 34.0 35.4 area ratio check
  run bash -c 'make_coupler_mosaic --atmos_mosaic C256_mosaic_32.0.nc \
                --ocean_mosaic ocean_mosaic.nc \
                --mosaic_name mosaic \
                --ocean_topog ocean_topog.nc \
                --verbose --check | tail -n 3 | head -n 1 | rev | cut -c2-9 | rev'
  tiling_error_32_0=$output | awk '{ print sprintf("%.9f", $1); }'
  echo "$output" 
  [[ ${tiling_error_32_0} < 0.01 ]]

  run bash -c 'make_coupler_mosaic --atmos_mosaic C256_mosaic_34.0.nc \
                --ocean_mosaic ocean_mosaic.nc \
                --mosaic_name mosaic \
                --ocean_topog ocean_topog.nc \
                --verbose --check | tail -n 3 | head -n 1 | rev | cut -c2-9 | rev'
  tiling_error_34_0=$output | awk '{ print sprintf("%.9f", $1); }'
  echo "$output"
  [[ ${tiling_error_34_0} < 0.01 ]]


  run bash -c 'make_coupler_mosaic --atmos_mosaic C256_mosaic_35.4.nc \
                --ocean_mosaic ocean_mosaic.nc \
                --mosaic_name mosaic \
                --ocean_topog ocean_topog.nc \
                --verbose --check | tail -n 3 | head -n 1 | rev | cut -c2-9 | rev'
  tiling_error_35_4=$output | awk '{ print sprintf("%.9f", $1); }'
  echo "$output"
  [[ ${tiling_error_35_4} < 0.01 ]]

  cd ..
  rm -rf Test30
}
