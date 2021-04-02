#!/usr/bin/env bats
# test no stretched grid 

@test "Test no stretched grid data lats 32.0 34.0 35.4" {

  if [ ! -d "Test29" ] 
  then
  		mkdir Test29
  fi

  cd Test29
  cp $top_srcdir/t/Test28-input/ocean_hgrid.nc . 
  cp $top_srcdir/t/Test28-input/ocean_mosaic.nc .
  cp $top_srcdir/t/Test28-input/topog.nc .

  for i in {1..6}
  do
    $top_srcdir/t/ncgenerator.py out.tile"$i".nc
  done

#Make no stretched grid 
  run command make_hgrid \
               --grid_type gnomonic_ed \
               --target_lon 262.4 \
               --target_lat 32.0 \
               --nlon 512 \
               --grid_name C384_grid_32.0
  [ "$status" -eq 0 ]

  run command make_hgrid \
               --grid_type gnomonic_ed \
               --target_lon 262.4 \
               --target_lat 32.0 \
               --nlon 512 \
               --grid_name C384_grid_34.0
  [ "$status" -eq 0 ]

  run command make_hgrid \
               --grid_type gnomonic_ed \
               --target_lon 262.4 \
               --target_lat 32.0 \
               --nlon 512 \
               --grid_name C384_grid_35.4
  [ "$status" -eq 0 ]

#Create no stretched grid mosaic
  run command make_solo_mosaic \
                --num_tiles 6 \
                --dir ./ \
                --mosaic_name C384_mosaic_32.0 --tile_file C384_grid_32.0.tile1.nc,C384_grid_32.0.tile2.nc,C384_grid_32.0.tile3.nc,C384_grid_32.0.tile4.nc,C384_grid_32.0.tile5.nc,C384_grid_32.0.tile6.nc
  [ "$status" -eq 0 ]

  run command make_solo_mosaic \
                --num_tiles 6 \
                --dir ./ \
                --mosaic_name C384_mosaic_34.0 --tile_file C384_grid_34.0.tile1.nc,C384_grid_34.0.tile2.nc,C384_grid_34.0.tile3.nc,C384_grid_34.0.tile4.nc,C384_grid_34.0.tile5.nc,C384_grid_34.0.tile6.nc
  [ "$status" -eq 0 ]

  run command make_solo_mosaic \
                --num_tiles 6 \
                --dir ./ \
                --mosaic_name C384_mosaic_35.4 --tile_file C384_grid_35.4.tile1.nc,C384_grid_35.4.tile2.nc,C384_grid_35.4.tile3.nc,C384_grid_35.4.tile4.nc,C384_grid_35.4.tile5.nc,C384_grid_35.4.tile6.nc
  [ "$status" -eq 0 ]

# no stretched grid lats 32.0 34.0 35.4
  result_32_0="$(fregrid \
                --input_mosaic C384_mosaic_32.0.nc \
                --input_file out \
                --scalar_field o3 \
                --nlon 640 \
                --nlat 400 \
                --latBegin 15.0 \
                --latEnd 65.0 \
                --lonBegin 230.0 \
                --lonEnd 310.0 \
                --remap_file fregrid_remap_file_640_by_400_32.0.nc \
                --output_file out_32.0.nc --check_conserve)"
  declare -f conserve_ratio_32_0 = echo $result_32_0 | head -n 5 | rev | cut -c39-49 | rev

  result_34_0="$(fregrid \
                --input_mosaic C384_mosaic_34.0.nc \
                --input_file out \
                --scalar_field o3 \
                --nlon 640 \
                --nlat 400 \
                --latBegin 15.0 \
                --latEnd 65.0 \
                --lonBegin 230.0 \
                --lonEnd 310.0 \
                --remap_file fregrid_remap_file_640_by_400_34.0.nc \
                --output_file out_34.0.nc --check_conserve)"
  declare -f conserve_ratio_34_0 = echo $result_34_0 | head -n 5 | rev | cut -c39-49 | rev

  result_35_4="$(fregrid \
                --input_mosaic C384_mosaic_35.4.nc \
                --input_file out \
                --scalar_field o3 \
                --nlon 640 \
                --nlat 400 \
                --latBegin 15.0 \
                --latEnd 65.0 \
                --lonBegin 230.0 \
                --lonEnd 310.0 \
                --remap_file fregrid_remap_file_640_by_400_35.4.nc \
                --output_file out_35.4.nc --check_conserve)"
  declare -f conserve_ratio_35_4 = echo $result_35_4 | head -n 5 | rev | cut -c39-49 | rev

  [[ ${conserve_ratio_32_0} < 1. ]]
  [[ ${conserve_ratio_34_0} < 1. ]]
  [[ ${conserve_ratio_35_4} < 1. ]]

  cd ..
  rm -rf Test29
}
