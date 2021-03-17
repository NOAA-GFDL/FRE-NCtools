#!/usr/bin/env bats
# test stretched grid 

@test "Test no stretched grid data lat 34.0" {

  if [ ! -d "Test36" ] 
  then
  		mkdir Test36
  fi

  cd Test36
  cp $top_srcdir/t/Test32-input/ocean_hgrid.nc . 
  cp $top_srcdir/t/Test32-input/ocean_mosaic.nc .
  cp $top_srcdir/t/Test32-input/topog.nc .
  cp $top_srcdir/t/Test32-input/fv_tracer.res.tile*.nc .

#Make no streetched grid 
  run command make_hgrid \
               --grid_type gnomonic_ed \
               --target_lon 262.4 \
               --target_lat 34.0 \
               --nlon 512 \
               --grid_name C384_grid
  [ "$status" -eq 0 ]

#Create no stretched grid mosaic
  run command make_solo_mosaic \
                --num_tiles 6 \
                --dir ./ \
                --mosaic_name C384_mosaic --tile_file C384_grid.tile1.nc,C384_grid.tile2.nc,C384_grid.tile3.nc,C384_grid.tile4.nc,C384_grid.tile5.nc,C384_grid.tile6.nc
  [ "$status" -eq 0 ]

# no stretched grid lat 34.0
  result="$(mpirun -np 16 fregrid \
                --input_mosaic C384_mosaic.nc \
                --input_file fv_tracer.res \
                --scalar_field o3 \
                --nlon 640 \
                --nlat 400 \
                --latBegin 15.0 \
                --latEnd 65.0 \
                --lonBegin 230.0 \
                --lonEnd 310.0 \
                --remap_file fregrid_remap_file_640_by_400.nc \
                --output_file out.nc \
                --check_conserve)"
#  declare -f conserve_interp_area1 = echo $result | head -n 5 | rev | cut -c20-30 | rev
#  declare -f conserve_interp_area2 = echo $result | head -n 5 | rev | cut -c1-11 | rev
#  declare -f diff_dat_area = echo $result | head -n 6 | rev | cut -c1-2 | rev
  declare -f conserve_ratio = echo $result | head -n 5 | rev | cut -c39-49 | rev

#  [[ ${diff_dat_area} < 1. ]]
  [[ ${conserve_ratio} < 1. ]]
  cd ..
#  rm -rf Test36

}
