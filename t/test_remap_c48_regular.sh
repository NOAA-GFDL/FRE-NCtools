#!/usr/bin/env bats

@test "remap data from C48 to regular lat-lon grid" {
  mkdir work_dir_3
  cd work_dir_3
  cp $top_srcdir/t/grid_coupled_model/* .

  mkdir output

  run command fregrid \
		--input_mosaic C48_mosaic.nc \
		--input_file 19800101.atmos_daily \
		--scalar_field zsurf,temp,t_surf \
		--nlon 144 \
		--nlat 90 \
		--interp_method conserve_order2 \
		--output_dir output \
		--output_file 19800101.atmos_daily.nc \
		--check_conserve \
		--remap_file C48_to_N45_remap.nc 
  [ "$status" -eq 0 ]

  cd ..
  rm -rf work_dir_3
}

