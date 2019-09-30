#!/usr/bin/env bats
# Test remap data onto cm2m ocean grid with extrapolation and vertical interpolation

@test "Test remap data onto cm2m ocean grid with extrapolation and vertical interpolation" {

  if [ ! -d "Test06" ] 
  then
  		mkdir Test06
  fi

  cd Test06
  cp $top_srcdir/t/Test06-input/* .

  run command make_hgrid \
		--grid_type regular_lonlat_grid \
		--nxbnd 2 \
		--nybnd 2 \
		--xbnd 0,360 \
		--ybnd -90,90 \
		--nlon 720 \
		--nlat 360 \
		--grid_name levitus_grid
  [ "$status" -eq 0 ]

  run command make_solo_mosaic \
		--num_tiles 1 \
		--dir . \
		--mosaic_name levitus_mosaic \
		--tile_file levitus_grid.nc \
		--periodx 360
  [ "$status" -eq 0 ]

  run command fregrid \
		--input_mosaic levitus_mosaic.nc \
		--input_file WOA09_ann_theta.nc \
		--scalar_field POTENTIAL_TEMP \
		--output_file WOA09_ann_theta_cm2g_extrap.nc \
		--output_mosaic ocean_mosaic.nc \
		--extrapolate \
		--dst_vgrid ocean_vgrid.nc \
		--check_conserve
  [ "$status" -eq 0 ]

  cd ..
  rm -rf Test06

}
