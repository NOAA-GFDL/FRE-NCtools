#!/usr/bin/env bats
# Test make_regional_mosaic: create mosaic file for regional output and do remapping
#             for regional output

@test "Test make_regional_mosaic" {
  skip "TO DO: Files are too large" 

  if [ ! -d "Test11" ] 
  then
  		mkdir Test11
  fi

  cd Test11
  cp $top_srcdir/t/Test11-input/* . 

  ncgen -o rregionatmos_4xdaily_eq_avg.tile2.nc rregionatmos_4xdaily_eq_avg.tile2.cdl
  ncgen -o rregionatmos_4xdaily_eq_avg.tile3.nc rregionatmos_4xdaily_eq_avg.tile3.cdl
  ncgen -o rregionatmos_4xdaily_eq_avg.tile5.nc rregionatmos_4xdaily_eq_avg.tile5.cdl
  ncgen -o rregionatmos_4xdaily_eq_avg.tile6.nc rregionatmos_4xdaily_eq_avg.tile6.cdl

  run command make_regional_mosaic \
		--global_mosaic C192-rot-120.-0.-r1_mosaic.nc \
		--regional_file rregionatmos_4xdaily_eq_avg.tile2.nc
  [ "$status" -eq 0 ] 

  run command make_regional_mosaic \
		--global_mosaic C192-rot-120.-0.-r1_mosaic.nc \
		--regional_file rregionatmos_4xdaily_eq_avg.tile3.nc
  [ "$status" -eq 0 ] 

  run command make_regional_mosaic  \
		--global_mosaic C192-rot-120.-0.-r1_mosaic.nc  \
		--regional_file rregionatmos_4xdaily_eq_avg.tile5.nc
  [ "$status" -eq 0 ] 

  run command make_regional_mosaic  \
		--global_mosaic C192-rot-120.-0.-r1_mosaic.nc  \
		--regional_file rregionatmos_4xdaily_eq_avg.tile6.nc
  [ "$status" -eq 0 ] 

#Create regional mosaic
  run command make_solo_mosaic  \
		--num_tiles 4  \
		--dir ./  \
		--mosaic regional_mosaic  \
		--tile_file regional_grid.tile2.nc,regional_grid.tile3.nc,regional_grid.tile5.nc,regional_grid.tile6.nc
  [ "$status" -eq 0 ] 

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

#Remap data
  run command fregrid  \
		--input_mosaic regional_mosaic.nc  \
		--output_mosaic latlon_mosaic.nc  \
		--input_file rregionatmos_4xdaily_eq_avg  \
		--output_file atmos_4xdaily_latlon.nc  \
		--scalar_field ps  \
		--check_conserve
  [ "$status" -eq 0 ] 

  cd ..
  rm -rf Test11

}


