#!/usr/bin/env bats
# test remap ocean restart file from CM2.1 to CM2.5

@test "Test remap ocean restart file" {
  skip "Skipping this for now because the input files needed are too large" 

  if [ ! -d "Test14" ] 
  then
  		mkdir Test14
  fi

  cd Test14
  cp $top_srcdir/t/Test14-input/* . 

#Create a OM3-like regular lat-lon grid.
  run command make_hgrid \
		--grid_type regular_lonlat_grid \
		--nxbnd 2 \
		--nybnd 7 \
		--xbnd -280,80 \
		--ybnd -82,-30,-10,0,10,30,90 \
		--dlon 1.0,1.0  \
		--dlat 1.0,1.0,0.6666667,0.3333333,0.6666667,1.0,1.0  \
		--grid_name latlon_grid  \
		--center c_cell
  [ "$status" -eq 0 ]

#Create the lat-lon solo mosaic
  run command make_solo_mosaic  \
		--num_tiles 1  \
		--dir ./  \
		--mosaic_name latlon_mosaic  \
		--tile_file latlon_grid.nc  \
		--periodx 360
  [ "$status" -eq 0 ]

#Remap data from CM2.1 ocean grid onto regular lat-lon grid.
  run command fregrid   \
		--input_mosaic CM2.1_mosaic.nc   \
		--input_file ocean_temp_salt.res.nc   \
		--scalar_field temp,salt  \
		--output_file ocean_temp_salt.res.latlon.nc   \
		--output_mosaic latlon_mosaic.nc   \
		--check_conserve
  [ "$status" -eq 0 ]

#TO DO: This fails, it can probably work with npes>1 and in parallel 
  run command fregrid    \
		--input_mosaic latlon_mosaic.nc    \
		--input_dir ./   \
		--input_file ocean_temp_salt.res.latlon.nc    \
		--scalar_field temp,salt   \
		--output_file ocean_temp_salt.res.CM2.5.nc    \
		--output_mosaic CM2.5_mosaic.nc    \
		--extrapolate     \
		--check_conserve 
  [ "$status" -eq 0 ]

  cd ..
  rm -rf Test14
}
