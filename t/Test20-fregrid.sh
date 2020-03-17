#!/usr/bin/env bats
# test regrid ocean restart file 

@test "Test remap ocean restart file" {
  skip "TO DO: the input files needed are too large" 

  if [ ! -d "Test20" ] 
  then
  		mkdir Test20
  fi

  cd Test20
  cp /archive/fms/tools/input/CM2.1_to_CM2.5/CM2.1_mosaic.nc .
  cp /archive/fms/tools/input/CM2.1_to_CM2.5/ocean_temp_salt.res.nc . 

#Remap data from CM2.1 ocean grid onto regular lat-lon grid.
  run command fregrid   \
		--input_mosaic CM2.1_mosaic.nc   \
		--input_file ocean_temp_salt.res.nc   \
		--scalar_field temp,salt  \
		--output_file ocean_temp_salt.res.latlon.nc   \
		--output_mosaic latlon_mosaic.nc   \
		--check_conserve
  [ "$status" -eq 0 ]

  cd ..
  rm -rf Test20
}
