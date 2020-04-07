#!/usr/bin/env bats
# test regrid ocean restart file 

@test "Test fregrid ocean data" {

  if [ ! -d "Test24" ] 
  then
  		mkdir Test24
  fi

  cd Test24
  cp /archive/fms/tools/input/CM2.1_to_CM2.5/CM2.1_mosaic.nc .
  cp /archive/fms/tools/input/CM2.1_to_CM2.5/CM2.1_grid.nc .
  cp /archive/fms/tools/input/CM2.1_to_CM2.5/ocean_temp_salt.res.nc . 

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

#Remap data from CM2.1 ocean grid onto regular lat-lon grid.
  run command fregrid   \
		--input_mosaic CM2.1_mosaic.nc   \
		--input_file ocean_temp_salt.res.nc   \
		--scalar_field temp,salt  \
		--output_file ocean_temp_salt.res.latlon.nc   \
		--output_mosaic latlon_mosaic.nc   \
		--check_conserve
  [ "$status" -eq 0 ]

   run nccmp -d  ocean_temp_salt.res.latlon.nc  $top_srcdir/t/Test20-reference/ocean_temp_salt.res.latlon.nc
   [ "$status" -eq 0 ]
   run nccmp -d  latlon_mosaic.nc  $top_srcdir/t/Test20-reference/latlon_mosaic.nc
   [ "$status" -eq 0 ]

  cd ..
  rm -rf Test24
}
