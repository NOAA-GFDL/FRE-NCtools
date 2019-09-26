#!/usr/bin/env bats
# Test remap_land: remap land restart files.

@test "Test remap_land: remap land restart files" {
  mkdir work_dir_9
  cd work_dir_9
  mkdir src_restart
  cd src_restart
  cp $top_srcdir/t/remap_land/src_restart/* . 
  cd ..
  mkdir dst_cold_restart
  cd dst_cold_restart
  cp $top_srcdir/t/remap_land/dst_cold_restart/* . 
  cd ..
  mkdir C48_mosaic
  cd C48_mosaic
  cp $top_srcdir/t/remap_land/C48_mosaic/* . 
  cd ..
  mkdir C192_mosaic
  cd C192_mosaic
  cp $top_srcdir/t/remap_land/C192_mosaic/* . 
  cd ..

  run command remap_land \
		--file_type land  \
		--src_mosaic C48_mosaic/C48_mosaic.nc \
		--dst_mosaic C192_mosaic/C192_mosaic.nc \
		--src_restart src_restart/land.res \
		--dst_restart land.res \
		--dst_cold_restart dst_cold_restart/land.res \
		--remap_file remap_file_C48_to_C192 --print_memory
  [ "$status" -eq 0 ]

#Try to remap the soil restart
  set filetype = soil
  set restart = soil 
  run command remap_land \
		--file_type soil \
		--src_mosaic C48_mosaic/C48_mosaic.nc \
		--dst_mosaic C192_mosaic/C192_mosaic.nc \
		--src_restart src_restart/soil.res \
		--dst_restart soil.res \
		--dst_cold_restart dst_cold_restart/soil.res \
		--land_src_restart src_restart/land.res \
		--land_cold_restart dst_cold_restart/land.res \
		--remap_file remap_file_C48_to_C192 --print_memory
  [ "$status" -eq 0 ] 

  cd ..
  rm -rf work_dir_9

}
