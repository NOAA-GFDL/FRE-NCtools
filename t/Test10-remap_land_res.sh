#!/usr/bin/env bats
# Test remap_land: remap land restart files.

@test "Test remap_land: remap land restart files" {
  if [ ! -d "Test10" ] 
  then
  		mkdir Test10
  fi

  cd Test10
  mkdir src_restart
  cd src_restart
  cp $top_srcdir/t/Test10-input/src_restart/* . 
  files=`ls`
  for file in $files
  do
	 ncgen -o ${file:0:-1} ${file}
  done
  cd ..

  mkdir dst_cold_restart
  cd dst_cold_restart
  cp $top_srcdir/t/Test10-input/dst_cold_restart/* . 
  files=`ls`
  for file in $files
  do
	 ncgen -o ${file:0:-1} ${file}
  done
  cd ..

  mkdir C48_mosaic
  cd C48_mosaic
  cp $top_srcdir/t/Test10-input/C48_mosaic/* . 
  files=`ls`
  for file in $files
  do
	 ncgen -o ${file:0:-1} ${file}
  done

  cd ..
  mkdir C192_mosaic
  cd C192_mosaic
  cp $top_srcdir/t/Test10-input/C192_mosaic/* . 
  files=`ls`
  for file in $files
  do
	 ncgen -o ${file:0:-1} ${file}
  done
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


#TO DO: Commenting this out because it requires too many files 
:'
#TO DO: Try this in a loop? 

#Try to remap the soil restart
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

#Try to remap the snow restart
  run command remap_land \
		--file_type snow \
		--src_mosaic C48_mosaic/C48_mosaic.nc \
		--dst_mosaic C192_mosaic/C192_mosaic.nc \
		--src_restart src_restart/snow.res \
		--dst_restart snow.res \
		--dst_cold_restart dst_cold_restart/snow.res \
		--land_src_restart src_restart/land.res \
		--land_cold_restart dst_cold_restart/land.res \
		--remap_file remap_file_C48_to_C192 --print_memory
  [ "$status" -eq 0 ] 

#Try to remap the cana restart
  run command remap_land \
		--file_type cana \
		--src_mosaic C48_mosaic/C48_mosaic.nc \
		--dst_mosaic C192_mosaic/C192_mosaic.nc \
		--src_restart src_restart/cana.res \
		--dst_restart cana.res \
		--dst_cold_restart dst_cold_restart/cana.res \
		--land_src_restart src_restart/land.res \
		--land_cold_restart dst_cold_restart/land.res \
		--remap_file remap_file_C48_to_C192 --print_memory
  [ "$status" -eq 0 ] 

#Try to remap the glac restart
  run command remap_land \
		--file_type glac \
		--src_mosaic C48_mosaic/C48_mosaic.nc \
		--dst_mosaic C192_mosaic/C192_mosaic.nc \
		--src_restart src_restart/glac.res \
		--dst_restart glac.res \
		--dst_cold_restart dst_cold_restart/glac.res \
		--land_src_restart src_restart/land.res \
		--land_cold_restart dst_cold_restart/land.res \
		--remap_file remap_file_C48_to_C192 --print_memory
  [ "$status" -eq 0 ] 

#Try to remap the lake restart
  run command remap_land \
		--file_type lake \
		--src_mosaic C48_mosaic/C48_mosaic.nc \
		--dst_mosaic C192_mosaic/C192_mosaic.nc \
		--src_restart src_restart/lake.res \
		--dst_restart lake.res \
		--dst_cold_restart dst_cold_restart/lake.res \
		--land_src_restart src_restart/land.res \
		--land_cold_restart dst_cold_restart/land.res \
		--remap_file remap_file_C48_to_C192 --print_memory
  [ "$status" -eq 0 ] 

#Try to remap the vegn1 restart
  run command remap_land \
		--file_type vegn \
		--src_mosaic C48_mosaic/C48_mosaic.nc \
		--dst_mosaic C192_mosaic/C192_mosaic.nc \
		--src_restart src_restart/vegn1.res \
		--dst_restart vegn1.res \
		--dst_cold_restart dst_cold_restart/vegn1.res \
		--land_src_restart src_restart/land.res \
		--land_cold_restart dst_cold_restart/land.res \
		--remap_file remap_file_C48_to_C192 --print_memory
  [ "$status" -eq 0 ] 

#Try to remap the vegn2 restart
  run command remap_land \
		--file_type vegn \
		--src_mosaic C48_mosaic/C48_mosaic.nc \
		--dst_mosaic C192_mosaic/C192_mosaic.nc \
		--src_restart src_restart/vegn2.res \
		--dst_restart vegn2.res \
		--dst_cold_restart dst_cold_restart/vegn2.res \
		--land_src_restart src_restart/land.res \
		--land_cold_restart dst_cold_restart/land.res \
		--remap_file remap_file_C48_to_C192 --print_memory
  [ "$status" -eq 0 ] 
'
  cd ..
  rm -rf Test10

}
