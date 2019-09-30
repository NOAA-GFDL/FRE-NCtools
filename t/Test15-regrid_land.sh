#!/usr/bin/env bats
# Test regrid land data with cell_measures and cell_methods attribute

@test "Test regrid land data" {

  if [ ! -d "Test15" ] 
  then
  		mkdir Test15
  fi

  cd Test15
  cp $top_srcdir/t/Test15-input/* .

  run command fregrid \
		--input_mosaic C180_mosaic.nc \
		--interp_method conserve_order1 \
		--nlon 144 \
		--nlat 90 \
		--remap_file remap_file.nc
  [ "$status" -eq 0 ] 

#remap static field
  run command fregrid  \
		--input_mosaic C180_mosaic.nc  \
		--interp_method conserve_order1  \
		--nlon 144  \
		--nlat 90  \
		--input_file 00050101.land_static  \
		--scalar_field soil_frac,lake_frac,glac_frac,area,soil_area,lake_area,glac_area  \
		--output_file out.nc  \
		--remap_file remap_file.nc 
  [ "$status" -eq 0 ] 

#TO DO: Add a parallel call and use nccmp to compare

# remap other fields
# Commented this part out because the input file is too large
#  run command fregrid  \
#		--input_mosaic C180_mosaic.nc  \
#		--interp_method conserve_order1  \
#		--nlon 144  \
#		--nlat 90  \
#		--input_file 00050101.land_month  \
#		--scalar_field evap_land,evap_soil,evap_glac,evap_lake  \
#		--output_file 00050101.land_month.nc  \
#		--remap_file remap_file.nc
#  [ "$status" -eq 0 ] 

  cd ..
  rm -rf Test15
}
