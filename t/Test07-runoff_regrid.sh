#!/usr/bin/env bats
# Test remap runoff data from regular lat-lon grid onto cm2m grid

@test "Test  remap runoff data from regular lat-lon grid onto cm2m grid" {

  if [ ! -d "Test07" ] 
  then
  		mkdir Test07
  fi

  cd Test07
  cp $top_srcdir/t/Test07-input/* .

  run command runoff_regrid \
		--input_file runoff.daitren.iaf.nc \
		--input_fld_name runoff \
		--output_mosaic ocean_mosaic.nc \
		--output_topog topog.nc \
		--output_file runoff.cm2m.nc
  [ "$status" -eq 0 ]

  cd ..
  rm -rf Test07
}
