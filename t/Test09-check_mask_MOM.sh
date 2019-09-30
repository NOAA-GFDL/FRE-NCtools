#!/usr/bin/env bats
# Test check_mask for MOM6 grid 

@test "Test check_mask for MOM6 grid " {
skip "Skip this for now because the input files is too big" 

  if [ ! -d "Test09" ] 
  then
  		mkdir Test09
  fi

  cd Test09
  cp $top_srcdir/t/Test09-input/* .

  run command check_mask \
		--grid_file ocean_mosaic.nc \
		--ocean_topog topog.nc \
		--layout 45,72
  [ "$status" -eq 0 ]

  cd ..
  rm -rf Test09

}
