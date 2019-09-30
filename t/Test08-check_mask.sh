#!/usr/bin/env bats
# Test check_mask: create mask_table to mask out some all-land domain
# to save processor usage for a sea-ice model, baltic1 experiment

@test "Test check_ mask for baltic1 experiment" {

  if [ ! -d "Test08" ] 
  then
  		mkdir Test08
  fi

  cd Test08
  cp $top_srcdir/t/Test08-input/* .

  run command check_mask \
		--grid_file baltic1_grid_spec.nc \
		--min_pe 60 \
		--max_pe 80 
  [ "$status" -eq 0 ]

  cd ..
  rm -rf Test08

}
