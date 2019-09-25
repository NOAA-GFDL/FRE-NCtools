#!/usr/bin/env bats
# Test check_mask: create mask_table to mask out some all-land domain
# to save processor usage for a sea-ice model, baltic1 experiment

@test "Test check_ mask for baltic1 experiment" {

mkdir work_dir_7
cd work_dir_7
cp $top_srcdir/t/checkmask_baltic1/* .

  run command check_mask \
		--grid_file baltic1_grid_spec.nc \
		--min_pe 60 \
		--max_pe 80 
  [ "$status" -eq 0 ]

cd ..
rm -rf work_dir_7

}
