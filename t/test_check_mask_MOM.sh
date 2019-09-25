#!/usr/bin/env bats
# Test check_mask: create mask_table to mask out some all-land domain
# to save processor usage for a coupled model with cm2m grid.

@test "Test check_ mask for baltic1 experiment" {
skip "Skip this for now because the input files is too big" 

mkdir work_dir_8
cd work_dir_8
cp $top_srcdir/t/checkmask_MOM6/* .

  run command check_mask \
		--grid_file ocean_mosaic.nc \
		--ocean_topog topog.nc \
		--layout 45,72
  [ "$status" -eq 0 ]

cd ..
rm -rf work_dir_8

}
