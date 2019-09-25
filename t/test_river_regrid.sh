#!/usr/bin/env bats
# Test river_regrid: remap data from lat-lon onto C48 grid

@test "Test  remap runoff data from regular lat-lon grid onto cm2m grid" {

mkdir work_dir_6
cd work_dir_6
cp $top_srcdir/t/river_regrid/* .

  run command river_regrid \
		--mosaic grid_spec.nc \
		--river_src z1l_river_output_M45_tripolar_aug24.nc \
		--output river_data_C48 \
		--land_thresh 0.000001 \
  [ "$status" -eq 0 ]

cd ..
rm -rf work_dir_6
}
