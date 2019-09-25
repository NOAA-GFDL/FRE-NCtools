#!/usr/bin/env bats
# Test remap runoff data from regular lat-lon grid onto cm2m grid

@test "Test  remap runoff data from regular lat-lon grid onto cm2m grid" {

mkdir work_dir_5
cd work_dir_5
cp $top_srcdir/t/runoff_regrid/* .

  run command runoff_regrid \
		--input_file runoff.daitren.iaf.nc \
		--input_fld_name runoff \
		--output_mosaic ocean_mosaic.nc \
		--output_topog topog.nc \
		--output_file runoff.cm2m.nc
  [ "$status" -eq 0 ]

cd ..
rm -rf work_dir_5
}
