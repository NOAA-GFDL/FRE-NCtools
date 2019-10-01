#!/usr/bin/env bats
# Test remap runoff data from regular lat-lon grid onto cm2m grid

@test "Test  remap runoff data from regular lat-lon grid onto cm2m grid" {

  if [ ! -d "Test07" ] 
  then
  		mkdir Test07
  fi

  cd Test07
  cp $top_srcdir/t/Test06-input/ocean_hgrid.nc ocean_hgrid.nc
  cp $top_srcdir/t/Test06-input/ocean_mosaic.nc ocean_mosaic.nc
  cp $top_srcdir/t/Test06-input/ocean_vgrid.nc ocean_vgrid.nc
  cp $top_srcdir/t/Test03-input/C48_grid.tile1.nc C48_grid.tile1.nc
  cp $top_srcdir/t/Test03-input/C48_grid.tile2.nc C48_grid.tile2.nc
  cp $top_srcdir/t/Test03-input/C48_grid.tile3.nc C48_grid.tile3.nc
  cp $top_srcdir/t/Test03-input/C48_grid.tile4.nc C48_grid.tile4.nc
  cp $top_srcdir/t/Test03-input/C48_grid.tile5.nc C48_grid.tile5.nc
  cp $top_srcdir/t/Test03-input/C48_grid.tile6.nc C48_grid.tile6.nc
  cp $top_srcdir/t/Test03-input/C48_mosaic.nc C48_mosaic.nc
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
