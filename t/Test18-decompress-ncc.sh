#!/usr/bin/env bats

@test "decompress input netcdf files" {
  if [ ! -d "Test18" ] 
  then
  mkdir Test18
  fi

  cd Test18

  cp $top_srcdir/t/Test18-input/decompress-ncc.atmos_daily.nc.copy .

  #Decompress compressed netcdf file(s) into 1
  run command decompress-ncc \
      decompress-ncc.atmos_daily.nc.copy \
      decompress-ncc_output.nc 
  [ "$status" -eq 0 ]
  [ -e decompress-ncc_output.nc ]
  run ncdump -h decompress-ncc_output.nc
  [ "$status" -eq 0 ]

  #Clean up 
  cd ..
  rm -rf Test18
}
