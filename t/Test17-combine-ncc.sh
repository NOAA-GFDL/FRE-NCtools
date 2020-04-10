#!/usr/bin/env bats

@test "combine-ncc combines compressed netcdf files" {
  if [ ! -d "Test17" ] 
  then
  mkdir Test17
  fi

  cd Test17

  #ncks atmos daily *.nc file
 
  ncks $top_srcdir/t/Test17-input/combine-ncc.atmos_daily.nc combine-ncc.atmos_daily.nc.copy

  #Combine netcdf copy file 
  run command combine-ncc \
      combine-ncc.atmos_daily.nc.copy \
      combine-ncc_output.nc 
  [ "$status" -eq 0 ]
  [ -e combine-ncc_output.nc ]
  run ncdump -h combine-ncc_output.nc
  [ "$status" -eq 0 ]

  #Clean up 
  cd ..
  rm -rf Test17
}
