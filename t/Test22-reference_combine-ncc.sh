#!/usr/bin/env bats

@test "reference combine-ncc combines comparison to bronx-16 stored copy" {
  if [ ! -d "Test22" ] 
  then
  mkdir Test22
  fi

  cd Test22

  #ncks atmos daily *.nc file
 
  cp $top_srcdir/t/Test17-input/combine-ncc.atmos_daily.nc.copy .

  #Combine netcdf copy file 
  run command combine-ncc \
      combine-ncc.atmos_daily.nc.copy \
      combine-ncc_output.nc 
  [ "$status" -eq 0 ]
  [ -e combine-ncc_output.nc ]
  run ncdump -h combine-ncc_output.nc
  [ "$status" -eq 0 ]

  run nccmp -d combine-ncc_output.nc  $top_srcdir/t/Test17-reference/combine-ncc_output.nc
  [ "$status" -eq 0 ]

  #Clean up 
  cd ..
  rm -rf Test22
}
