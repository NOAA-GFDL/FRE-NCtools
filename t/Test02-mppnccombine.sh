#!/usr/bin/env bats

@test "mppnccombine combines" {
  if [ ! -d "Test02" ] 
  then
  mkdir Test02
  fi

  cd Test02

  #Generate the netcdf files from the .ncl 
  for f in $top_srcdir/t/Test02-input/*.ncl.????
  do
    ncgen -o mppnccombine.nc.${f##*.} $f
  done

  #Combine the files into 1
  run command mppnccombine \
      mppnccombine_output.nc \
      mppnccombine.nc.????
  [ "$status" -eq 0 ]
  [ -e mppnccombine_output.nc ]
  run ncdump -h mppnccombine_output.nc
  [ "$status" -eq 0 ]

  #Clean up 
  cd ..
  rm -rf Test02
}
