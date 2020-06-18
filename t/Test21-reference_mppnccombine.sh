#!/usr/bin/env bats

@test " refernce mppnccombine combines comparison to bronx-16 stored copy" {
  if [ ! -d "Test21" ] 
  then
  mkdir Test21
  fi

  cd Test21

  #Generate the netcdf files from the .ncl 
  for f in $top_srcdir/t/Test02-input/*.ncl.????
  do
    ncgen -o mppnccombine.nc.${f##*.} $f
    [ -e mppnccombine.nc.${f##*.} ]
  done

  #Combine the files into 1
  run command mppnccombine \
      mppnccombine_output.nc \
      mppnccombine.nc.????
  [ "$status" -eq 0 ]
  [ -e mppnccombine_output.nc ]
  run ncdump -h mppnccombine_output.nc
  [ "$status" -eq 0 ]

  [ -e $top_srcdir/t/Test02-reference/mppnccombine_output.nc ]

  run nccmp -V
  [ "$status" -eq 0 ]

  run nccmp -d mppnccombine_output.nc  $top_srcdir/t/Test02-reference/mppnccombine_output.nc
  [ "$status" -eq 0 ]

  #Clean up 
  cd ..
  rm -rf Test21
}
