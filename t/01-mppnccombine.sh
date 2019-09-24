#!/usr/bin/env bats

@test "mppnccombine combines" {
  mkdir work_dir_2
  cd work_dir_2

  #Generate the netcdf files from the .ncl 
  for f in $top_srcdir/t/mppnccombine/*.ncl.????
  do
    ncgen -o mppnccombine.nc.${f##*.} $f
  done

  run command mppnccombine \
      mppnccombine_output.nc \
      mppnccombine.nc.????
  [ "$status" -eq 0 ]
  [ -e mppnccombine_output.nc ]
  run ncdump -h mppnccombine_output.nc
  [ "$status" -eq 0 ]

  cd ..
  rm -rf work_dir_2
}
