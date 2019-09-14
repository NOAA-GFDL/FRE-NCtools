#!/usr/bin/env bats
setup () {
  for f in $top_srcdir/t/mppnccombine/*.ncl.????
  do
    ncgen -o mppnccombine.nc.${f##*.} $f
  done
}

teardown () {
  rm -f *.nc *.nc.????
}

@test "mppnccombine exists and is executable" {
  run command -v mppnccombine
  [ "$status" -eq 0 ]
  run mppnccombine -h
  [ "$status" -eq 1 ]
}

@test "mppnccombine combines" {
  run command mppnccombine \
      mppnccombine_output.nc \
      mppnccombine.nc.????
  [ "$status" -eq 0 ]
  [ -e mppnccombine_output.nc ]
  run ncdump -h mppnccombine_output.nc
  [ "$status" -eq 0 ]
}
