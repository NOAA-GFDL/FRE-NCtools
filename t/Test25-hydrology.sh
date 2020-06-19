#!/usr/bin/env bats

@test "Wrapper complete hydrology test" {

  if [ ! -d "Test25" ]
  then
                mkdir Test25
  fi

  cd Test25

  #Run wrapper hydrology script
  $top_srcdir/tools/simple_hydrog/make_simple_hydrog.csh -f 0. -t 1.e-5 -m $top_srcdir/t/Test25-input/grid_spec.nc
  status=$(echo "$?")
  [ "$status" -eq 0 ]
  [ -e lake_frac.tile6.nc ]
  run ncdump -h lake_frac.tile6.nc
  [ "$status" -eq 0 ]
  [ -e  hydrography.tile6.nc ]
  run ncdump -h  hydrography.tile6.nc
  [ "$status" -eq 0 ]

  #Clean up
  cd ..
  rm -rf Test25

}

