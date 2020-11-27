#!/usr/bin/env bats

@test "Wrapper complete hydrology test" {

  # Not all systems have /bin/csh, skip if it doesn't exist
  test ! -e /bin/csh && skip 'System does not have /bin/csh'
  if [ ! -d "Test25" ]
  then
                mkdir Test25
  fi

  cd Test25

  run $top_srcdir/tools/simple_hydrog/share/make_simple_hydrog.csh
  [ "$status" -eq 1 ]

  mkdir tmpdir
  export TMPDIR=tmpdir

  #Run wrapper hydrology script
  run $top_srcdir/tools/simple_hydrog/share/make_simple_hydrog.csh -f 0. -t 1.e-5 -m $top_srcdir/t/Test25-input/grid_spec.nc
  [ "$status" -eq 0 ]

  #Clean up
  cd ..
  rm -rf Test25

}
