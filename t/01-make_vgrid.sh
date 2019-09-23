#!/usr/bin/env bats
teardown () {
  echo "$output"
  rm -f ocean_vgrid.nc 
}

@test "make_vgrid exists and is executable" {
  run command -v make_vgrid
  [ "$status" -eq 0 ]
  run make_vgrid -h
  [ "$status" -eq 1 ]
}

@test "make_vgrid creates ocean_vgrid" {
  run command make_vgrid \
      --nbnds 3 \
      --bnds 0.,220.,5500. \
      --dbnds 10.,10.,367.14286 \
      --center c_cell \
      --grid_name ocean_vgrid
  [ "$status" -eq 0 ]
}
