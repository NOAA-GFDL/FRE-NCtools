#!/usr/bin/env bats
teardown () {
  echo "$output"
  rm -f ocean_hgrid.nc
}

@test "make_hgrid exists and is executable" {
  run command -v make_hgrid
  [ "$status" -eq 0 ]
  run make_hgrid -h
  [ "$status" -eq 2 ]
}

@test "make_hgrid creates ocean_hgrid" {
  run command make_hgrid \
      --grid_type tripolar_grid \
      --nxbnd 2 \
      --nybnd 7 \
      --xbnd -280,80 \
      --ybnd -82,-30,-10,0,10,30,90 \
      --dlon 1.0,1.0 \
      --dlat 1.0,1.0,0.6666667,0.3333333,0.6666667,1.0,1.0 \
      --grid_name ocean_hgrid \
      --center c_cell
  [ "$status" -eq 0 ]
}
