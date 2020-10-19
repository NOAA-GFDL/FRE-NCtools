#!/usr/bin/env bats

@test "Check make_hgrid exists and is executable" {
  run command -v make_hgrid
  [ "$status" -eq 0 ]
  run make_hgrid -h
  [ "$status" -eq 2 ]
}

@test "Check make_vgrid exists and is executable" {
  run command -v make_vgrid
  [ "$status" -eq 0 ]
  run make_vgrid -h
  [ "$status" -eq 1 ]
}

@test "Check make_solo_mosaic exists and is executable" {
  run command -v make_solo_mosaic
  [ "$status" -eq 0 ]
  run make_solo_mosaic -h
  [ "$status" -eq 2 ]
}

@test "Check make_topog exists and is executable" {
  run command -v make_topog
  [ "$status" -eq 0 ]
  run make_topog -h
  [ "$status" -eq 1 ]
}

@test "Check make_topog_parallel exists and is executable" {
  skip "TO DO: The test fails on travis"
  run command -v make_topog_parallel
  [ "$status" -eq 0 ]
}

@test "Check coupler_mosaic exists and is executable" {
  run command -v make_coupler_mosaic
  [ "$status" -eq 0 ]
  run make_coupler_mosaic -h
  [ "$status" -eq 2 ]
}

@test "Check coupler_mosaic_parallel exists and is executable" {
  skip "TO DO: The test fails on travis"
  run command -v make_coupler_mosaic_parallel
  [ "$status" -eq 0 ]
}

@test "Check fregrid exists and is executable" {
  run command -v fregrid
  [ "$status" -eq 0 ]
  run fregrid -h
  [ "$status" -eq 2 ]
}

@test "Check fregrid_parallel exists and is executable" {
  skip "TO DO: The test fails on travis"
  run command -v fregrid_parallel
  [ "$status" -eq 0 ]
}

@test "Check runoff_regrid exists and is executable" {
  run command -v runoff_regrid
  [ "$status" -eq 0 ]
  run runoff_regrid -h
  [ "$status" -eq 2 ]
}

@test "Check river_regrid exists and is executable" {
  run command -v river_regrid
  [ "$status" -eq 0 ]
  run river_regrid -h
  [ "$status" -eq 2 ]
}

@test "Check check_mask exists and is executable" {
  run command -v check_mask
  [ "$status" -eq 0 ]
  run check_mask -h
  [ "$status" -eq 2 ]
}

@test "Check remap_land exists and is executable" {
  run command -v remap_land
  [ "$status" -eq 0 ]
  run remap_land -h
  [ "$status" -eq 1 ]
}

@test "Check remap_land_parallel exists and is executable" {
  skip "TO DO: The test fails on travis"
  run command -v remap_land_parallel
  [ "$status" -eq 0 ]
}

@test "Check make_regional_mosaic exists and is executable" {
  run command -v make_regional_mosaic
  [ "$status" -eq 0 ]
  run make_regional_mosaic -h
  [ "$status" -eq 2 ]
}

@test "Check mppncscatter exists and is executable" {
  run command -v mppncscatter
  [ "$status" -eq 0 ]
}

@test "Check mppnccombine exists and is executable" {
  run command -v mppnccombine
  [ "$status" -eq 0 ]
  run mppnccombine -h
  [ "$status" -eq 1 ]
}

@test "Check combine-ncc exists and is executable" {
  run command -v combine-ncc
  [ "$status" -eq 0 ]
}

@test "Check decompress-ncc exists and is executable" {
  run command -v decompress-ncc
  [ "$status" -eq 0 ]
}

@test "Check cr_lake_files exists and is executable" {
     run command -v cr_lake_files
     [ "$status" -eq 0 ]
}

@test "Check cp_river_vars exists and is executable" {
     run command -v cp_river_vars
     [ "$status" -eq 0 ]
}

@test "Check rmv_parallel_rivers exists and is executable" {
     run command -v rmv_parallel_rivers
     [ "$status" -eq 0 ]
}

