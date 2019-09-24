#!/usr/bin/env bats

@test "Test 01: Check make_hgrid exists and is executable" {
  run command -v make_hgrid
  [ "$status" -eq 0 ]
  run make_hgrid -h
  [ "$status" -eq 2 ]
}

@test "Test 02: Check make_vgrid exists and is executable" {
  run command -v make_vgrid
  [ "$status" -eq 0 ]
  run make_vgrid -h
  [ "$status" -eq 1 ]
}

@test "Test 03: Check make_solo_mosaic exists and is executable" {
  run command -v make_solo_mosaic
  [ "$status" -eq 0 ]
  run make_solo_mosaic -h
  [ "$status" -eq 2 ]
}

@test "Test 04: Check make_topog exists and is executable" {
  run command -v make_topog
  [ "$status" -eq 0 ]
  run make_topog -h
  [ "$status" -eq 1 ]
}

@test "Test 05: Check make_topog_parallel exists and is executable" {
  skip "The test fails on travis"
  run command -v make_topog_parallel
  [ "$status" -eq 0 ]
}

@test "Test 06: Check coupler_mosaic exists and is executable" {
  run command -v make_coupler_mosaic
  [ "$status" -eq 0 ]
  run make_coupler_mosaic -h
  [ "$status" -eq 2 ]
}

@test "Test 07: Check coupler_mosaic_parallel exists and is executable" {
  skip "The test fails on travis"
  run command -v make_coupler_mosaic_parallel
  [ "$status" -eq 0 ]
}

@test "Test 08: Check fregrid exists and is executable" {
  run command -v fregrid
  [ "$status" -eq 0 ]
  run fregrid -h
  [ "$status" -eq 2 ]
}

@test "Test 09: Check fregrid_parallel exists and is executable" {
  skip "The test fails on travis"
  run command -v fregrid_parallel
  [ "$status" -eq 0 ]
}

@test "Test 10: Check runoff_regrid exists and is executable" {
  run command -v runoff_regrid
  [ "$status" -eq 0 ]
  run runoff_regrid -h
  [ "$status" -eq 2 ]
}

@test "Test 11: Check river_regrid exists and is executable" {
  run command -v river_regrid
  [ "$status" -eq 0 ]
  run river_regrid -h
  [ "$status" -eq 2 ]
}

@test "Test 12: Check check_mask exists and is executable" {
  run command -v check_mask
  [ "$status" -eq 0 ]
  run check_mask -h
  [ "$status" -eq 2 ]
}

@test "Test 13: Check remap_land exists and is executable" {
  run command -v remap_land
  [ "$status" -eq 0 ]
  run remap_land -h
  [ "$status" -eq 1 ]
}

@test "Test 14: Check remap_land_parallel exists and is executable" {
  skip "The test fails on travis"
  run command -v remap_land_parallel
  [ "$status" -eq 0 ]
}

@test "Test 15: Check make_regional_mosaic exists and is executable" {
  run command -v make_regional_mosaic
  [ "$status" -eq 0 ]
  run make_regional_mosaic -h
  [ "$status" -eq 2 ]
}

@test "Test 16: Check mppncscatter exists and is executable" {
  run command -v mppncscatter
  [ "$status" -eq 0 ]

}

