#!/usr/bin/env bats

#***********************************************************************
#                   GNU Lesser General Public License
#
# This file is part of the GFDL FRE NetCDF tools package (FRE-NCTools).
#
# FRE-NCTools is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.
#
# FRE-NCTools is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
# for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with FRE-NCTools.  If not, see
# <http://www.gnu.org/licenses/>.
#***********************************************************************
load test_utils

@test "Check make_hgrid exists and is executable" {
  run_and_check -v make_hgrid
  run make_hgrid -h
  [ "$status" -eq 2 ]
}

@test "Check make_vgrid exists and is executable" {
  run_and_check -v make_vgrid
  run make_vgrid -h
  [ "$status" -eq 1 ]
}

@test "Check make_solo_mosaic exists and is executable" {
  run_and_check -v make_solo_mosaic
  run make_solo_mosaic -h
  [ "$status" -eq 2 ]
}

@test "Check make_topog exists and is executable" {
  run_and_check -v make_topog
  run make_topog -h
  [ "$status" -eq 1 ]
}
# only run with mpi
@test "Check make_topog_parallel exists and is executable" {
  [ ! -z $skip_mpi ] && skip "not built with MPI"
  run_and_check -v make_topog_parallel

}

@test "Check coupler_mosaic exists and is executable" {
  run_and_check -v make_coupler_mosaic
  run make_coupler_mosaic -h
  [ "$status" -eq 2 ]
}

# only run with mpi
@test "Check coupler_mosaic_parallel exists and is executable" {
  [ ! -z $skip_mpi ] && skip "not built with MPI"
  run_and_check -v make_coupler_mosaic_parallel

}

@test "Check fregrid exists and is executable" {
  run_and_check -v fregrid
  run fregrid -h
  [ "$status" -eq 2 ]
}

# only run with mpi
@test "Check fregrid_parallel exists and is executable" {
  [ ! -z $skip_mpi ] && skip "not built with MPI"
  run_and_check -v fregrid_parallel

}

@test "Check runoff_regrid exists and is executable" {
  run_and_check -v runoff_regrid
  run runoff_regrid -h
  [ "$status" -eq 2 ]
}

@test "Check river_regrid exists and is executable" {
  run_and_check -v river_regrid
  run river_regrid -h
  [ "$status" -eq 2 ]
}

@test "Check check_mask exists and is executable" {
  run_and_check -v check_mask
  run check_mask -h
  [ "$status" -eq 2 ]
}

@test "Check remap_land exists and is executable" {
  run_and_check -v remap_land
  run remap_land -h
  [ "$status" -eq 1 ]
}

# only run with mpi
@test "Check remap_land_parallel exists and is executable" {
  [ ! -z $skip_mpi ] && skip "not built with MPI"
  run_and_check -v remap_land_parallel

}

@test "Check make_regional_mosaic exists and is executable" {
  run_and_check -v make_regional_mosaic
  run make_regional_mosaic -h
  [ "$status" -eq 2 ]
}

@test "Check mppncscatter exists and is executable" {
  run_and_check -v mppncscatter
}

@test "Check mppnccombine exists and is executable" {
  run_and_check -v mppnccombine
  run mppnccombine -h
  [ "$status" -eq 1 ]
}

@test "Check cr_lake_files exists and is executable" {
     run_and_check -v cr_lake_files

}

@test "Check cp_river_vars exists and is executable" {
     run_and_check -v cp_river_vars
}

@test "Check rmv_parallel_rivers exists and is executable" {
     run_and_check -v rmv_parallel_rivers
}

