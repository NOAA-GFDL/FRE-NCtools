#!/usr/bin/env bats

# Test grid for coupled model (land and atmosphere are C48 and ocean is 1 degree tripolar grid)

@test "Test grid for coupled model (land and atmosphere are C48 and ocean is 1 degree tripolar grid)" {
cd grid_coupled_model
#Make_hgrid: create ocean_hgrid"
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
