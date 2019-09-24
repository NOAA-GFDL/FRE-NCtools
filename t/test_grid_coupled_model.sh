#!/usr/bin/env bats

# Test grid for coupled model (land and atmosphere are C48 and ocean is 1 degree tripolar grid)

@test "Test grid for coupled model (land and atmosphere are C48 and ocean is 1 degree tripolar grid)" {

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

#Make_vgrid: create ocean_vgrid
  run command make_vgrid \
		--nbnds 3 \
		--bnds 0.,220.,5500. \
		--dbnds 10.,10.,367.14286 \
		--center c_cell \
		--grid_name ocean_vgrid 
  [ "$status" -eq 0 ]

#Make_solo_mosaic: create ocean solo mosaic
  run command make_solo_mosaic \
		--num_tiles 1 \
		--dir . \
		--mosaic_name ocean_mosaic \
		--tile_file ocean_hgrid.nc \
		--periodx 360
  [ "$status" -eq 0 ]
}
