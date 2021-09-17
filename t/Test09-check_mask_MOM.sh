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

# Test check_mask for MOM6 grid 
load test_utils

@test "Test check_mask for MOM6 grid " {

  ncgen -o OCCAM_p5degree.nc $top_srcdir/t/Test03-input/OCCAM_p5degree.ncl

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

#Make_topog: create ocean topography data
  run command make_topog \
		--mosaic ocean_mosaic.nc \
		--topog_type realistic \
		--topog_file OCCAM_p5degree.nc \
		--topog_field TOPO \
		--scale_factor -1 \
		--vgrid ocean_vgrid.nc \
		--output topog.nc 
  [ "$status" -eq 0 ]

#Run check_mask
  run command check_mask \
		--grid_file ocean_mosaic.nc \
		--ocean_topog topog.nc \
		--layout 45,72
  [ "$status" -eq 0 ]

}
