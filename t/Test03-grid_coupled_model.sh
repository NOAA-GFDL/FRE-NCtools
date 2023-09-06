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
# Test grid for coupled model (land and atmosphere are C48 and ocean is 1 degree tripolar grid)
# set setup to generate ncls from test's input directory
load test_utils

@test "Test grid for coupled model (land and atmosphere are C48 and ocean is 1 degree tripolar grid)" {

  generate_all_from_ncl

#Make_hgrid: create ocean_hgrid"
    make_hgrid \
		--grid_type tripolar_grid \
		--nxbnd 2 \
		--nybnd 7 \
		--xbnd -280,80 \
		--ybnd -82,-30,-10,0,10,30,90 \
		--dlon 1.0,1.0 \
		--dlat 1.0,1.0,0.6666667,0.3333333,0.6666667,1.0,1.0 \
		--grid_name ocean_hgrid \
		--center c_cell

#Make_vgrid: create ocean_vgrid
   make_vgrid \
		--nbnds 3 \
		--bnds 0.,220.,5500. \
		--dbnds 10.,10.,367.14286 \
		--center c_cell \
		--grid_name ocean_vgrid

#Make_solo_mosaic: create ocean solo mosaic
   make_solo_mosaic \
		--num_tiles 1 \
		--dir . \
		--mosaic_name ocean_mosaic \
		--tile_file ocean_hgrid.nc \
		--periodx 360

#Make_topog: create ocean topography data
   make_topog \
		--mosaic ocean_mosaic.nc \
		--topog_type realistic \
		--topog_file OCCAM_p5degree.nc \
		--topog_field TOPO \
		--scale_factor -1 \
		--vgrid ocean_vgrid.nc \
		--output topog.nc

# MPI only tests
  if [ -z "$skip_mpi" ] ; then
       mpirun -n 2 make_topog_parallel \
		--mosaic ocean_mosaic.nc \
		--topog_type realistic \
		--topog_file OCCAM_p5degree.nc \
		--topog_field TOPO \
		--scale_factor -1 \
		--vgrid ocean_vgrid.nc \
		--output topog_parallel.nc

    run_and_check  nccmp -md topog.nc topog_parallel.nc
  fi


#Make_hgrid: create C48 grid for atmos/land
   make_hgrid \
		--grid_type gnomonic_ed \
		--nlon 96 \
		--grid_name C48_grid \

#Make_solo_mosaic: create C48 solo mosaic for atmos/land
   make_solo_mosaic \
		--num_tiles 6 \
		--dir . \
		--mosaic C48_mosaic \
		--tile_file C48_grid.tile1.nc,C48_grid.tile2.nc,C48_grid.tile3.nc,C48_grid.tile4.nc,C48_grid.tile5.nc,C48_grid.tile6.nc

#Make_coupler_mosaic: coupler_mosaic with and without parallel and compare
   make_coupler_mosaic \
		--atmos_mosaic C48_mosaic.nc \
		--ocean_mosaic ocean_mosaic.nc \
		--ocean_topog  topog.nc \
		--check \
		--area_ratio_thresh 1.e-10 \
		--mosaic_name grid_spec

#Make coupler mosaic with parallel
  if [ -z "$skip_mpi" ] ; then
      [ ! -d parallel ] && mkdir parallel
      cd parallel

       mpirun -n 4 make_coupler_mosaic_parallel \
		--atmos_mosaic ../C48_mosaic.nc \
		--ocean_mosaic ../ocean_mosaic.nc \
		--ocean_topog  ../topog.nc \
		--mosaic_name grid_spec  \
		--area_ratio_thresh 1.e-10

      # compare any created files to non-parallel (exclude directory differences)
      run_and_check nccmp -md --exclude=atm_mosaic_dir --exclude=lnd_mosaic_dir --exclude=ocn_mosaic_dir \
                --exclude=ocn_topog_dir grid_spec.nc ../grid_spec.nc
  fi
}
