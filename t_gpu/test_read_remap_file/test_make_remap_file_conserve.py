#!/usr/bin/env python3

#***********************************************************************
#                   GNU Lesser General Public License
#
# This file is part of the GFDL FRE NetCDF tools package (FRE-NCTools).
#
# FRE-NCtools is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at
# your option) any loner version.
#
# FRE-NCtools is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
# for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with FRE-NCTools.  If not, see
# <http://www.gnu.org/licenses/>.
# **********************************************************************

# This python script generates made-up remap files for first order
# and second order interpolation scheme where the number
# output grid tile2 = 2.  In addition, this script generate answer
# files in ASCII format.

import xarray as xr
import numpy as np
import random
import sys

output_grid_ntiles = 2;

input_grid_ncells_array = [random.randint(100,500), random.randint(100,500)]
output_grid_ncells_array = input_grid_ncells_array

nxcells_per_output_tile  = [ input_grid_ncells_array[i] * output_grid_ncells_array[i]
                             for i in range(0,output_grid_ntiles) ]
nxcells2_per_output_tile = [ 2*nxcells_per_output_tile[i] for i in range(0,output_grid_ntiles) ]

interp_method = int(sys.argv[1]) #1 or 2

def write_answers(myfile, answer) :
  myfile.write(" ".join(answer))
  myfile.write("\n")


for otile in (1,output_grid_ntiles) :

  input_grid_ncells, output_grid_ncells  = input_grid_ncells_array[otile-1], output_grid_ncells_array[otile-1]
  nxcells, nxcells2  = nxcells_per_output_tile[otile-1], nxcells2_per_output_tile[otile-1]

  #tile1
  tile1_data = [ random.randint(1,6) for i in range(1,nxcells+1) ]
  tile1_data_attrs = dict(description="made up numbers",
                          standard_name="tile_number_in_mosaic1")
  tile1 = xr.DataArray( data = tile1_data, dims=('ncells'), attrs=tile1_data_attrs )


  #tile1_cell
  i_in = [ random.randint(1,input_grid_ncells) for i in range(nxcells) ]
  j_in = [ random.randint(1,input_grid_ncells) for i in range(nxcells) ]
  tile1_cell_data = np.array( [i_in,j_in] ).flatten('F')
  tile1_cell_attrs=dict(description="made up numbers",
                        standard_name="parent_cell_indices_in_mosaic1")
  tile1_cell = xr.DataArray( data=tile1_cell_data, dims=('ncells2'), attrs=tile1_cell_attrs )


  #tile2_cell
  i_out = [ random.randint(1,output_grid_ncells) for i in range(nxcells) ]
  j_out = [ random.randint(1,output_grid_ncells) for i in range(nxcells) ]
  tile2_cell_data = np.array( [i_out,j_out] ).flatten('F')
  tile2_cell_attrs=dict(description="made up numbers",
                        standard_name="parent_cell_indices_in_mosaic1")
  tile2_cell = xr.DataArray( data=tile2_cell_data,
                             dims=('ncells2'),
                             attrs=tile2_cell_attrs )


  #xgrid_area
  xgrid_area_data = np.array([ random.randint(1,10000)+0.1 for i in range(nxcells) ], dtype=np.float64)
  xgrid_area_attrs = dict(description="made up numbers",
                          standard_name="exchange_grid_area",
                          units="m2" )
  xgrid_area = xr.DataArray( data=xgrid_area_data, dims=('ncells'), attrs=xgrid_area_attrs )


  if( interp_method == 2 ) :
    #tile1_distance
    di_in = np.array([ random.randint(1,10000)+0.2 for i in range(nxcells) ], dtype=np.float64)
    dj_in = np.array([ random.randint(1,10000)+0.3 for i in range(nxcells) ], dtype=np.float64)
    tile1_distance_data = np.array( [di_in, dj_in] ).flatten('F')
    tile1_distance_attrs=dict(description="Made up numbers",
                              standard_name="distance_from_parent1_cell_centroid")
    tile1_distance = xr.DataArray( data=tile1_distance_data, dims=('ncells2'), attrs=tile1_distance_attrs )


  #dataset conserve_order1
  data_vars = {"tile1": tile1, "tile1_cell":tile1_cell, "tile2_cell":tile2_cell, "xgrid_area":xgrid_area}
  if( interp_method == 2 ) : data_vars["tile1_distance"] = tile1_distance
  remap = xr.Dataset(data_vars=data_vars)
  #write netcdf file
  remap.to_netcdf(f"remap_conserve{interp_method}.tile{otile}.nc", format="NETCDF4_CLASSIC")

  #write out answer file
  myfile = open(f"answers_conserve{interp_method}.tile{otile}.txt", 'w')

  #answer for nxcells
  write_answers(myfile,[str(nxcells)]) #total_nxcells
  write_answers(myfile, [str(tile1_data.count(i)) for i in range(1,7) ]) #nxcells per gridin tile

  #answer for interp[n].per_intile[m].i_in
  # -1 to be consistent with zhi....
  for itile in range(1,7) :
    answer = [ str(i_in[i]-1) for i in range(nxcells) if tile1_data[i]==itile ]
    write_answers(myfile, answer)

  #answer for interp[n].per_intile[m].j_in
  # -1 to be consistent with zhi....
  for itile in range(1,7) :
    answer = [ str(j_in[i]-1) for i in range(nxcells) if tile1_data[i]==itile ]
    write_answers(myfile, answer)

  #answer for interp[n].per_intile[m].i_out
  # -1 to be consistent with zhi....
  for itile in range(1,7) :
    answer = [ str(i_out[i]-1) for i in range(nxcells) if tile1_data[i]==itile ]
    write_answers(myfile, answer)

  #answer for interp[n].per_intile[m].j_out
  # -1 to be consistent with zhi....
  for itile in range(1,7) :
    answer = [ str(j_out[i]-1) for i in range(nxcells) if tile1_data[i]==itile ]
    write_answers(myfile, answer)

  #answer for interp[n].per_intile[m].area
  for itile in range(1,7) :
    answer = [ str(xgrid_area_data[i]) for i in range(nxcells) if tile1_data[i]==itile ]
    write_answers(myfile, answer)

  if( interp_method == 2 ) :
    #answer for interp[n].per_intile[m].di_in
    for itile in range(1,7) :
      answer = [ str(di_in[i]) for i in range(nxcells) if tile1_data[i]==itile ]
      write_answers(myfile, answer)
    #answer for interp[n].per_intile[m].di_in
    for itile in range(1,7) :
      answer = [ str(dj_in[i]) for i in range(nxcells) if tile1_data[i]==itile ]
      write_answers(myfile, answer)

  myfile.close()
