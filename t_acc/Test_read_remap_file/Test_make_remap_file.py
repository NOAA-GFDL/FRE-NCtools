#!/usr/bin/env python3

import xarray as xr
import numpy as np
import random

nout_tiles = 2;

ngrid_in_array = [random.randint(100,500), random.randint(100,500)]
ngrid_out_array = ngrid_in_array

ncells_array  = [ ngrid_in_array[i] * ngrid_out_array[i] for i in range(0,nout_tiles) ]
ncells2_array = [ 2*ncells_array[i] for i in range(0,nout_tiles) ]

def write_answers(myfile, answer) :
  myfile.write(" ".join(answer))
  myfile.write("\n")


for out_itile in (1,nout_tiles) :

  ngrid_in, ngrid_out  = ngrid_in_array[out_itile-1], ngrid_in_array[out_itile-1]
  ncells, ncells2  = ncells_array[out_itile-1], ncells2_array[out_itile-1]

  #tile1
  tile1_data = [ random.randint(1,6) for i in range(1,ncells+1) ]
  tile1_data_attrs = dict(description="made up numbers",
                          standard_name="tile_number_in_mosaic1")
  tile1 = xr.DataArray( data = tile1_data, dims=('ncells'), attrs=tile1_data_attrs )


  #tile1_cell
  i_in = [ random.randint(1,ngrid_in) for i in range(ngrid_in) for j in range(ngrid_in) ]
  j_in = [ random.randint(1,ngrid_in) for i in range(ngrid_in) for j in range(ngrid_in) ]
  tile1_cell_data = np.array( [i_in,j_in] ).flatten('F')
  tile1_cell_attrs=dict(description="made up numbers",
                        standard_name="parent_cell_indices_in_mosaic1")
  tile1_cell = xr.DataArray( data=tile1_cell_data, dims=('ncells2'), attrs=tile1_cell_attrs )


  #tile2_cell
  i_out = [ random.randint(1,ngrid_in) for i in range(ngrid_in) for j in range(ngrid_in) ]
  j_out = [ random.randint(1,ngrid_in) for i in range(ngrid_in) for j in range(ngrid_in) ]
  tile2_cell_data = np.array( [i_out,j_out] ).flatten('F')
  tile2_cell_attrs=dict(description="made up numbers",
                        standard_name="parent_cell_indices_in_mosaic1")
  tile2_cell = xr.DataArray( data=tile2_cell_data,
                             dims=('ncells2'),
                             attrs=tile2_cell_attrs )


  #xgrid_area
  xgrid_area_data = np.array([ random.randint(1,10000)+0.1 for i in range(ncells) ], dtype=np.float64)
  xgrid_area_attrs = dict(description="made up numbers",
                          standard_name="exchange_grid_area",
                          units="m2" )
  xgrid_area = xr.DataArray( data=xgrid_area_data, dims=('ncells'), attrs=xgrid_area_attrs )


  #tile1_distance
  di_in = np.array([ random.randint(1,10000)+0.2 for i in range(ncells) ], dtype=np.float64)
  dj_in = np.array([ random.randint(1,10000)+0.3 for i in range(ncells) ], dtype=np.float64)
  tile1_distance_data = np.array( [di_in, dj_in] ).flatten('F')
  tile1_distance_attrs=dict(description="Made up numbers",
                            standard_name="distance_from_parent1_cell_centroid")
  tile1_distance = xr.DataArray( data=tile1_distance_data, dims=('ncells2'), attrs=tile1_distance_attrs )


  #dataset conserve_order1
  data_vars= {"tile1": tile1, "tile1_cell":tile1_cell, "tile2_cell":tile2_cell, "xgrid_area":xgrid_area}
  attrs=dict(description=f"""Made up test case where input and output grids cells are the same.
                              Grids are of size {ngrid_in} x {ngrid_out}, gridin with 6 tiles""")
  remap = xr.Dataset(data_vars=data_vars, attrs=attrs )
  #write netcdf file
  remap.to_netcdf(f"remap_conserve1.tile{out_itile}.nc", format="NETCDF4_CLASSIC")


  #dataset conserve_order2
  data_vars["tile1_distance"] = tile1_distance
  attrs=dict(description=f"""Made up test case where input and output grids cells are the same.
                             Grids are of size {ngrid_in} x {ngrid_out}, gridin with 6 tiles""")
  remap = xr.Dataset(data_vars=data_vars, attrs=attrs )
  #write netcdf file
  remap.to_netcdf(f"remap_conserve2.tile{out_itile}.nc", format="NETCDF4_CLASSIC")


  #write out answer file
  myfile1 = open(f"answers_conserve1.tile{out_itile}.txt", 'w')
  myfile2 = open(f"answers_conserve2.tile{out_itile}.txt", 'w')

  #answer for interp[n].nxgrid
  write_answers(myfile1,[str(ncells)])
  write_answers(myfile2,[str(ncells)])

  #answer for interp[n].interp_mini[m].nxgrid
  write_answers(myfile1, [str(tile1_data.count(i)) for i in range(1,7) ])
  write_answers(myfile2, [str(tile1_data.count(i)) for i in range(1,7) ])

  #answer for interp[n].interp_mini[m].i_in
  # -1 to be consistent with zhi....
  for itile in range(1,7) :
    answer = [ str(i_in[i]-1) for i in range(ncells) if tile1_data[i]==itile ]
    write_answers(myfile1, answer)
    write_answers(myfile2, answer)

  #answer for interp[n].interp_mini[m].j_in
  # -1 to be consistent with zhi....
  for itile in range(1,7) :
    answer = [ str(j_in[i]-1) for i in range(ncells) if tile1_data[i]==itile ]
    write_answers(myfile1, answer)
    write_answers(myfile2, answer)

  #answer for interp[n].interp_mini[m].i_out
  # -1 to be consistent with zhi....
  for itile in range(1,7) :
    answer = [ str(i_out[i]-1) for i in range(ncells) if tile1_data[i]==itile ]
    write_answers(myfile1, answer)
    write_answers(myfile2, answer)

  #answer for interp[n].interp_mini[m].j_out
  # -1 to be consistent with zhi....
  for itile in range(1,7) :
    answer = [ str(j_out[i]-1) for i in range(ncells) if tile1_data[i]==itile ]
    write_answers(myfile1, answer)
    write_answers(myfile2, answer)

  #answer for interp[n].interp_mini[m].area
  for itile in range(1,7) :
    answer = [ str(xgrid_area_data[i]) for i in range(ncells) if tile1_data[i]==itile ]
    write_answers(myfile1, answer)
    write_answers(myfile2, answer)

  #answer for interp[n].interp_mini[m].di_in
  for itile in range(1,7) :
    answer = [ str(di_in[i]) for i in range(ncells) if tile1_data[i]==itile ]
    write_answers(myfile2, answer)

  #answer for interp[n].interp_mini[m].di_in
  for itile in range(1,7) :
    answer = [ str(dj_in[i]) for i in range(ncells) if tile1_data[i]==itile ]
    write_answers(myfile2, answer)

  myfile1.close()
  myfile2.close()
