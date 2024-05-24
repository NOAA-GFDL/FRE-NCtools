/***********************************************************************
 *                   GNU Lesser General Public License
 *
 * This file is part of the GFDL FRE NetCDF tools package (FRE-NCTools).
 *
 * FRE-NCtools is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at
 * your option) any later version.
 *
 * FRE-NCtools is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with FRE-NCTools.  If not, see
 * <http://www.gnu.org/licenses/>.
 **********************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <netcdf.h>
#include <math.h>
#include "globals_acc.h"
#include "conserve_interp_acc.h"
#include "interp_utils_acc.h"
#include "create_xgrid_acc.h"
#include "create_xgrid_utils_acc.h"
#include "general_utils_acc.h"
#include "mpp.h"
#include "mpp_io.h"
#include "read_mosaic.h"

/*******************************************************************************
  void setup_conserve_interp
  Setup the interpolation weight for conservative interpolation
*******************************************************************************/
void setup_conserve_interp_acc(int ntiles_input_grid, const Grid_config *input_grid, int ntiles_output_grid,
			   Grid_config *output_grid, Xgrid_config *xgrid, unsigned int opcode)
{
  int n, m, i, ii, jj, nlon_input_cells, nlat_input_cells;
  int nlon_output_cells, nlat_output_cells, tile;
  size_t nxcells, nxcells_prev;
  int *input_parent_lon_indices=NULL, *input_parent_lat_indices=NULL;
  int *output_parent_lon_indices=NULL, *output_parent_lat_indices=NULL;
  int *tmp_input_parent_tile=NULL, *tmp_input_parent_lon_indices=NULL;
  int *tmp_input_parent_lat_indices=NULL, *tmp_output_parent_lon_indices=NULL, *tmp_output_parent_lat_indices=NULL;
  double *tmp_dcentroid_lon, *tmp_dcentroid_lat;
  double *xcell_area=NULL, *tmp_area=NULL, *xcell_centroid_lon=NULL, *xcell_centroid_lat=NULL;

  Grid_cells_struct_config *input_cells;

  if( opcode & READ) {
    read_remap_file_acc(ntiles_input_grid, ntiles_output_grid, xgrid, opcode);
    copy_xgrid_to_device_acc(ntiles_input_grid, ntiles_output_grid, xgrid, opcode);
  }
  else {
    input_parent_lon_indices       = (int    *)malloc(MAXXGRID   * sizeof(int   ));
    input_parent_lat_indices       = (int    *)malloc(MAXXGRID   * sizeof(int   ));
    output_parent_lon_indices      = (int    *)malloc(MAXXGRID   * sizeof(int   ));
    output_parent_lat_indices      = (int    *)malloc(MAXXGRID   * sizeof(int   ));
    xcell_area = (double *)malloc(MAXXGRID   * sizeof(double));
    xcell_centroid_lon = (double *)malloc(MAXXGRID   * sizeof(double));
    xcell_centroid_lat = (double *)malloc(MAXXGRID   * sizeof(double));
    input_cells    = (Grid_cells_struct_config *)malloc(ntiles_input_grid * sizeof(Grid_cells_struct_config));
    for(m=0; m<ntiles_input_grid; m++) {
      nlon_input_cells = input_grid[m].nx;
      nlat_input_cells = input_grid[m].ny;
      input_cells[m].recomputed_area = (double *)malloc(nlon_input_cells*nlat_input_cells*sizeof(double));
      input_cells[m].centroid_lon = (double *)malloc(nlon_input_cells*nlat_input_cells*sizeof(double));
      input_cells[m].centroid_lat = (double *)malloc(nlon_input_cells*nlat_input_cells*sizeof(double));
      for(n=0; n<nlon_input_cells*nlat_input_cells; n++) {
        input_cells[m].recomputed_area[n] = 0;
        input_cells[m].centroid_lon[n] = 0;
        input_cells[m].centroid_lat[n] = 0;
      }
    }
    for(n=0; n<ntiles_output_grid; n++) {

      nlon_output_cells = output_grid[n].nxc;
      nlat_output_cells = output_grid[n].nyc;
      xgrid[n].nxcells = 0;
      for(m=0; m<ntiles_input_grid; m++){
        int jlat_overlap_starts_offset, nlat_overlapping_cells;
        int start_from_this_corner;

        nlon_input_cells = input_grid[m].nx;
        nlat_input_cells = input_grid[m].ny;

        get_skip_cells(nlon_input_cells*nlat_input_cells, &(input_cells[m].skip_cells));

        if(opcode & GREAT_CIRCLE) {
          nxcells = create_xgrid_great_circle_acc(&nlon_input_cells, &nlat_input_cells,
                                                  &nlon_output_cells, &nlat_output_cells,
                                                  input_grid[m].lonc, input_grid[m].latc,
                                                  output_grid[n].lonc, output_grid[n].latc, input_cells[m].skip_cells,
                                                  input_parent_lon_indices, input_parent_lat_indices,
                                                  output_parent_lon_indices, output_parent_lat_indices,
                                                  xcell_area, xcell_centroid_lon, xcell_centroid_lat);
        }
        else {
          //get the input grid portion (bounding indices) that overlaps with the output grid in the latitudonal direction.
          get_bounding_indices(nlon_output_cells, nlat_output_cells, nlon_input_cells, nlat_input_cells,
                               output_grid[n].latc, input_grid[m].latc, &jlat_overlap_starts_offset, &nlat_overlapping_cells);
          start_from_this_corner = jlat_overlap_starts_offset*(nlon_input_cells+1);

          if(opcode & CONSERVE_ORDER1) {
            nxcells = create_xgrid_2dx2d_order1_acc(&nlon_input_cells, &nlat_overlapping_cells,
                                                    &nlon_output_cells, &nlat_output_cells,
                                                    input_grid[m].lonc+start_from_this_corner,
                                                    input_grid[m].latc+start_from_this_corner,
                                                    output_grid[n].lonc, output_grid[n].latc, input_cells[m].skip_cells,
                                                    input_parent_lon_indices, input_parent_lat_indices,
                                                    output_parent_lon_indices, output_parent_lat_indices, xcell_area);
            for(i=0; i<nxcells; i++) input_parent_lat_indices[i] += jlat_overlap_starts_offset;
          }
          else if(opcode & CONSERVE_ORDER2) {
            int g_nxcells;
            int    *g_input_parent_lon_indices, *g_input_parent_lat_indices;
            double *g_xcell_area, *g_xcell_centroid_lon, *g_xcell_centroid_lat;

            nxcells = create_xgrid_2dx2d_order2_acc(&nlon_input_cells, &nlat_overlapping_cells,
                                                    &nlon_output_cells, &nlat_output_cells,
                                                    input_grid[m].lonc+start_from_this_corner,
                                                    input_grid[m].latc+start_from_this_corner,
                                                    output_grid[n].lonc,  output_grid[n].latc, input_cells[m].skip_cells,
                                                    input_parent_lon_indices, input_parent_lat_indices,
                                                    output_parent_lon_indices, output_parent_lat_indices, xcell_area,
                                                    xcell_centroid_lon, xcell_centroid_lat);
            for(i=0; i<nxcells; i++) input_parent_lat_indices[i] += jlat_overlap_starts_offset;

	    /* For the purpose of bitiwise reproducing, the following operation is needed. */
	    g_nxcells = nxcells;
	    mpp_sum_int(1, &g_nxcells);
	    if(g_nxcells > 0) {
	      g_input_parent_lon_indices = (int *)malloc(g_nxcells*sizeof(int));
	      g_input_parent_lat_indices = (int *)malloc(g_nxcells*sizeof(int));
	      g_xcell_area = (double *)malloc(g_nxcells*sizeof(double));
	      g_xcell_centroid_lon = (double *)malloc(g_nxcells*sizeof(double));
	      g_xcell_centroid_lat = (double *)malloc(g_nxcells*sizeof(double));
	      mpp_gather_field_int   (nxcells, input_parent_lon_indices,  g_input_parent_lon_indices);
	      mpp_gather_field_int   (nxcells, input_parent_lat_indices,  g_input_parent_lat_indices);
	      mpp_gather_field_double(nxcells, xcell_area, g_xcell_area);
	      mpp_gather_field_double(nxcells, xcell_centroid_lon, g_xcell_centroid_lon);
	      mpp_gather_field_double(nxcells, xcell_centroid_lat, g_xcell_centroid_lat);
	      for(i=0; i<g_nxcells; i++) {
          ii = g_input_parent_lat_indices[i]*nlon_input_cells+g_input_parent_lon_indices[i];
          input_cells[m].recomputed_area[ii] += g_xcell_area[i];
          input_cells[m].centroid_lon[ii] += g_xcell_centroid_lon[i];
          input_cells[m].centroid_lat[ii] += g_xcell_centroid_lat[i];
	      }
	      free(g_input_parent_lon_indices);
	      free(g_input_parent_lat_indices);
	      free(g_xcell_area);
	      free(g_xcell_centroid_lon);
	      free(g_xcell_centroid_lat);
	    }
	  }
	  else
	    mpp_error("conserve_interp: interp_method should be CONSERVE_ORDER1 or CONSERVE_ORDER2");
	}

	if(nxcells > 0) {
	  nxcells_prev = xgrid[n].nxcells;
	  xgrid[n].nxcells += nxcells;
	  if(nxcells_prev == 0 ) {
	    xgrid[n].input_parent_lon_indices   = (int    *)malloc(xgrid[n].nxcells*sizeof(int   ));
	    xgrid[n].input_parent_lat_indices   = (int    *)malloc(xgrid[n].nxcells*sizeof(int   ));
	    xgrid[n].output_parent_lon_indices  = (int    *)malloc(xgrid[n].nxcells*sizeof(int   ));
	    xgrid[n].output_parent_lat_indices  = (int    *)malloc(xgrid[n].nxcells*sizeof(int   ));
	    xgrid[n].xcell_area        = (double *)malloc(xgrid[n].nxcells*sizeof(double));
	    xgrid[n].input_parent_tile   = (int    *)malloc(xgrid[n].nxcells*sizeof(int   ));
	    for(i=0; i<xgrid[n].nxcells; i++) {
	      xgrid[n].input_parent_tile [i] = m;
	      xgrid[n].input_parent_lon_indices [i] = input_parent_lon_indices [i];
	      xgrid[n].input_parent_lat_indices [i] = input_parent_lat_indices [i];
	      xgrid[n].output_parent_lon_indices[i] = output_parent_lon_indices[i];
	      xgrid[n].output_parent_lat_indices[i] = output_parent_lat_indices[i];
	      xgrid[n].xcell_area[i]  = xcell_area[i];
	    }
	    if(opcode & CONSERVE_ORDER2) {
	      xgrid[n].dcentroid_lon   = (double *)malloc(xgrid[n].nxcells*sizeof(double));
	      xgrid[n].dcentroid_lat   = (double *)malloc(xgrid[n].nxcells*sizeof(double));
	      for(i=0; i<xgrid[n].nxcells; i++) {
		xgrid[n].dcentroid_lon [i] = xcell_centroid_lon[i]/xcell_area[i];
		xgrid[n].dcentroid_lat [i] = xcell_centroid_lat[i]/xcell_area[i];
	      }
	    }
	  }
	  else {
	    tmp_input_parent_lon_indices  = xgrid[n].input_parent_lon_indices;
	    tmp_input_parent_lat_indices  = xgrid[n].input_parent_lat_indices;
	    tmp_output_parent_lon_indices = xgrid[n].output_parent_lon_indices;
	    tmp_output_parent_lat_indices = xgrid[n].output_parent_lat_indices;
	    tmp_area  = xgrid[n].xcell_area;
	    tmp_input_parent_tile  = xgrid[n].input_parent_tile;
	    xgrid[n].input_parent_lon_indices   = (int    *)malloc(xgrid[n].nxcells*sizeof(int   ));
	    xgrid[n].input_parent_lat_indices   = (int    *)malloc(xgrid[n].nxcells*sizeof(int   ));
	    xgrid[n].output_parent_lon_indices  = (int    *)malloc(xgrid[n].nxcells*sizeof(int   ));
	    xgrid[n].output_parent_lat_indices  = (int    *)malloc(xgrid[n].nxcells*sizeof(int   ));
	    xgrid[n].xcell_area   = (double *)malloc(xgrid[n].nxcells*sizeof(double));
	    xgrid[n].input_parent_tile   = (int    *)malloc(xgrid[n].nxcells*sizeof(int   ));
	    for(i=0; i<nxcells_prev; i++) {
	      xgrid[n].input_parent_tile [i] = tmp_input_parent_tile [i];
	      xgrid[n].input_parent_lon_indices [i] = tmp_input_parent_lon_indices [i];
	      xgrid[n].input_parent_lat_indices [i] = tmp_input_parent_lat_indices [i];
	      xgrid[n].output_parent_lon_indices[i] = tmp_output_parent_lon_indices[i];
	      xgrid[n].output_parent_lat_indices[i] = tmp_output_parent_lat_indices[i];
	      xgrid[n].xcell_area [i] = tmp_area [i];
	    }
	    for(i=0; i<nxcells; i++) {
	      ii = i + nxcells_prev;
	      xgrid[n].input_parent_tile [ii] = m;
	      xgrid[n].input_parent_lon_indices [ii] = input_parent_lon_indices [i];
	      xgrid[n].input_parent_lat_indices [ii] = input_parent_lat_indices [i];
	      xgrid[n].output_parent_lon_indices[ii] = output_parent_lon_indices[i];
	      xgrid[n].output_parent_lat_indices[ii] = output_parent_lat_indices[i];
	      xgrid[n].xcell_area [ii] = xcell_area[i];
	    }
	    if(opcode & CONSERVE_ORDER2) {
	      tmp_dcentroid_lon  = xgrid[n].dcentroid_lon;
	      tmp_dcentroid_lat  = xgrid[n].dcentroid_lat;
	      xgrid[n].dcentroid_lon   = (double *)malloc(xgrid[n].nxcells*sizeof(double));
	      xgrid[n].dcentroid_lat   = (double *)malloc(xgrid[n].nxcells*sizeof(double));
	      for(i=0; i<nxcells_prev; i++) {
		xgrid[n].dcentroid_lon [i] = tmp_dcentroid_lon [i];
		xgrid[n].dcentroid_lat [i] = tmp_dcentroid_lat [i];
	      }
	      for(i=0; i<nxcells; i++) {
		ii = i + nxcells_prev;
		jj = input_parent_lat_indices [i]*nlon_input_cells+input_parent_lon_indices [i];
		xgrid[n].dcentroid_lon [ii] = xcell_centroid_lon[i]/xcell_area[i];
		xgrid[n].dcentroid_lat [ii] = xcell_centroid_lat[i]/xcell_area[i];
	      }
	      free(tmp_dcentroid_lon);
	      free(tmp_dcentroid_lat);
	    }
	    free(tmp_input_parent_tile);
	    free(tmp_input_parent_lon_indices);
	    free(tmp_input_parent_lat_indices);
	    free(tmp_output_parent_lon_indices);
	    free(tmp_output_parent_lat_indices);
	    free(tmp_area);
	  }
	}  /* if(nxcells>0) */
      }
    }
    if(opcode & CONSERVE_ORDER2) {
      /* subtract the input_grid centroid_lon and xgrid_centroid_lat to get the distance between xgrid and input_grid */
      for(n=0; n<ntiles_input_grid; n++) {
        double x1_in[50], y1_in[50], lon_in_avg, centroid_lon, centroid_lat;
        int    j, n0, n1, n2, n3, n1_in;
        /* calcualte cell area */
        nlon_input_cells = input_grid[n].nx;
        nlat_input_cells = input_grid[n].ny;
        for(j=0; j<nlat_input_cells; j++) for(i=0; i<nlon_input_cells; i++) {
            ii = j*nlon_input_cells + i;
            if(input_cells[n].recomputed_area[ii] > 0) {
              if( fabs(input_cells[n].recomputed_area[ii]-input_grid[n].cell_area[ii])/input_grid[n].cell_area[ii] < AREA_RATIO ) {
                input_cells[n].centroid_lon[ii] /= input_cells[n].recomputed_area[ii];
                input_cells[n].centroid_lat[ii] /= input_cells[n].recomputed_area[ii];
              }
              else {
                n0 = j*(nlon_input_cells+1)+i;       n1 = j*(nlon_input_cells+1)+i+1;
                n2 = (j+1)*(nlon_input_cells+1)+i+1; n3 = (j+1)*(nlon_input_cells+1)+i;
                x1_in[0] = input_grid[n].lonc[n0]; y1_in[0] = input_grid[n].latc[n0];
                x1_in[1] = input_grid[n].lonc[n1]; y1_in[1] = input_grid[n].latc[n1];
                x1_in[2] = input_grid[n].lonc[n2]; y1_in[2] = input_grid[n].latc[n2];
                x1_in[3] = input_grid[n].lonc[n3]; y1_in[3] = input_grid[n].latc[n3];
                n1_in = fix_lon_acc(x1_in, y1_in, 4, M_PI);
                lon_in_avg = avgval_double_acc(n1_in, x1_in);
                centroid_lon = poly_ctrlon_acc(x1_in, y1_in, n1_in, lon_in_avg);
                centroid_lat = poly_ctrlat_acc(x1_in, y1_in, n1_in );
                input_cells[n].centroid_lon[ii] = centroid_lon/input_grid[n].cell_area[ii];
                input_cells[n].centroid_lat[ii] = centroid_lat/input_grid[n].cell_area[ii];
              }
            }
          }
      }
      for(n=0; n<ntiles_output_grid; n++) {
        for(i=0; i<xgrid[n].nxcells; i++) {
          tile = xgrid[n].input_parent_tile[i];
          ii   = xgrid[n].input_parent_lat_indices[i] * input_grid[tile].nx + xgrid[n].input_parent_lon_indices[i];
          xgrid[n].dcentroid_lon[i] -= input_cells[tile].centroid_lon[ii];
          xgrid[n].dcentroid_lat[i] -= input_cells[tile].centroid_lat[ii];
        }
      }

      /* free the memory */
      for(n=0; n<ntiles_input_grid; n++) {
        free(input_cells[n].recomputed_area);
        free(input_cells[n].centroid_lon);
        free(input_cells[n].centroid_lat);
      }
      free(input_cells);
    }
    if( opcode & WRITE) { /* write out remapping information */
      for(n=0; n<ntiles_output_grid; n++) {
	int nxcells;

	nxcells = xgrid[n].nxcells;
	mpp_sum_int(1, &nxcells);
	if(nxcells > 0) {
	  size_t start[4], nwrite[4];
	  int    fid, dim_string, dim_ncells, dim_two, dims[4];
	  int    id_xcell_area, id_tile1_dist;
	  int    id_tile1_cell, id_tile2_cell, id_tile1;
	  int    *gdata_int, *ldata_int;
	  double *gdata_dbl;

	  fid = mpp_open( xgrid[n].remap_file, MPP_WRITE);
	  dim_string = mpp_def_dim(fid, "string", STRING);
	  dim_ncells = mpp_def_dim(fid, "ncells", nxcells);
	  dim_two    = mpp_def_dim(fid, "two", 2);
	  dims[0] = dim_ncells; dims[1] = dim_two;
	  id_tile1      = mpp_def_var(fid, "tile1",      NC_INT, 1, &dim_ncells, 1,
				      "standard_name", "tile_number_in_mosaic1");
	  id_tile1_cell = mpp_def_var(fid, "tile1_cell", NC_INT, 2, dims, 1,
				      "standard_name", "parent_cell_indices_in_mosaic1");
	  id_tile2_cell = mpp_def_var(fid, "tile2_cell", NC_INT, 2, dims, 1,
				      "standard_name", "parent_cell_indices_in_mosaic2");
	  id_xcell_area = mpp_def_var(fid, "xcell_area", NC_DOUBLE, 1, &dim_ncells, 2,
				      "standard_name", "exchange_grid_area", "units", "m2");
	  if(opcode & CONSERVE_ORDER2) id_tile1_dist = mpp_def_var(fid, "tile1_distance", NC_DOUBLE, 2, dims, 1,
								   "standard_name", "distance_from_parent1_cell_centroid");
	  mpp_end_def(fid);
	  for(i=0; i<4; i++) {
	    start[i] = 0; nwrite[i] = 1;
	  }
	  nwrite[0] = nxcells;
	  gdata_int = (int *)malloc(nxcells*sizeof(int));
	  if(xgrid[n].nxcells>0) ldata_int = (int *)malloc(xgrid[n].nxcells*sizeof(int));
          mpp_gather_field_int(xgrid[n].nxcells, xgrid[n].input_parent_tile, gdata_int);
	  for(i=0; i<nxcells; i++) gdata_int[i]++;
	  mpp_put_var_value(fid, id_tile1, gdata_int);

	  mpp_gather_field_int(xgrid[n].nxcells, xgrid[n].input_parent_lon_indices, gdata_int);
	  for(i=0; i<nxcells; i++) gdata_int[i]++;
	  mpp_put_var_value_block(fid, id_tile1_cell, start, nwrite, gdata_int);

          for(i=0; i<xgrid[n].nxcells; i++) ldata_int[i] = xgrid[n].output_parent_lon_indices[i] + output_grid[n].isc + 1;
	  mpp_gather_field_int(xgrid[n].nxcells, ldata_int, gdata_int);
	  mpp_put_var_value_block(fid, id_tile2_cell, start, nwrite, gdata_int);

	  mpp_gather_field_int(xgrid[n].nxcells, xgrid[n].input_parent_lat_indices, gdata_int);
	  for(i=0; i<nxcells; i++) gdata_int[i]++;
	  start[1] = 1;
	  mpp_put_var_value_block(fid, id_tile1_cell, start, nwrite, gdata_int);

          for(i=0; i<xgrid[n].nxcells; i++) ldata_int[i] = xgrid[n].output_parent_lat_indices[i] + output_grid[n].jsc + 1;
	  mpp_gather_field_int(xgrid[n].nxcells, ldata_int, gdata_int);
	  mpp_put_var_value_block(fid, id_tile2_cell, start, nwrite, gdata_int);

	  free(gdata_int);
	  if(xgrid[n].nxcells>0)free(ldata_int);

	  gdata_dbl = (double *)malloc(nxcells*sizeof(double));
	  mpp_gather_field_double(xgrid[n].nxcells, xgrid[n].xcell_area, gdata_dbl);
	  mpp_put_var_value(fid, id_xcell_area, gdata_dbl);

	  if(opcode & CONSERVE_ORDER2) {
	    start[1] = 0;
	    mpp_gather_field_double(xgrid[n].nxcells, xgrid[n].dcentroid_lon, gdata_dbl);
	    mpp_put_var_value_block(fid, id_tile1_dist, start, nwrite, gdata_dbl);
	    start[1] = 1;
	    mpp_gather_field_double(xgrid[n].nxcells, xgrid[n].dcentroid_lat, gdata_dbl);
	    mpp_put_var_value_block(fid, id_tile1_dist, start, nwrite, gdata_dbl);
	  }

	  free(gdata_dbl);
	  mpp_close(fid);
	}
      }
    }
    if(mpp_pe() == mpp_root_pe())printf("NOTE: done calculating index and weight for conservative interpolation\n");
  }

  /* check the input area match exchange grid area */
  if(opcode & CHECK_CONSERVE) {
    int nx1, ny1, max_i, max_j, i, j;
    double max_ratio, ratio_change;
    double *recomputed_output_area;

    /* sum over exchange grid to get the area of input_grid */
    nx1  = output_grid[0].nxc;
    ny1  = output_grid[0].nyc;

    recomputed_output_area = (double *)malloc(nx1*ny1*sizeof(double));

    for(n=0; n<ntiles_output_grid; n++) {
      for(i=0; i<nx1*ny1; i++) recomputed_output_area[i] = 0;
      for(i=0; i<xgrid[n].nxcells; i++) {
	ii = xgrid[n].output_parent_lat_indices[i]*nx1 + xgrid[n].output_parent_lon_indices[i];
	recomputed_output_area[ii] +=  xgrid[n].xcell_area[i];
      }
      max_ratio = 0;
      max_i = 0;
      max_j = 0;
      /* comparing area1 and recomputed_output_area */
      for(j=0; j<ny1; j++) for(i=0; i<nx1; i++) {
	ii = j*nx1+i;
	ratio_change = fabs(output_grid[n].cell_area[ii]-recomputed_output_area[ii])/output_grid[n].cell_area[ii];
	if(ratio_change > max_ratio) {
	  max_ratio = ratio_change;
	  max_i = i;
	  max_j = j;
	}
	if( ratio_change > 1.e-4 ) {
	  printf("(i,j)=(%d,%d), change = %g, area1=%g, recomputed_output_area=%g\n",
           i, j, ratio_change, output_grid[n].cell_area[ii],recomputed_output_area[ii]);
	}
      }
      ii = max_j*nx1+max_i;
      printf("The maximum ratio change at (%d,%d) = %g, area1=%g, recomputed_output_area=%g\n",
             max_i, max_j, max_ratio, output_grid[n].cell_area[ii],recomputed_output_area[ii]);

    }

    free(recomputed_output_area);

  }

  free(input_parent_lon_indices);
  free(input_parent_lat_indices);
  free(output_parent_lon_indices);
  free(output_parent_lat_indices);
  free(xcell_area);
  if(xcell_centroid_lon) free(xcell_centroid_lon);
  if(xcell_centroid_lat) free(xcell_centroid_lat);

}; /* setup_conserve_interp */


/*******************************************************************************
 void read_remap_file
 Reads in the weight/remap file if provided and copies the data to the device
*******************************************************************************/
void read_remap_file_acc(int ntiles_input_grid, int ntiles_output_grid, Xgrid_config *xgrid, unsigned int opcode)
{

  int *input_parent_lon_indices=NULL, *input_parent_lat_indices=NULL, *output_parent_lon_indices=NULL, *output_parent_lat_indices=NULL;
  double *xcell_area=NULL, *xcell_centroid_lon=NULL, *xcell_centroid_lat=NULL;

  size_t nxcells;
  int nxcells_acc;

  for(int n=0; n<ntiles_output_grid; n++) {
    if( xgrid[n].file_exist ) {
      int *t_in, *ind_acc;
      int fid, vid;

      nxcells = read_mosaic_xgrid_size(xgrid[n].remap_file);
      xgrid[n].nxcells = nxcells;

      t_in  = (int *)malloc(nxcells*sizeof(int));
      input_parent_lon_indices  = (int *)malloc(nxcells*sizeof(int));
      input_parent_lat_indices  = (int *)malloc(nxcells*sizeof(int));
      output_parent_lon_indices = (int *)malloc(nxcells*sizeof(int));
      output_parent_lat_indices = (int *)malloc(nxcells*sizeof(int));
      xcell_area = (double *)malloc(nxcells*sizeof(double));
      if( opcode & CONSERVE_ORDER2) {
        xcell_centroid_lon = (double *)malloc(nxcells*sizeof(double));
        xcell_centroid_lat = (double *)malloc(nxcells*sizeof(double));
      }

      if(opcode & CONSERVE_ORDER1)
        read_mosaic_xgrid_order1(xgrid[n].remap_file, input_parent_lon_indices, input_parent_lat_indices,
                                 output_parent_lon_indices, output_parent_lat_indices, xcell_area);
      else
        read_mosaic_xgrid_order2(xgrid[n].remap_file, input_parent_lon_indices, input_parent_lat_indices,
                                 output_parent_lon_indices, output_parent_lat_indices, xcell_area,
                                 xcell_centroid_lon, xcell_centroid_lat);

      //rescale the xgrid area
      for(int i=0; i<xgrid[n].nxcells; i++) xcell_area[i] *= GAREA;

      //read in tile number of input parent cells
      fid = mpp_open(xgrid[n].remap_file, MPP_READ);
      vid = mpp_get_varid(fid, "tile1");
      mpp_get_var_value(fid, vid, t_in);
      mpp_close(fid);

      //tile number starts from 1, not 0, in the weight files
      for(int i=0 ; i<nxcells ; i++) t_in[i]--;

      //get number of nxcells per input tile
      for(int m=0 ; m<ntiles_input_grid ; m++) xgrid[n].per_intile[m].nxcells=0;
      for(int i=0 ; i<nxcells ; i++) xgrid[n].per_intile[ t_in[i] ].nxcells++;

      for(int m=0 ; m<ntiles_input_grid ; m++) {
        nxcells_acc = xgrid[n].per_intile[m].nxcells;
        xgrid[n].per_intile[m].input_parent_lon_indices  = (int *)malloc(nxcells_acc*sizeof(int));
        xgrid[n].per_intile[m].input_parent_lat_indices  = (int *)malloc(nxcells_acc*sizeof(int));
        xgrid[n].per_intile[m].output_parent_lon_indices = (int *)malloc(nxcells_acc*sizeof(int));
        xgrid[n].per_intile[m].output_parent_lat_indices = (int *)malloc(nxcells_acc*sizeof(int));
        xgrid[n].per_intile[m].xcell_area  = (double *)malloc(nxcells_acc*sizeof(double));
        if(opcode & CONSERVE_ORDER2) {
          xgrid[n].per_intile[m].dcentroid_lon = (double *)malloc(nxcells_acc*sizeof(double));
          xgrid[n].per_intile[m].dcentroid_lat = (double *)malloc(nxcells_acc*sizeof(double));
        }
      }

      ind_acc = (int *)calloc(ntiles_input_grid, sizeof(int));
      for(int i=0 ; i<nxcells ; i++) {
        int itile, ii;
        itile = t_in[i];
        ii=ind_acc[itile];
        xgrid[n].per_intile[itile].input_parent_lon_indices[ii] = input_parent_lon_indices[i];
        xgrid[n].per_intile[itile].input_parent_lat_indices[ii] = input_parent_lat_indices[i];
        xgrid[n].per_intile[itile].output_parent_lon_indices[ii] = output_parent_lon_indices[i];
        xgrid[n].per_intile[itile].output_parent_lat_indices[ii] = output_parent_lat_indices[i];
        xgrid[n].per_intile[itile].xcell_area[ii] = xcell_area[i];
        if( opcode & CONSERVE_ORDER2) {
          xgrid[n].per_intile[itile].dcentroid_lon[ii] = xcell_centroid_lon[i];
          xgrid[n].per_intile[itile].dcentroid_lat[ii] = xcell_centroid_lat[i];
        }
        ind_acc[itile]++;
      }

      free(t_in) ; free(ind_acc);
      free(input_parent_lon_indices) ; free(input_parent_lat_indices);
      free(output_parent_lon_indices); free(output_parent_lat_indices);
      free(xcell_area); free(xcell_centroid_lon); free(xcell_centroid_lat);

    }//if file exists
  } //ntiles out

  printf("NOTE: Finish reading index and weight for conservative interpolation from file.\n");

} //end read_remap_file



/*******************************************************************************
 void do_scalar_conserve_interp( )
 doing conservative interpolation
*******************************************************************************/
void do_scalar_conserve_interp_acc(Xgrid_config *xgrid, int varid, int ntiles_input_grid, const Grid_config *input_grid,
                                   int ntiles_output_grid, const Grid_config *output_grid, const Field_config *field_in,
                                   Field_config *field_out, unsigned int opcode, int nz)
{
  int nx1, ny1, nx2, ny2, i1, j1, i2, j2, tile, n, m, i, j, n1, n2;
  int k, n0;
  int has_missing, halo, interp_method;
  int weight_exist;
  int cell_measures, cell_methods;
  double area, missing, di, dj, area_missing;
  double *out_area;
  int    *out_miss;
  double gsum_out;
  int monotonic;
  int target_grid;
  Monotone_config *monotone_data;

  gsum_out = 0;
  interp_method = field_in->var[varid].interp_method;
  halo = 0;
  monotonic = 0;
  if(interp_method == CONSERVE_ORDER2) {
    halo = 1;
    monotonic = opcode & MONOTONIC;
  }

  area_missing = field_in->var[varid].area_missing;
  has_missing = field_in->var[varid].has_missing;
  weight_exist = input_grid[0].weight_exist;
  cell_measures = field_in->var[varid].cell_measures;
  cell_methods = field_in->var[varid].cell_methods;
  target_grid = opcode & TARGET;
  if( field_in->var[varid].use_volume ) target_grid = 0;

  missing = -MAXVAL;
  if(has_missing) missing = field_in->var[varid].missing;

  if( nz>1 && has_missing ) mpp_error("conserve_interp: has_missing should be false when nz > 1");
  if( nz>1 && cell_measures ) mpp_error("conserve_interp: cell_measures should be false when nz > 1");
  if( nz>1 && cell_methods == CELL_METHODS_SUM ) mpp_error("conserve_interp: cell_methods should not be sum when nz > 1");
  /*  if( nz>1 && monotonic ) mpp_error("conserve_interp: monotonic should be false when nz > 1"); */

  if(monotonic) monotone_data = (Monotone_config *)malloc(ntiles_input_grid*sizeof(Monotone_config));

  for(m=0; m<ntiles_output_grid; m++) {
    nx2 = output_grid[m].nxc;
    ny2 = output_grid[m].nyc;
    out_area = (double *)malloc(nx2*ny2*nz*sizeof(double));
    out_miss = (int *)malloc(nx2*ny2*nz*sizeof(int));
    for(i=0; i<nx2*ny2*nz; i++) {
      field_out[m].data[i] = 0.0;
      out_area[i] = 0.0;
      out_miss[i] = 0;
    }
    if(interp_method == CONSERVE_ORDER1) {
      if(has_missing) {
	for(n=0; n<xgrid[m].nxcells; n++) {
	  i2   = xgrid[m].output_parent_lon_indices[n];
	  j2   = xgrid[m].output_parent_lat_indices[n];
	  i1   = xgrid[m].input_parent_lon_indices [n];
	  j1   = xgrid[m].input_parent_lat_indices [n];
	  tile = xgrid[m].input_parent_tile [n];
	  area = xgrid[m].xcell_area [n];
	  nx1  = input_grid[tile].nx;
	  ny1  = input_grid[tile].ny;
          if(weight_exist) area *= input_grid[tile].weight[j1*nx1+i1];
	  n1 = j1*nx1+i1;
	  n0 = j2*nx2+i2;

	  if( field_in[tile].data[n1] != missing ) {
	    if( cell_methods == CELL_METHODS_SUM )
	      area /= input_grid[tile].cell_area[n1];
            else if( cell_measures ) {
	      if(field_in[tile].area[n1] == area_missing) {
	         printf("name=%s,tile=%d,i1,j1=%d,%d,i2,j2=%d,%d\n",field_in->var[varid].name,tile,i1,j1,i2,j2);
	         mpp_error("conserve_interp: data is not missing but area is missing");
	      }
	      area *= (field_in[tile].area[n1]/input_grid[tile].cell_area[n1]);
	    }
  	    field_out[m].data[n0] += (field_in[tile].data[n1]*area);
            out_area[n0] += area;
	    out_miss[n0] = 1;
          }
        }
      }
      else {
	for(n=0; n<xgrid[m].nxcells; n++) {
	  i2   = xgrid[m].output_parent_lon_indices[n];
	  j2   = xgrid[m].output_parent_lat_indices[n];
	  i1   = xgrid[m].input_parent_lon_indices [n];
	  j1   = xgrid[m].input_parent_lat_indices [n];
	  tile = xgrid[m].input_parent_tile [n];
	  area = xgrid[m].xcell_area [n];
	  nx1  = input_grid[tile].nx;
	  ny1  = input_grid[tile].ny;
	  if(weight_exist) area *= input_grid[tile].weight[j1*nx1+i1];
	  for(k=0; k<nz; k++) {
	    n1 = k*nx1*ny1 + j1*nx1+i1;
	    n0 = k*nx2*ny2 + j2*nx2+i2;
	    if(  cell_methods == CELL_METHODS_SUM )
	      area /= input_grid[tile].cell_area[n1];
            else if( cell_measures )
	      area *= (field_in[tile].area[n1]/input_grid[tile].cell_area[n1]);
  	    field_out[m].data[n0] += (field_in[tile].data[n1]*area);
	    out_area[n0] += area;
	    out_miss[n0] = 1;
	  }
	}
      }
    }
    else if(monotonic) {
      int ii, jj;
      double f_bar;
      double *xdata;
      for(n=0; n<ntiles_input_grid; n++) {
	nx1 =  input_grid[n].nx;
	ny1 =  input_grid[n].ny;
	monotone_data[n].f_bar_max = (double *)malloc(nx1*ny1*sizeof(double));
	monotone_data[n].f_bar_min = (double *)malloc(nx1*ny1*sizeof(double));
	monotone_data[n].f_max     = (double *)malloc(nx1*ny1*sizeof(double));
	monotone_data[n].f_min     = (double *)malloc(nx1*ny1*sizeof(double));
	for(j=0; j<ny1; j++) for(i=0; i<nx1; i++) {
	  n1 = j*nx1+i;

	  monotone_data[n].f_bar_max[n1] = -MAXVAL;
	  monotone_data[n].f_bar_min[n1] = MAXVAL;
	  monotone_data[n].f_max[n1]     = -MAXVAL;
	  monotone_data[n].f_min[n1]     = MAXVAL;
	  n1 = j*nx1+i;
	  for(jj=j-1; jj<=j+1; jj++) for(ii=i-1; ii<=i+1; ii++) {
	    n2 = (jj+1)*(nx1+2)+ii+1;
	    if( field_in[n].data[n2] != missing ) {
	      if( field_in[n].data[n2] > monotone_data[n].f_bar_max[n1] ) monotone_data[n].f_bar_max[n1] = field_in[n].data[n2];
	      if( field_in[n].data[n2] < monotone_data[n].f_bar_min[n1] ) monotone_data[n].f_bar_min[n1] = field_in[n].data[n2];
	    }
	  }
	}
      }

      xdata = (double *)malloc(xgrid[m].nxcells*sizeof(double));
      for(n=0; n<xgrid[m].nxcells; n++) {
	i1   = xgrid[m].input_parent_lon_indices [n];
	j1   = xgrid[m].input_parent_lat_indices [n];
	di   = xgrid[m].dcentroid_lon[n];
	dj   = xgrid[m].dcentroid_lat[n];
	tile = xgrid[m].input_parent_tile [n];
	n1 = j1*nx1+i1;
        n2 = (j1+1)*(nx1+2)+i1+1;
	if( field_in[tile].data[n2] != missing ) {
	  if( field_in[tile].grad_mask[n1] ) { /* use zero gradient */
	    xdata[n] = field_in[tile].data[n2];
	  }
	  else {
	    xdata[n] = field_in[tile].data[n2]+field_in[tile].grad_x[n1]*di+field_in[tile].grad_y[n1]*dj;
	  }
	  if(monotonic) {
	    if( xdata[n] > monotone_data[tile].f_max[n1]) monotone_data[tile].f_max[n1] = xdata[n];
	    if( xdata[n] < monotone_data[tile].f_min[n1]) monotone_data[tile].f_min[n1] = xdata[n];
	  }
	}
	else
	  xdata[n] = missing;
      }

      /* get the global f_max and f_min */
      if(mpp_npes() >1) {
	for(n=0; n<ntiles_input_grid; n++) {
	  mpp_min_double(nx1*ny1, monotone_data[n].f_min);
	  mpp_max_double(nx1*ny1, monotone_data[n].f_max);
	}
      }

      /* adjust the exchange grid cell data to make it monotonic */
      for(n=0; n<xgrid[m].nxcells; n++) {
	i1   = xgrid[m].input_parent_lon_indices [n];
	j1   = xgrid[m].input_parent_lat_indices [n];
	tile = xgrid[m].input_parent_tile [n];
	n1 = j1*nx1+i1;
	n2 = (j1+1)*(nx1+2)+i1+1;
	f_bar = field_in[tile].data[n2];
	if(xdata[n] == missing) continue;

	if( monotone_data[tile].f_max[n1] > monotone_data[tile].f_bar_max[n1] ) {
	  /* z1l: Due to truncation error, we might get xdata[n] > f_bar_max[n1]. So
	     we allow some tolerance. What is the suitable tolerance? */
	  xdata[n] = f_bar + ((xdata[n]-f_bar)/(monotone_data[tile].f_max[n1]-f_bar))
	    * (monotone_data[tile].f_bar_max[n1]-f_bar);
	  if( xdata[n] > monotone_data[tile].f_bar_max[n1]) {
	    if(xdata[n] - monotone_data[tile].f_bar_max[n1] < TOLERANCE ) xdata[n] = monotone_data[tile].f_bar_max[n1];
	    if( xdata[n] > monotone_data[tile].f_bar_max[n1]) {
	      printf(" n = %d, n1 = %d, xdata = %f, f_bar_max=%f\n", n, n1, xdata[n], monotone_data[tile].f_bar_max[n1]);
	      mpp_error(" xdata is greater than f_bar_max ");
	    }
	  }
	}
	else if( monotone_data[tile].f_min[n1] < monotone_data[tile].f_bar_min[n1] ) {
	  /* z1l: Due to truncation error, we might get xdata[n] < f_bar_min[n1]. So
	     we allow some tolerance. What is the suitable tolerance? */
	  xdata[n] = f_bar + ((xdata[n]-f_bar)/(monotone_data[tile].f_min[n1]-f_bar)) * (monotone_data[tile].f_bar_min[n1]-f_bar);
	  if( xdata[n] < monotone_data[tile].f_bar_min[n1]) {
	    if(monotone_data[tile].f_bar_min[n1] - xdata[n]< TOLERANCE ) xdata[n] = monotone_data[tile].f_bar_min[n1];
	    if( xdata[n] < monotone_data[tile].f_bar_min[n1]) {
	      printf(" n = %d, n1 = %d, xdata = %f, f_bar_min=%f\n", n, n1, xdata[n], monotone_data[tile].f_bar_min[n1]);
	      mpp_error(" xdata is less than f_bar_min ");
	    }
	  }
	}
      }
      for(n=0; n<ntiles_input_grid; n++) {
	free(monotone_data[n].f_bar_max);
	free(monotone_data[n].f_bar_min);
	free(monotone_data[n].f_max);
	free(monotone_data[n].f_min);
      }

      /* remap onto destination grid */
      for(n=0; n<xgrid[m].nxcells; n++) {
	i2   = xgrid[m].output_parent_lon_indices[n];
	j2   = xgrid[m].output_parent_lat_indices[n];
	i1   = xgrid[m].input_parent_lon_indices [n];
	j1   = xgrid[m].input_parent_lat_indices [n];
	tile = xgrid[m].input_parent_tile [n];
	area = xgrid[m].xcell_area [n];
	if(xdata[n] == missing) continue;
	if(weight_exist) area *= input_grid[tile].weight[j1*nx1+i1];
	n1 = j1*nx1+i1;
	n0 = j2*nx2+i2;
	if( cell_methods == CELL_METHODS_SUM )
	  area /= input_grid[tile].cell_area[n1];
	else if( cell_measures )
	  area *= (field_in[tile].area[n1]/input_grid[tile].cell_area[n1]);
	field_out[m].data[n0] += xdata[n]*area;
	out_area[n0] += area;
      }
      free(xdata);
    }
    else {
      if(has_missing) {
	for(n=0; n<xgrid[m].nxcells; n++) {
	  i2   = xgrid[m].output_parent_lon_indices[n];
	  j2   = xgrid[m].output_parent_lat_indices[n];
	  i1   = xgrid[m].input_parent_lon_indices [n];
	  j1   = xgrid[m].input_parent_lat_indices [n];
	  di   = xgrid[m].dcentroid_lon[n];
	  dj   = xgrid[m].dcentroid_lat[n];
	  tile = xgrid[m].input_parent_tile [n];
	  area = xgrid[m].xcell_area [n];
	  nx1  = input_grid[tile].nx;
	  ny1  = input_grid[tile].ny;

	  if(weight_exist) area *= input_grid[tile].weight[j1*nx1+i1];
	  n2 = (j1+1)*(nx1+2)+i1+1;
	  n0 = j2*nx2+i2;
	  if( field_in[tile].data[n2] != missing ) {
	    n1 = j1*nx1+i1;
	    n0 = j2*nx2+i2;
            if( cell_methods == CELL_METHODS_SUM )
	      area /= input_grid[tile].cell_area[n1];
            else if( cell_measures ) {
	      if(field_in[tile].area[n1] == area_missing) {
                printf("name=%s,tile=%d,i1,j1=%d,%d,i2,j2=%d,%d\n",field_in->var[varid].name,tile,i1,j1,i2,j2);
	        mpp_error("conserve_interp: data is not missing but area is missing");
              }
	      area *= (field_in[tile].area[n1]/input_grid[tile].cell_area[n1]);
            }
	    if(field_in[tile].grad_mask[n1]) { /* use zero gradient */
	      field_out[m].data[n0] += field_in[tile].data[n2]*area;
	    }
	    else {
	      field_out[m].data[n0] += (field_in[tile].data[n2]+field_in[tile].grad_x[n1]*di
					+field_in[tile].grad_y[n1]*dj)*area;
	    }
	    out_area[n0] += area;
	    out_miss[n0] = 1;
	  }
	}
      }
      else {
	for(n=0; n<xgrid[m].nxcells; n++) {
	  i2   = xgrid[m].output_parent_lon_indices[n];
	  j2   = xgrid[m].output_parent_lat_indices[n];
	  i1   = xgrid[m].input_parent_lon_indices [n];
	  j1   = xgrid[m].input_parent_lat_indices [n];
	  di   = xgrid[m].dcentroid_lon[n];
	  dj   = xgrid[m].dcentroid_lat[n];
	  tile = xgrid[m].input_parent_tile [n];
	  area = xgrid[m].xcell_area [n];

	  nx1  = input_grid[tile].nx;
	  ny1  = input_grid[tile].ny;
	  if(weight_exist) area *= input_grid[tile].weight[j1*nx1+i1];
	  for(k=0; k<nz; k++) {
	    n0 = k*nx2*ny2 + j2*nx2+i2;
	    n1 = k*nx1*ny1+j1*nx1+i1;
	    n2 = k*(nx1+2)*(ny1+2)+(j1+1)*(nx1+2)+i1+1;
	    if( cell_methods == CELL_METHODS_SUM )
	      area /= input_grid[tile].cell_area[n1];
	    else if( cell_measures )
	      area *= (field_in[tile].area[n1]/input_grid[tile].cell_area[n1]);
	    field_out[m].data[n0] += (field_in[tile].data[n2]+field_in[tile].grad_x[n1]*di
				      +field_in[tile].grad_y[n1]*dj)*area;
	    out_area[n0] += area;
            out_miss[n0] = 1;
	  }
	}
      }
    }

    if(opcode & CHECK_CONSERVE) {
      for(i=0; i<nx2*ny2*nz; i++) {
	if(out_area[i] > 0) gsum_out += field_out[m].data[i];
      }
    }

    if ( cell_methods == CELL_METHODS_SUM ) {
      for(i=0; i<nx2*ny2*nz; i++) {
        if(out_area[i] == 0) {
          if(out_miss[i] == 0)
            for(k=0; k<nz; k++) field_out[m].data[k*nx2*ny2+i] = missing;
          else
            for(k=0; k<nz; k++) field_out[m].data[k*nx2*ny2+i] = 0.0;
        }
      }
    }
    else {
      for(i=0; i<nx2*ny2*nz; i++) {
	if(out_area[i] > 0)
	  field_out[m].data[i] /= out_area[i];
	else if(out_miss[i] == 1)
	  field_out[m].data[i] = 0.0;
	else
	  field_out[m].data[i] = missing;
      }

      if( (target_grid) ) {
	for(i=0; i<nx2*ny2; i++) out_area[i] = 0.0;
	for(n=0; n<xgrid[m].nxcells; n++) {
	  i2   = xgrid[m].output_parent_lon_indices[n];
	  j2   = xgrid[m].output_parent_lat_indices[n];
	  i1   = xgrid[m].input_parent_lon_indices [n];
	  j1   = xgrid[m].input_parent_lat_indices [n];
	  tile = xgrid[m].input_parent_tile [n];
	  area = xgrid[m].xcell_area [n];
	  nx1  = input_grid[tile].nx;
	  ny1  = input_grid[tile].ny;
	  n0 = j2*nx2+i2;
	  n1 = j1*nx1+i1;
	  if(cell_measures )
	    out_area[n0] += (area*field_in[tile].area[n1]/input_grid[tile].cell_area[n1]);
	  else
	    out_area[n0] += area;
	}
	for(i=0; i<nx2*ny2*nz; i++) {
	  if(field_out[m].data[i] != missing) {
	    i2 = i%(nx2*ny2);
	    field_out[m].data[i] *=  (out_area[i2]/output_grid[m].cell_area[i2]);
	  }
	}
      }
    }

    free(out_area);
    free(out_miss);
  }


  /* conservation check if needed */
  if(opcode & CHECK_CONSERVE) {
    double gsum_in, dd;
    gsum_in = 0;
    for(n=0; n<ntiles_input_grid; n++) {

      nx1  = input_grid[n].nx;
      ny1  = input_grid[n].ny;


      if( cell_measures ) {
        for(j=0; j<ny1; j++) for(i=0; i<nx1; i++) {
  	  dd = field_in[n].data[(j+halo)*(nx1+2*halo)+i+halo];
	  if(dd != missing) gsum_in += dd*field_in[n].area[j*nx1+i];
        }
      }
      else if ( cell_methods == CELL_METHODS_SUM ) {
        for(j=0; j<ny1; j++) for(i=0; i<nx1; i++) {
	  dd = field_in[n].data[(j+halo)*(nx1+2*halo)+i+halo];
	  if(dd != missing) gsum_in += dd;
        }
      }
      else {
        for(k=0; k<nz; k++) for(j=0; j<ny1; j++) for(i=0; i<nx1; i++) {
	  dd = field_in[n].data[k*(nx1+2*halo)*(ny1+2*halo)+(j+halo)*(nx1+2*halo)+i+halo];
	  if(dd != missing) gsum_in += dd*input_grid[n].cell_area[j*nx1+i];
        }
      }
    }
    mpp_sum_double(1, &gsum_out);

    if(mpp_pe() == mpp_root_pe()) printf("the flux(data*area) sum of %s: input = %g, output = %g, diff = %g. \n",
					 field_in->var[varid].name, gsum_in, gsum_out, gsum_out-gsum_in);

  }


}; /* do_scalar_conserve_interp */


/*******************************************************************************
 void do_vector_conserve_interp( )
 doing conservative interpolation
*******************************************************************************/
void do_vector_conserve_interp_acc(Xgrid_config *xgrid, int varid, int ntiles_input_grid, const Grid_config *input_grid, int ntiles_output_grid,
                               const Grid_config *output_grid, const Field_config *u_in,  const Field_config *v_in,
                               Field_config *u_out, Field_config *v_out, unsigned int opcode)
{
  int          nx1, ny1, nx2, ny2, i1, j1, i2, j2, tile, n, m, i;
  double       area, missing, tmp_x, tmp_y;
  double       *out_area;

  missing = u_in->var[varid].missing;
  /* first rotate input data */
  for(n = 0; n < ntiles_input_grid; n++) {
    if(input_grid[n].rotate) {
      nx1 = input_grid[n].nx;
      ny1 = input_grid[n].ny;
      for(i=0; i<nx1*ny1; i++) {
	tmp_x = u_in[n].data[i];
	tmp_y = v_in[n].data[i];
	if( tmp_x != missing && tmp_y != missing) {
	  u_in[n].data[i] = tmp_x * input_grid[n].cosrot[i] - tmp_y * input_grid[n].sinrot[i];
	  v_in[n].data[i] = tmp_x * input_grid[n].sinrot[i] + tmp_y * input_grid[n].cosrot[i];
	}
      }
    }
  }

  for(m=0; m<ntiles_output_grid; m++) {
    nx2 = output_grid[m].nxc;
    ny2 = output_grid[m].nyc;
    out_area = (double *)malloc(nx2*ny2*sizeof(double));

    for(i=0; i<nx2*ny2; i++) {
      u_out[m].data[i] = 0.0;
      v_out[m].data[i] = 0.0;
    }
    for(i=0; i<nx2*ny2; i++) out_area[i] = 0.0;

    for(n=0; n<xgrid[m].nxcells; n++) {
      i2   = xgrid[m].output_parent_lon_indices[n];
      j2   = xgrid[m].output_parent_lat_indices[n];
      i1   = xgrid[m].input_parent_lon_indices [n];
      j1   = xgrid[m].input_parent_lat_indices [n];
      tile = xgrid[m].input_parent_tile [n];
      area = xgrid[m].xcell_area [n];
      nx1  = input_grid[tile].nx;
      ny1  = input_grid[tile].ny;
      tmp_x = u_in[tile].data[j1*nx1+i1];
      tmp_y = v_in[tile].data[j1*nx1+i1];
      if( tmp_x != missing && tmp_y != missing ) {
	u_out[m].data[j2*nx2+i2] += u_in[tile].data[j1*nx1+i1]*area;
	v_out[m].data[j2*nx2+i2] += v_in[tile].data[j1*nx1+i1]*area;
	out_area[j2*nx2+i2] += area;
      }
    }
    if(opcode & TARGET) {
      for(i=0; i<nx2*ny2; i++) {
	if(out_area[i] > 0) {
	  u_out[m].data[i] /= output_grid[m].area[i];
	  v_out[m].data[i] /= output_grid[m].area[i];
	}
	else {
	  u_out[m].data[i] = missing;
	  v_out[m].data[i] = missing;
	}
      }
    }
    else {
      for(i=0; i<nx2*ny2; i++) {
	if(out_area[i] > 0) {
	  u_out[m].data[i] /= out_area[i];
	  v_out[m].data[i] /= out_area[i];
	}
	else {
	  u_out[m].data[i] = missing;
	  v_out[m].data[i] = missing;
	}
      }
    }
    /* rotate the data if needed */
    if(output_grid[m].rotate) {
      for(i=0; i<nx2*ny2; i++) {
	tmp_x = u_out[m].data[i];
	tmp_y = v_out[m].data[i];
	if( tmp_x != missing && tmp_y != missing) {
	  u_out[m].data[i] =  tmp_x * output_grid[m].cosrot[i] + tmp_y * output_grid[m].sinrot[i];
	  v_out[m].data[i] = -tmp_x * output_grid[m].sinrot[i] + tmp_y * output_grid[m].cosrot[i];
	}
      }
    }
    free(out_area);
  }

}; /* do_vector_conserve_interp */
