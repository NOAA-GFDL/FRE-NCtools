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
#include <time.h>
#include "constant.h"
#include "globals.h"
#include "create_xgrid.h"
#include "create_xgrid_acc.h"
#include "create_xgrid_util.h"
#include "mosaic_util.h"
#include "conserve_interp.h"
#include "conserve_interp_util.h"
#include "fregrid_util.h"
#include "mpp.h"
#include "mpp_io.h"
#include "read_mosaic.h"

#define DEBUG 0

#define  AREA_RATIO (1.e-3)
#define  MAXVAL (1.e20)
#define  TOLERANCE  (1.e-10)
/*******************************************************************************
  void setup_conserve_interp
  Setup the interpolation weight for conservative interpolation
*******************************************************************************/
void setup_conserve_interp(int ntiles_in, const Grid_config *grid_in, int ntiles_out,
         Grid_config *grid_out, Interp_config *interp, unsigned int opcode)
{
  int    n, m, i, ii, jj, nx_in, ny_in, nx_out, ny_out, tile;
  size_t nxgrid, nxgrid_prev;
  int mxxgrid, zero=0;

  CellStruct *cell_in;

  double time_nxgrid=0;
  clock_t time_start, time_end;

  Minmaxavg_lists out_minmaxavg_lists;
  out_minmaxavg_lists.lon_list=NULL;
  out_minmaxavg_lists.lat_list=NULL;
  out_minmaxavg_lists.lon_min_list=NULL;
  out_minmaxavg_lists.lat_min_list=NULL;
  out_minmaxavg_lists.lon_max_list=NULL;
  out_minmaxavg_lists.lat_max_list=NULL;
  out_minmaxavg_lists.n_list=NULL;
  out_minmaxavg_lists.lon_avg=NULL;

  if( opcode & READ) {
    read_remap_file(ntiles_in,ntiles_out, grid_out, interp, opcode);
  }
  else {

    //only needed for order2?
    cell_in    = (CellStruct *)malloc(ntiles_in * sizeof(CellStruct));
    for(m=0; m<ntiles_in; m++) {
      nx_in = grid_in[m].nx;
      ny_in = grid_in[m].ny;
      cell_in[m].area = (double *)calloc(nx_in*ny_in,sizeof(double));
      cell_in[m].clon = (double *)calloc(nx_in*ny_in,sizeof(double));
      cell_in[m].clat = (double *)calloc(nx_in*ny_in,sizeof(double));
    }

    //START NTILES_OUT
    for(n=0; n<ntiles_out; n++) {

      nx_out = grid_out[n].nxc;
      ny_out = grid_out[n].nyc;
      interp[n].nxgrid = 0;

#pragma acc enter data copyin(grid_out[n].lonc[0:(nx_out+1)*(ny_out+1)], \
                              grid_out[n].latc[0:(nx_out+1)*(ny_out+1)])


      //allocate memory for the lists
      malloc_minmaxavg_lists(nx_out*ny_out, &out_minmaxavg_lists);

#define MAX_V 8
#pragma acc enter data create(out_minmaxavg_lists)
#pragma acc enter data create(out_minmaxavg_lists.lon_list[0:MAX_V*nx_out*ny_out], \
                              out_minmaxavg_lists.lat_list[0:MAX_V*nx_out*ny_out], \
                              out_minmaxavg_lists.lon_min_list[0:nx_out*ny_out], \
                              out_minmaxavg_lists.lon_max_list[0:nx_out*ny_out], \
                              out_minmaxavg_lists.lat_min_list[0:nx_out*ny_out], \
                              out_minmaxavg_lists.lat_max_list[0:nx_out*ny_out], \
                              out_minmaxavg_lists.n_list[0:nx_out*ny_out], \
                              out_minmaxavg_lists.lon_avg[0:nx_out*ny_out] )

      //compute the list values
      get_minmaxavg_lists(nx_out, ny_out, grid_out[n].lonc, grid_out[n].latc, &out_minmaxavg_lists);

      //START NTILES_IN
      for(m=0; m<ntiles_in; m++) {

        nx_in = grid_in[m].nx;
        ny_in = grid_in[m].ny;

        if(opcode & GREAT_CIRCLE) {
          do_great_circle(n, m, grid_in, grid_out, interp, opcode) ;
        }
        else {
          if(opcode & CONSERVE_ORDER1) {
            do_create_xgrid_order1(n, m, grid_in, grid_out, interp, opcode);
          }
          else if(opcode & CONSERVE_ORDER2) {
            do_create_xgrid_order2(n, m, grid_in, grid_out, &out_minmaxavg_lists, cell_in, interp, opcode);
          }
          else
            mpp_error("conserve_interp: interp_method should be CONSERVE_ORDER1 or CONSERVE_ORDER2");
        } // opcode GREAT_CIRCLE or CONSERVE_ORDERs
      } // ntiles_in

      malloc_minmaxavg_lists(zero, &out_minmaxavg_lists);
#pragma acc exit data delete(out_minmaxavg_lists)
#pragma acc exit data delete(grid_out[n].latc, grid_out[n].lonc)

    } // ntimes_out

    if(DEBUG) print_time("time_nxgrid", time_nxgrid);

    if(opcode & CONSERVE_ORDER2) {
      /* subtrack the grid_in clon and clat to get the distance between xgrid and grid_in */
      for(n=0; n<ntiles_in; n++) {
        double x1_in[50], y1_in[50], lon_in_avg, clon, clat;
        int    j, n0, n1, n2, n3, n1_in;
        /* calcualte cell area */
        nx_in = grid_in[n].nx;
        ny_in = grid_in[n].ny;
        for(j=0; j<ny_in; j++) for(i=0; i<nx_in; i++) {
            ii = j*nx_in + i;
            if(cell_in[n].area[ii] > 0) {
              if( fabs(cell_in[n].area[ii]-grid_in[n].cell_area[ii])/grid_in[n].cell_area[ii] < AREA_RATIO ) {
                cell_in[n].clon[ii] /= cell_in[n].area[ii];
                cell_in[n].clat[ii] /= cell_in[n].area[ii];
              }
              else {
                n0 = j*(nx_in+1)+i;       n1 = j*(nx_in+1)+i+1;
                n2 = (j+1)*(nx_in+1)+i+1; n3 = (j+1)*(nx_in+1)+i;
                x1_in[0] = grid_in[n].lonc[n0]; y1_in[0] = grid_in[n].latc[n0];
                x1_in[1] = grid_in[n].lonc[n1]; y1_in[1] = grid_in[n].latc[n1];
                x1_in[2] = grid_in[n].lonc[n2]; y1_in[2] = grid_in[n].latc[n2];
                x1_in[3] = grid_in[n].lonc[n3]; y1_in[3] = grid_in[n].latc[n3];
                n1_in = fix_lon(x1_in, y1_in, 4, M_PI);
                lon_in_avg = avgval_double(n1_in, x1_in);
                clon = poly_ctrlon(x1_in, y1_in, n1_in, lon_in_avg);
                clat = poly_ctrlat (x1_in, y1_in, n1_in );
                cell_in[n].clon[ii] = clon/grid_in[n].cell_area[ii];
                cell_in[n].clat[ii] = clat/grid_in[n].cell_area[ii];
              }
            }
          }
      }
      for(n=0; n<ntiles_out; n++) {
        for(i=0; i<interp[n].nxgrid; i++) {
          tile = interp[n].t_in[i];
          ii   = interp[n].j_in[i] * grid_in[tile].nx + interp[n].i_in[i];
          interp[n].di_in[i] -= cell_in[tile].clon[ii];
          interp[n].dj_in[i] -= cell_in[tile].clat[ii];
        }
      }

      /* free the memory */
      for(n=0; n<ntiles_in; n++) {
        free(cell_in[n].area);
        free(cell_in[n].clon);
        free(cell_in[n].clat);
      }
      free(cell_in);
    }
    if( opcode & WRITE) { /* write out remapping information */
      for(n=0; n<ntiles_out; n++) {
        int nxgrid;

        nxgrid = interp[n].nxgrid;
        mpp_sum_int(1, &nxgrid);
        if(nxgrid > 0) {
          size_t start[4], nwrite[4];
          int    fid, dim_string, dim_ncells, dim_two, dims[4];
          int    id_xgrid_area, id_tile1_dist;
          int    id_tile1_cell, id_tile2_cell, id_tile1;
          int    *gdata_int, *ldata_int;
          double *gdata_dbl;

          fid = mpp_open( interp[n].remap_file, MPP_WRITE);
          dim_string = mpp_def_dim(fid, "string", STRING);
          dim_ncells = mpp_def_dim(fid, "ncells", nxgrid);
          dim_two    = mpp_def_dim(fid, "two", 2);
          dims[0] = dim_ncells; dims[1] = dim_two;
          id_tile1      = mpp_def_var(fid, "tile1",      NC_INT, 1, &dim_ncells, 1,
                                      "standard_name", "tile_number_in_mosaic1");
          id_tile1_cell = mpp_def_var(fid, "tile1_cell", NC_INT, 2, dims, 1,
                                      "standard_name", "parent_cell_indices_in_mosaic1");
          id_tile2_cell = mpp_def_var(fid, "tile2_cell", NC_INT, 2, dims, 1,
                                      "standard_name", "parent_cell_indices_in_mosaic2");
          id_xgrid_area = mpp_def_var(fid, "xgrid_area", NC_DOUBLE, 1, &dim_ncells, 2,
                                      "standard_name", "exchange_grid_area", "units", "m2");
          if(opcode & CONSERVE_ORDER2) id_tile1_dist = mpp_def_var(fid, "tile1_distance", NC_DOUBLE, 2, dims, 1,
                                                                   "standard_name", "distance_from_parent1_cell_centroid");
          mpp_end_def(fid);
          for(i=0; i<4; i++) {
            start[i] = 0; nwrite[i] = 1;
          }
          nwrite[0] = nxgrid;
          gdata_int = (int *)malloc(nxgrid*sizeof(int));
          if(interp[n].nxgrid>0) ldata_int = (int *)malloc(interp[n].nxgrid*sizeof(int));
          mpp_gather_field_int(interp[n].nxgrid, interp[n].t_in, gdata_int);
          for(i=0; i<nxgrid; i++) gdata_int[i]++;
          mpp_put_var_value(fid, id_tile1, gdata_int);

          mpp_gather_field_int(interp[n].nxgrid, interp[n].i_in, gdata_int);
          for(i=0; i<nxgrid; i++) gdata_int[i]++;
          mpp_put_var_value_block(fid, id_tile1_cell, start, nwrite, gdata_int);

          for(i=0; i<interp[n].nxgrid; i++) ldata_int[i] = interp[n].i_out[i] + grid_out[n].isc + 1;
          mpp_gather_field_int(interp[n].nxgrid, ldata_int, gdata_int);
          mpp_put_var_value_block(fid, id_tile2_cell, start, nwrite, gdata_int);

          mpp_gather_field_int(interp[n].nxgrid, interp[n].j_in, gdata_int);
          for(i=0; i<nxgrid; i++) gdata_int[i]++;
          start[1] = 1;
          mpp_put_var_value_block(fid, id_tile1_cell, start, nwrite, gdata_int);

          for(i=0; i<interp[n].nxgrid; i++) ldata_int[i] = interp[n].j_out[i] + grid_out[n].jsc + 1;
          mpp_gather_field_int(interp[n].nxgrid, ldata_int, gdata_int);
          mpp_put_var_value_block(fid, id_tile2_cell, start, nwrite, gdata_int);

          free(gdata_int);
          if(interp[n].nxgrid>0)free(ldata_int);

          gdata_dbl = (double *)malloc(nxgrid*sizeof(double));
          mpp_gather_field_double(interp[n].nxgrid, interp[n].area, gdata_dbl);
          mpp_put_var_value(fid, id_xgrid_area, gdata_dbl);

          if(opcode & CONSERVE_ORDER2) {
            start[1] = 0;
            mpp_gather_field_double(interp[n].nxgrid, interp[n].di_in, gdata_dbl);
            mpp_put_var_value_block(fid, id_tile1_dist, start, nwrite, gdata_dbl);
            start[1] = 1;
            mpp_gather_field_double(interp[n].nxgrid, interp[n].dj_in, gdata_dbl);
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
    double *area2;

    /* sum over exchange grid to get the area of grid_in */
    nx1  = grid_out[0].nxc;
    ny1  = grid_out[0].nyc;

    area2 = (double *)malloc(nx1*ny1*sizeof(double));

    for(n=0; n<ntiles_out; n++) {
      for(i=0; i<nx1*ny1; i++) area2[i] = 0;
      for(i=0; i<interp[n].nxgrid; i++) {
        ii = interp[n].j_out[i]*nx1 + interp[n].i_out[i];
        area2[ii] +=  interp[n].area[i];
      }
      max_ratio = 0;
      max_i = 0;
      max_j = 0;
      /* comparing area1 and area2 */
      for(j=0; j<ny1; j++) for(i=0; i<nx1; i++) {
  ii = j*nx1+i;
  ratio_change = fabs(grid_out[n].cell_area[ii]-area2[ii])/grid_out[n].cell_area[ii];
  if(ratio_change > max_ratio) {
    max_ratio = ratio_change;
    max_i = i;
    max_j = j;
  }
  if( ratio_change > 1.e-4 ) {
    printf("(i,j)=(%d,%d), change = %g, area1=%g, area2=%g\n", i, j, ratio_change, grid_out[n].cell_area[ii],area2[ii]);
  }
      }
      ii = max_j*nx1+max_i;
      printf("The maximum ratio change at (%d,%d) = %g, area1=%g, area2=%g\n", max_i, max_j, max_ratio, grid_out[n].cell_area[ii],area2[ii]);

    }

    free(area2);

  }

}; /* setup_conserve_interp */


/*******************************************************************************
 void do_scalar_conserve_interp( )
 doing conservative interpolation
*******************************************************************************/
void do_scalar_conserve_interp(Interp_config *interp, int varid, int ntiles_in, const Grid_config *grid_in,
             int ntiles_out, const Grid_config *grid_out, const Field_config *field_in,
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
  weight_exist = grid_in[0].weight_exist;
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

  if(monotonic) monotone_data = (Monotone_config *)malloc(ntiles_in*sizeof(Monotone_config));

  for(m=0; m<ntiles_out; m++) {
    nx2 = grid_out[m].nxc;
    ny2 = grid_out[m].nyc;
    out_area = (double *)malloc(nx2*ny2*nz*sizeof(double));
    out_miss = (int *)malloc(nx2*ny2*nz*sizeof(int));
    for(i=0; i<nx2*ny2*nz; i++) {
      field_out[m].data[i] = 0.0;
      out_area[i] = 0.0;
      out_miss[i] = 0;
    }
    if(interp_method == CONSERVE_ORDER1) {
      if(has_missing) {
  for(n=0; n<interp[m].nxgrid; n++) {
    i2   = interp[m].i_out[n];
    j2   = interp[m].j_out[n];
    i1   = interp[m].i_in [n];
    j1   = interp[m].j_in [n];
    tile = interp[m].t_in [n];
    area = interp[m].area [n];
    nx1  = grid_in[tile].nx;
    ny1  = grid_in[tile].ny;
          if(weight_exist) area *= grid_in[tile].weight[j1*nx1+i1];
    n1 = j1*nx1+i1;
    n0 = j2*nx2+i2;

    if( field_in[tile].data[n1] != missing ) {
      if( cell_methods == CELL_METHODS_SUM )
        area /= grid_in[tile].cell_area[n1];
            else if( cell_measures ) {
        if(field_in[tile].area[n1] == area_missing) {
           printf("name=%s,tile=%d,i1,j1=%d,%d,i2,j2=%d,%d\n",field_in->var[varid].name,tile,i1,j1,i2,j2);
           mpp_error("conserve_interp: data is not missing but area is missing");
        }
        area *= (field_in[tile].area[n1]/grid_in[tile].cell_area[n1]);
      }
        field_out[m].data[n0] += (field_in[tile].data[n1]*area);
            out_area[n0] += area;
      out_miss[n0] = 1;
          }
        }
      }
      else {
  for(n=0; n<interp[m].nxgrid; n++) {
    i2   = interp[m].i_out[n];
    j2   = interp[m].j_out[n];
    i1   = interp[m].i_in [n];
    j1   = interp[m].j_in [n];
    tile = interp[m].t_in [n];
    area = interp[m].area [n];
    nx1  = grid_in[tile].nx;
    ny1  = grid_in[tile].ny;
    if(weight_exist) area *= grid_in[tile].weight[j1*nx1+i1];
    for(k=0; k<nz; k++) {
      n1 = k*nx1*ny1 + j1*nx1+i1;
      n0 = k*nx2*ny2 + j2*nx2+i2;
      if(  cell_methods == CELL_METHODS_SUM )
        area /= grid_in[tile].cell_area[n1];
            else if( cell_measures )
        area *= (field_in[tile].area[n1]/grid_in[tile].cell_area[n1]);
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
      for(n=0; n<ntiles_in; n++) {
  nx1 =  grid_in[n].nx;
  ny1 =  grid_in[n].ny;
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

      xdata = (double *)malloc(interp[m].nxgrid*sizeof(double));
      for(n=0; n<interp[m].nxgrid; n++) {
  i1   = interp[m].i_in [n];
  j1   = interp[m].j_in [n];
  di   = interp[m].di_in[n];
  dj   = interp[m].dj_in[n];
  tile = interp[m].t_in [n];
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
  for(n=0; n<ntiles_in; n++) {
    mpp_min_double(nx1*ny1, monotone_data[n].f_min);
    mpp_max_double(nx1*ny1, monotone_data[n].f_max);
  }
      }

      /* adjust the exchange grid cell data to make it monotonic */
      for(n=0; n<interp[m].nxgrid; n++) {
  i1   = interp[m].i_in [n];
  j1   = interp[m].j_in [n];
  tile = interp[m].t_in [n];
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
      for(n=0; n<ntiles_in; n++) {
  free(monotone_data[n].f_bar_max);
  free(monotone_data[n].f_bar_min);
  free(monotone_data[n].f_max);
  free(monotone_data[n].f_min);
      }

      /* remap onto destination grid */
      for(n=0; n<interp[m].nxgrid; n++) {
  i2   = interp[m].i_out[n];
  j2   = interp[m].j_out[n];
  i1   = interp[m].i_in [n];
  j1   = interp[m].j_in [n];
  tile = interp[m].t_in [n];
  area = interp[m].area [n];
  if(xdata[n] == missing) continue;
  if(weight_exist) area *= grid_in[tile].weight[j1*nx1+i1];
  n1 = j1*nx1+i1;
  n0 = j2*nx2+i2;
  if( cell_methods == CELL_METHODS_SUM )
    area /= grid_in[tile].cell_area[n1];
  else if( cell_measures )
    area *= (field_in[tile].area[n1]/grid_in[tile].cell_area[n1]);
  field_out[m].data[n0] += xdata[n]*area;
  out_area[n0] += area;
      }
      free(xdata);
    }
    else {
      if(has_missing) {
  for(n=0; n<interp[m].nxgrid; n++) {
    i2   = interp[m].i_out[n];
    j2   = interp[m].j_out[n];
    i1   = interp[m].i_in [n];
    j1   = interp[m].j_in [n];
    di   = interp[m].di_in[n];
    dj   = interp[m].dj_in[n];
    tile = interp[m].t_in [n];
    area = interp[m].area [n];
    nx1  = grid_in[tile].nx;
    ny1  = grid_in[tile].ny;

    if(weight_exist) area *= grid_in[tile].weight[j1*nx1+i1];
    n2 = (j1+1)*(nx1+2)+i1+1;
    n0 = j2*nx2+i2;
    if( field_in[tile].data[n2] != missing ) {
      n1 = j1*nx1+i1;
      n0 = j2*nx2+i2;
            if( cell_methods == CELL_METHODS_SUM )
        area /= grid_in[tile].cell_area[n1];
            else if( cell_measures ) {
        if(field_in[tile].area[n1] == area_missing) {
                printf("name=%s,tile=%d,i1,j1=%d,%d,i2,j2=%d,%d\n",field_in->var[varid].name,tile,i1,j1,i2,j2);
          mpp_error("conserve_interp: data is not missing but area is missing");
              }
        area *= (field_in[tile].area[n1]/grid_in[tile].cell_area[n1]);
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
  for(n=0; n<interp[m].nxgrid; n++) {
    i2   = interp[m].i_out[n];
    j2   = interp[m].j_out[n];
    i1   = interp[m].i_in [n];
    j1   = interp[m].j_in [n];
    di   = interp[m].di_in[n];
    dj   = interp[m].dj_in[n];
    tile = interp[m].t_in [n];
    area = interp[m].area [n];

    nx1  = grid_in[tile].nx;
    ny1  = grid_in[tile].ny;
    if(weight_exist) area *= grid_in[tile].weight[j1*nx1+i1];
    for(k=0; k<nz; k++) {
      n0 = k*nx2*ny2 + j2*nx2+i2;
      n1 = k*nx1*ny1+j1*nx1+i1;
      n2 = k*(nx1+2)*(ny1+2)+(j1+1)*(nx1+2)+i1+1;
      if( cell_methods == CELL_METHODS_SUM )
        area /= grid_in[tile].cell_area[n1];
      else if( cell_measures )
        area *= (field_in[tile].area[n1]/grid_in[tile].cell_area[n1]);
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
  for(n=0; n<interp[m].nxgrid; n++) {
    i2   = interp[m].i_out[n];
    j2   = interp[m].j_out[n];
    i1   = interp[m].i_in [n];
    j1   = interp[m].j_in [n];
    tile = interp[m].t_in [n];
    area = interp[m].area [n];
    nx1  = grid_in[tile].nx;
    ny1  = grid_in[tile].ny;
    n0 = j2*nx2+i2;
    n1 = j1*nx1+i1;
    if(cell_measures )
      out_area[n0] += (area*field_in[tile].area[n1]/grid_in[tile].cell_area[n1]);
    else
      out_area[n0] += area;
  }
  for(i=0; i<nx2*ny2*nz; i++) {
    if(field_out[m].data[i] != missing) {
      i2 = i%(nx2*ny2);
      field_out[m].data[i] *=  (out_area[i2]/grid_out[m].cell_area[i2]);
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
    for(n=0; n<ntiles_in; n++) {

      nx1  = grid_in[n].nx;
      ny1  = grid_in[n].ny;


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
    if(dd != missing) gsum_in += dd*grid_in[n].cell_area[j*nx1+i];
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
void do_vector_conserve_interp(Interp_config *interp, int varid, int ntiles_in, const Grid_config *grid_in, int ntiles_out,
                               const Grid_config *grid_out, const Field_config *u_in,  const Field_config *v_in,
                               Field_config *u_out, Field_config *v_out, unsigned int opcode)
{
  int          nx1, ny1, nx2, ny2, i1, j1, i2, j2, tile, n, m, i;
  double       area, missing, tmp_x, tmp_y;
  double       *out_area;

  missing = u_in->var[varid].missing;
  /* first rotate input data */
  for(n = 0; n < ntiles_in; n++) {
    if(grid_in[n].rotate) {
      nx1 = grid_in[n].nx;
      ny1 = grid_in[n].ny;
      for(i=0; i<nx1*ny1; i++) {
  tmp_x = u_in[n].data[i];
  tmp_y = v_in[n].data[i];
  if( tmp_x != missing && tmp_y != missing) {
    u_in[n].data[i] = tmp_x * grid_in[n].cosrot[i] - tmp_y * grid_in[n].sinrot[i];
    v_in[n].data[i] = tmp_x * grid_in[n].sinrot[i] + tmp_y * grid_in[n].cosrot[i];
  }
      }
    }
  }

  for(m=0; m<ntiles_out; m++) {
    nx2 = grid_out[m].nxc;
    ny2 = grid_out[m].nyc;
    out_area = (double *)malloc(nx2*ny2*sizeof(double));

    for(i=0; i<nx2*ny2; i++) {
      u_out[m].data[i] = 0.0;
      v_out[m].data[i] = 0.0;
    }
    for(i=0; i<nx2*ny2; i++) out_area[i] = 0.0;

    for(n=0; n<interp[m].nxgrid; n++) {
      i2   = interp[m].i_out[n];
      j2   = interp[m].j_out[n];
      i1   = interp[m].i_in [n];
      j1   = interp[m].j_in [n];
      tile = interp[m].t_in [n];
      area = interp[m].area [n];
      nx1  = grid_in[tile].nx;
      ny1  = grid_in[tile].ny;
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
    u_out[m].data[i] /= grid_out[m].area[i];
    v_out[m].data[i] /= grid_out[m].area[i];
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
    if(grid_out[m].rotate) {
      for(i=0; i<nx2*ny2; i++) {
  tmp_x = u_out[m].data[i];
  tmp_y = v_out[m].data[i];
  if( tmp_x != missing && tmp_y != missing) {
    u_out[m].data[i] =  tmp_x * grid_out[m].cosrot[i] + tmp_y * grid_out[m].sinrot[i];
    v_out[m].data[i] = -tmp_x * grid_out[m].sinrot[i] + tmp_y * grid_out[m].cosrot[i];
  }
      }
    }
    free(out_area);
  }

}; /* do_vector_conserve_interp */

void do_create_xgrid_order1( const int n, const int m,
                             const Grid_config *grid_in, const Grid_config *grid_out,
                             Interp_config *interp, unsigned int opcode )
{

  int    *i_in=NULL, *j_in=NULL, *i_out=NULL, *j_out=NULL ;
  double *xgrid_area=NULL, *xgrid_clon=NULL, *xgrid_clat=NULL ;
  double *mask;

  int mxxgrid, nxgrid;
  int nx_out, ny_out, nx_in, ny_in ;
  int jstart, jend, ny_now;

  int zero=0;

  nx_out = grid_out[n].nxc;
  ny_out = grid_out[n].nyc;
  nx_in = grid_in[m].nx;
  ny_in = grid_in[m].ny;

  mask = (double *)malloc(nx_in*ny_in*sizeof(double));
  for(int i=0; i<nx_in*ny_in; i++) mask[i] = 1.0;

  get_jstart_jend( nx_out, ny_out, nx_in, ny_in,
                   grid_out[n].latc, grid_in[m].latc, &jstart, &jend, &ny_now);

  mxxgrid=get_maxxgrid();
  malloc_xgrid_arrays(mxxgrid, &i_in, &j_in, &i_out, &j_out, &xgrid_area, &xgrid_clon, &xgrid_clat);
  nxgrid = create_xgrid_2dx2d_order1(&nx_in, &ny_now, &nx_out, &ny_out, grid_in[m].lonc+jstart*(nx_in+1),
                                     grid_in[m].latc+jstart*(nx_in+1),  grid_out[n].lonc,  grid_out[n].latc,
                                     mask, i_in, j_in, i_out, j_out, xgrid_area);
  for(int i=0; i<nxgrid; i++) j_in[i] += jstart;

  get_interp( opcode, nxgrid, interp, m, n, i_in, j_in, i_out, j_out,
              xgrid_clon, xgrid_clat, xgrid_area );
  malloc_xgrid_arrays(zero, &i_in, &j_in, &i_out, &j_out, &xgrid_area, &xgrid_clon , &xgrid_clat);
  free(mask);

}

void do_create_xgrid_order2( const int n, const int m, const Grid_config *grid_in, const Grid_config *grid_out,
                             Minmaxavg_lists *out_minmaxavg_lists, CellStruct *cell_in, Interp_config *interp,
                             unsigned int opcode )
{

  int    *i_in=NULL, *j_in=NULL, *i_out=NULL, *j_out=NULL ;
  double *xgrid_area=NULL, *xgrid_clon=NULL, *xgrid_clat=NULL ;
  double *mask;

  int nxgrid, approx_nxgrid;
  int nx_out, ny_out, nx_in, ny_in ;
  int jstart, jend, ny_now;

  int zero=0;
  clock_t time_start, time_end, time_nxgrid;

  int *counts_per_ij1=NULL, *ij2_start=NULL, *ij2_end=NULL;

  time_start = clock();

  nx_out = grid_out[n].nxc;
  ny_out = grid_out[n].nyc;
  nx_in = grid_in[m].nx;
  ny_in = grid_in[m].ny;

  get_jstart_jend( nx_out, ny_out, nx_in, ny_in,
                   grid_out[n].latc, grid_in[m].latc, &jstart, &jend, &ny_now);

  mask = (double *)malloc(nx_in*ny_in*sizeof(double));
  for(int i=0; i<nx_in*ny_in; i++) mask[i] = 1.0;

  counts_per_ij1 = (int *)malloc( nx_in*ny_in*sizeof(int) );
  ij2_start = (int *)malloc( nx_in*ny_in*sizeof(int) );
  ij2_end = (int *)malloc( nx_in*ny_in*sizeof(int) );

#pragma acc enter data copyin(grid_in[m].latc[0:(nx_in+1)*(ny_in+1)], \
                              grid_in[m].lonc[0:(nx_in+1)*(ny_in+1)], mask[0:nx_in*ny_in])
#pragma acc enter data create(counts_per_ij1[0:nx_in*ny_in], ij2_start[0:nx_in*ny_in], ij2_end[0:nx_in*ny_in])

#ifdef _OPENACC

  approx_nxgrid = prepare_create_xgrid_2dx2d_order2_acc(&nx_in, &ny_now, &nx_out, &ny_out, grid_in[m].lonc+jstart*(nx_in+1),
                                                 grid_in[m].latc+jstart*(nx_in+1), grid_out[n].lonc, grid_out[n].latc,
                                                 out_minmaxavg_lists, mask, counts_per_ij1, ij2_start, ij2_end );

  nxgrid = create_xgrid_2dx2d_order2_acc(&nx_in, &ny_now, &nx_out, &ny_out, grid_in[m].lonc+jstart*(nx_in+1),
                                         grid_in[m].latc+jstart*(nx_in+1), grid_out[n].lonc, grid_out[n].latc,
                                         out_minmaxavg_lists, mask, approx_nxgrid, counts_per_ij1, ij2_start, ij2_end,
                                         &i_in, &j_in, &i_out, &j_out, &xgrid_area, &xgrid_clon, &xgrid_clat);
#else

  nxgrid = create_xgrid_2dx2d_order2(&nx_in, &ny_now, &nx_out, &ny_out, grid_in[m].lonc+jstart*(nx_in+1),
                                     grid_in[m].latc+jstart*(nx_in+1),  grid_out[n].lonc,  grid_out[n].latc,
                                     mask, i_in, j_in, i_out, j_out, xgrid_area, xgrid_clon, xgrid_clat);

#endif
  if(DEBUG) printf("nxgrid, m, & n is: %d %d %d\n",nxgrid, m, n);
  time_end = clock();
  time_nxgrid += 1.0 * (time_end - time_start)/CLOCKS_PER_SEC;

  for(int i=0; i<nxgrid; i++) j_in[i] += jstart;

  /* For the purpose of bitiwise reproducing, the following operation is needed. */
  get_CellStruct(m, nx_in, nxgrid, i_in, j_in, xgrid_area, xgrid_clon, xgrid_clat, cell_in);
  get_interp( opcode, nxgrid, interp, m, n, i_in, j_in, i_out, j_out,
              xgrid_clon, xgrid_clat, xgrid_area );

  malloc_xgrid_arrays(zero, &i_in, &j_in, &i_out, &j_out, &xgrid_area, &xgrid_clon , &xgrid_clat);
  free(mask);
  free(counts_per_ij1);
  free(ij2_start);
  free(ij2_end);

#pragma acc exit data delete(grid_in[m].latc, grid_in[m].lonc)
#pragma acc exit data delete(mask)
#pragma acc exit data delete(counts_per_ij1, ij2_start, ij2_end)

}

void do_great_circle( const int n, const int m, const Grid_config *grid_in, const Grid_config *grid_out,
                      Interp_config *interp, unsigned int opcode)
{

  int    *i_in=NULL, *j_in=NULL, *i_out=NULL, *j_out=NULL ;
  double *xgrid_area=NULL, *xgrid_clon=NULL, *xgrid_clat=NULL ;
  double *mask;

  int nx_out, ny_out, nx_in, ny_in ;
  int mxxgrid, nxgrid;
  int zero=0;

  nx_out = grid_out[n].nxc;
  ny_out = grid_out[n].nyc;
  nx_in = grid_in[m].nx;
  ny_in = grid_in[m].ny;

  mask = (double *)malloc(nx_in*ny_in*sizeof(double));
  for(int i=0; i<nx_in*ny_in; i++) mask[i] = 1.0;
  mxxgrid=get_maxxgrid();
  malloc_xgrid_arrays(mxxgrid, &i_in, &j_in, &i_out, &j_out, &xgrid_area, &xgrid_clon, &xgrid_clat);
  nxgrid = create_xgrid_great_circle(&nx_in, &ny_in, &nx_out, &ny_out, grid_in[m].lonc,
                                     grid_in[m].latc,  grid_out[n].lonc,  grid_out[n].latc,
                                     mask, i_in, j_in, i_out, j_out, xgrid_area, xgrid_clon, xgrid_clat);
  get_interp( opcode, nxgrid, interp, m, n, i_in, j_in, i_out, j_out,
              xgrid_clon, xgrid_clat, xgrid_area );
  malloc_xgrid_arrays(zero, &i_in, &j_in, &i_out, &j_out, &xgrid_area, &xgrid_clon , &xgrid_clat);

}
