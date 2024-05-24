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
#include <math.h>
#include "general_utils_acc.h"
#include "create_xgrid_acc.h"
#include "create_xgrid_utils_acc.h"
#include "globals_acc.h"

/*******************************************************************************
void get_nxgrid_upbound_2dx2d_acc
This function computes the upperbound to nxgrid.  This upper bound will be used
to malloc arrays used in create_xgrid
*******************************************************************************/
int get_nxgrid_upbound_2dx2d_acc(const int nlon_input_cells, const int nlat_input_cells,
                                 const int nlon_output_cells, const int nlat_output_cells,
                                 const double *intput_grid_lon, const double *input_grid_lat,
                                 const double *output_grid_lon, const double *output_grid_lat,
                                 const double *skip_input_cells, const Grid_cells_struct_config *output_grid_cell,
                                 int *xcells_per_ij1, int *ij2_start, int *ij2_end)
{

  int input_grid_ncells  = nlon_input_cells*nlat_input_cells;
  int output_grid_ncells = nlon_output_cells*nlat_output_cells;
  int input_grid_npts    = (nlon_input_cells+1)*(nlat_input_cells+1);
  int output_grid_npts   = (nlon_output_cells+1)*(nlat_output_cells+1);

  int nxgrid_upbound=0;

#pragma acc data present(output_grid_lon[:output_grid_npts], output_grid_lat[:output_grid_npts], \
                         intput_grid_lon[:input_grid_npts], input_grid_lat[:input_grid_npts],    \
                         output_grid_cell[:1], xcells_per_ij1[0:input_grid_ncells],              \
                         ij2_start[0:input_grid_ncells], ij2_end[0:input_grid_ncells],           \
                         skip_input_cells[:input_grid_ncells])
#pragma acc data copyin(input_grid_ncells, output_grid_ncells) copy(nxgrid_upbound)
#pragma acc parallel loop independent reduction(+:nxgrid_upbound)
  for( int ij1=0 ; ij1<input_grid_ncells ; ij1++) {
    if( skip_input_cells[ij1] > MASK_THRESH ) {
      int nvertices;
      int i_xcells_per_ij1=0;
      int ij2_max=0;
      int ij2_min=input_grid_ncells;
      double input_cell_lat_min, input_cell_lat_max, input_cell_lon_min, input_cell_lon_max, input_cell_lon_cent;
      double x1[MV], y1[MV];

      xcells_per_ij1[ij1]=0;

      get_cell_vertices_acc(ij1, nlon_input_cells, intput_grid_lon, input_grid_lat, x1, y1);

      input_cell_lat_min = minval_double_acc(4, y1);
      input_cell_lat_max = maxval_double_acc(4, y1);
      nvertices = fix_lon_acc(x1, y1, 4, M_PI);
      input_cell_lon_min = minval_double_acc(nvertices, x1);
      input_cell_lon_max = maxval_double_acc(nvertices, x1);
      input_cell_lon_cent = avgval_double_acc(nvertices, x1);

#pragma acc loop independent reduction(+:nxgrid_upbound) reduction(+:i_xcells_per_ij1) \
  reduction(min:ij2_min) reduction(max:ij2_max)
      for(int ij2=0; ij2<output_grid_ncells; ij2++) {

        double dlon_cent, output_cell_lon_min, output_cell_lon_max;

        if(output_grid_cell->lat_min[ij2] >= input_cell_lat_max) continue;
        if(output_grid_cell->lat_max[ij2] <= input_cell_lat_min) continue;

        /* adjust according to intput_grid_lon_cent*/
        output_cell_lon_min = output_grid_cell->lon_min[ij2];
        output_cell_lon_max = output_grid_cell->lon_max[ij2];

        dlon_cent = output_grid_cell->lon_cent[ij2] - input_cell_lon_cent;
        if(dlon_cent < -M_PI ) {
          output_cell_lon_min += TPI;
          output_cell_lon_max += TPI;
        }
        else if (dlon_cent >  M_PI) {
          output_cell_lon_min -= TPI;
          output_cell_lon_max -= TPI;
        }

        /* x2_in should in the same range as x1_in after lon_fix, so no need to consider cyclic condition */
        if(output_cell_lon_min >= input_cell_lon_max ) continue;
        if(output_cell_lon_max <= input_cell_lon_min ) continue;

        //Note, the check for AREA_RATIO_THRESH has been removed
        //Thus, the computed value of nxgrid_upbound will be equal to or greater than nxgrid
        i_xcells_per_ij1++;
        nxgrid_upbound++;
        ij2_min = min(ij2_min, ij2);
        ij2_max = max(ij2_max, ij2);

      } //ij2
      xcells_per_ij1[ij1] = i_xcells_per_ij1;
      ij2_start[ij1] = ij2_min ;
      ij2_end[ij1]   = ij2_max;

    } //mask
  } //ij1

  return nxgrid_upbound;

}

/*******************************************************************************
  void create_xgrid_2DX2D_order1
  This routine generate exchange grids between two grids for the first order
  conservative interpolation. nlon_input_cells,ninput_grid_lat,nlon_output_cells,nlat_output_cells are the size of the grid cell
  and input_grid_lon,input_grid_lat, output_grid_lon,output_grid_lat are geographic grid location of grid cell bounds.
*******************************************************************************/
int create_xgrid_2dx2d_order1_acc(const int *nlon_input_cells, const int *nlat_input_cells,
                                  const int *nlon_output_cells, const int *nlat_output_cells,
                                  const double *input_grid_lon, const double *input_grid_lat,
                                  const double *output_grid_lon, const double *output_grid_lat,
                                  const double *skip_input_cells,
                                  int *i_in, int *j_in, int *i_out, int *j_out, double *xgrid_area)
{

  int nx1, nx2, ny1, ny2, nx1p, nx2p, nxgrid;
  double *area_in, *area_out;
  int  ij;
  double *output_grid_lon_min_list,*output_grid_lon_max_list,*output_grid_lon_avg;
  double *output_grid_lat_min_list,*output_grid_lat_max_list;
  double *output_grid_lon_list, *output_grid_lat_list;
  int    *n2_list;

  nx1 = *nlon_input_cells;
  ny1 = *nlat_input_cells;
  nx2 = *nlon_output_cells;
  ny2 = *nlat_output_cells;
  nx1p = nx1 + 1;
  nx2p = nx2 + 1;

  area_in  = (double *)malloc(nx1*ny1*sizeof(double));
  area_out = (double *)malloc(nx2*ny2*sizeof(double));
  get_grid_area_acc(nlon_input_cells, nlat_input_cells, input_grid_lon, input_grid_lat, area_in);
  get_grid_area_acc(nlon_output_cells, nlat_output_cells, output_grid_lon, output_grid_lat, area_out);

  output_grid_lon_min_list = (double *)malloc(nx2*ny2*sizeof(double));
  output_grid_lon_max_list = (double *)malloc(nx2*ny2*sizeof(double));
  output_grid_lat_min_list = (double *)malloc(nx2*ny2*sizeof(double));
  output_grid_lat_max_list = (double *)malloc(nx2*ny2*sizeof(double));
  output_grid_lon_avg = (double *)malloc(nx2*ny2*sizeof(double));
  n2_list     = (int *)malloc(nx2*ny2*sizeof(int));
  output_grid_lon_list = (double *)malloc(MAX_V*nx2*ny2*sizeof(double));
  output_grid_lat_list = (double *)malloc(MAX_V*nx2*ny2*sizeof(double));

  for(ij=0; ij<nx2*ny2; ij++){
    int i2, j2, n, n0, n1, n2, n3, n2_in, l;
    double x2_in[MV], y2_in[MV];
    i2 = ij%nx2;
    j2 = ij/nx2;
    n = j2*nx2+i2;
    n0 = j2*nx2p+i2; n1 = j2*nx2p+i2+1;
    n2 = (j2+1)*nx2p+i2+1; n3 = (j2+1)*nx2p+i2;
    x2_in[0] = output_grid_lon[n0]; y2_in[0] = output_grid_lat[n0];
    x2_in[1] = output_grid_lon[n1]; y2_in[1] = output_grid_lat[n1];
    x2_in[2] = output_grid_lon[n2]; y2_in[2] = output_grid_lat[n2];
    x2_in[3] = output_grid_lon[n3]; y2_in[3] = output_grid_lat[n3];

    output_grid_lat_min_list[n] = minval_double_acc(4, y2_in);
    output_grid_lat_max_list[n] = maxval_double_acc(4, y2_in);
    n2_in = fix_lon_acc(x2_in, y2_in, 4, M_PI);
    output_grid_lon_min_list[n] = minval_double_acc(n2_in, x2_in);
    output_grid_lon_max_list[n] = maxval_double_acc(n2_in, x2_in);
    output_grid_lon_avg[n] = avgval_double_acc(n2_in, x2_in);
    n2_list[n] = n2_in;
    for(l=0; l<n2_in; l++) {
      output_grid_lon_list[n*MAX_V+l] = x2_in[l];
      output_grid_lat_list[n*MAX_V+l] = y2_in[l];
    }
  }

  nxgrid = 0;

  for(int j1=0; j1<ny1; j1++) for(int i1=0; i1<nx1; i1++) {
      if(skip_input_cells[j1*nx1+i1] > MASK_THRESH)  {
        int n0, n1, n2, n3, l,n1_in;
        double input_grid_lat_min,input_grid_lat_max,input_grid_lon_min,input_grid_lon_max,input_grid_lon_avg;
        double x1_in[MV], y1_in[MV], x_out[MV], y_out[MV];

        n0 = j1*nx1p+i1;       n1 = j1*nx1p+i1+1;
        n2 = (j1+1)*nx1p+i1+1; n3 = (j1+1)*nx1p+i1;
        x1_in[0] = input_grid_lon[n0]; y1_in[0] = input_grid_lat[n0];
        x1_in[1] = input_grid_lon[n1]; y1_in[1] = input_grid_lat[n1];
        x1_in[2] = input_grid_lon[n2]; y1_in[2] = input_grid_lat[n2];
        x1_in[3] = input_grid_lon[n3]; y1_in[3] = input_grid_lat[n3];
        input_grid_lat_min = minval_double_acc(4, y1_in);
        input_grid_lat_max = maxval_double_acc(4, y1_in);
        n1_in = fix_lon_acc(x1_in, y1_in, 4, M_PI);
        input_grid_lon_min = minval_double_acc(n1_in, x1_in);
        input_grid_lon_max = maxval_double_acc(n1_in, x1_in);
        input_grid_lon_avg = avgval_double_acc(n1_in, x1_in);
        for(ij=0; ij<nx2*ny2; ij++) {
          int n_in, n_out, i2, j2, n2_in;
          double xarea, dx, output_grid_lon_min, output_grid_lon_max;
          double x2_in[MAX_V], y2_in[MAX_V];

          i2 = ij%nx2;
          j2 = ij/nx2;

          if(output_grid_lat_min_list[ij] >= input_grid_lat_max || output_grid_lat_max_list[ij] <= input_grid_lat_min ) continue;
          /* adjust x2_in according to input_grid_lon_avg*/
          n2_in = n2_list[ij];
          for(l=0; l<n2_in; l++) {
            x2_in[l] = output_grid_lon_list[ij*MAX_V+l];
            y2_in[l] = output_grid_lat_list[ij*MAX_V+l];
          }
          output_grid_lon_min = output_grid_lon_min_list[ij];
          output_grid_lon_max = output_grid_lon_max_list[ij];
          dx = output_grid_lon_avg[ij] - input_grid_lon_avg;
          if(dx < -M_PI ) {
            output_grid_lon_min += TPI;
            output_grid_lon_max += TPI;
            for (l=0; l<n2_in; l++) x2_in[l] += TPI;
          }
          else if (dx >  M_PI) {
            output_grid_lon_min -= TPI;
            output_grid_lon_max -= TPI;
            for (l=0; l<n2_in; l++) x2_in[l] -= TPI;
          }

          /* x2_in should in the same range as x1_in after lon_fix, so no need to
             consider cyclic condition
          */
          if(output_grid_lon_min >= input_grid_lon_max || output_grid_lon_max <= input_grid_lon_min ) continue;
          if (  (n_out = clip_2dx2d_acc( x1_in, y1_in, n1_in, x2_in, y2_in, n2_in, x_out, y_out )) > 0) {
            double min_area;
            int    nn;
            xarea = poly_area_acc(x_out, y_out, n_out );
            min_area = min(area_in[j1*nx1+i1], area_out[j2*nx2+i2]);
            if( xarea/min_area > AREA_RATIO_THRESH ) {
              xgrid_area[nxgrid] = xarea;
              i_in[nxgrid]       = i1;
              j_in[nxgrid]       = j1;
              i_out[nxgrid]      = i2;
              j_out[nxgrid]      = j2;
              nxgrid++;
            }

          }

        }
      }
  }

  free(area_in);
  free(area_out);
  free(output_grid_lon_min_list);
  free(output_grid_lon_max_list);
  free(output_grid_lat_min_list);
  free(output_grid_lat_max_list);
  free(output_grid_lon_avg);
  free(n2_list);
  free(output_grid_lon_list);
  free(output_grid_lat_list);

  return nxgrid;

};/* get_xgrid_2Dx2D_order1 */

/********************************************************************************
  void create_xgrid_2dx1d_order2
  This routine generate exchange grids between two grids for the second order
  conservative interpolation. nlon_input_cells,nlat_input_cells,nlon_output_cells,nlat_output_cells are the size of the grid cell
  and input_grid_lon,input_grid_lat, output_grid_lon,output_grid_lat are geographic grid location of grid cell bounds.
********************************************************************************/
int create_xgrid_2dx2d_order2_acc(const int *nlon_input_cells, const int *nlat_input_cells,
                                  const int *nlon_output_cells, const int *nlat_output_cells,
                                  const double *input_grid_lon, const double *input_grid_lat,
                                  const double *output_grid_lon, const double *output_grid_lat,
                                  const double *skip_input_cells,
                                  int *i_in, int *j_in, int *i_out, int *j_out, double *xgrid_area,
                                  double *xgrid_clon, double *xgrid_clat)
{

  int nx1, nx2, ny1, ny2, nx1p, nx2p, nxgrid;
  double *area_in, *area_out;
  int ij;
  double *output_grid_lon_min_list,*output_grid_lon_max_list;
  double *output_grid_lon_avg,*output_grid_lat_min_list,*output_grid_lat_max_list;
  double *output_grid_lon_list, *output_grid_lat_list;
  int    *n2_list;

  nx1 = *nlon_input_cells;
  ny1 = *nlat_input_cells;
  nx2 = *nlon_output_cells;
  ny2 = *nlat_output_cells;
  nx1p = nx1 + 1;
  nx2p = nx2 + 1;

  area_in  = (double *)malloc(nx1*ny1*sizeof(double));
  area_out = (double *)malloc(nx2*ny2*sizeof(double));
  get_grid_area_acc(nlon_input_cells, nlat_input_cells, input_grid_lon, input_grid_lat, area_in);
  get_grid_area_acc(nlon_output_cells, nlat_output_cells, output_grid_lon, output_grid_lat, area_out);

  output_grid_lon_min_list = (double *)malloc(nx2*ny2*sizeof(double));
  output_grid_lon_max_list = (double *)malloc(nx2*ny2*sizeof(double));
  output_grid_lat_min_list = (double *)malloc(nx2*ny2*sizeof(double));
  output_grid_lat_max_list = (double *)malloc(nx2*ny2*sizeof(double));
  output_grid_lon_avg = (double *)malloc(nx2*ny2*sizeof(double));
  n2_list     = (int *)malloc(nx2*ny2*sizeof(int));
  output_grid_lon_list = (double *)malloc(MAX_V*nx2*ny2*sizeof(double));
  output_grid_lat_list = (double *)malloc(MAX_V*nx2*ny2*sizeof(double));

  for(ij=0; ij<nx2*ny2; ij++){
    int i2, j2, n, n0, n1, n2, n3, n2_in, l;
    double x2_in[MV], y2_in[MV];
    i2 = ij%nx2;
    j2 = ij/nx2;
    n = j2*nx2+i2;
    n0 = j2*nx2p+i2; n1 = j2*nx2p+i2+1;
    n2 = (j2+1)*nx2p+i2+1; n3 = (j2+1)*nx2p+i2;
    x2_in[0] = output_grid_lon[n0]; y2_in[0] = output_grid_lat[n0];
    x2_in[1] = output_grid_lon[n1]; y2_in[1] = output_grid_lat[n1];
    x2_in[2] = output_grid_lon[n2]; y2_in[2] = output_grid_lat[n2];
    x2_in[3] = output_grid_lon[n3]; y2_in[3] = output_grid_lat[n3];

    output_grid_lat_min_list[n] = minval_double_acc(4, y2_in);
    output_grid_lat_max_list[n] = maxval_double_acc(4, y2_in);
    n2_in = fix_lon_acc(x2_in, y2_in, 4, M_PI);
    output_grid_lon_min_list[n] = minval_double_acc(n2_in, x2_in);
    output_grid_lon_max_list[n] = maxval_double_acc(n2_in, x2_in);
    output_grid_lon_avg[n] = avgval_double_acc(n2_in, x2_in);
    n2_list[n] = n2_in;
    for(l=0; l<n2_in; l++) {
      output_grid_lon_list[n*MAX_V+l] = x2_in[l];
      output_grid_lat_list[n*MAX_V+l] = y2_in[l];
    }
  }

  nxgrid = 0;

  for(int j1=0; j1<ny1; j1++) for(int i1=0; i1<nx1; i1++) {
      if( skip_input_cells[j1*nx1+i1] > MASK_THRESH ) {
        int n0, n1, n2, n3, l,n1_in;
        double input_grid_lat_min,input_grid_lat_max,input_grid_lon_min,input_grid_lon_max,input_grid_lon_avg;
        double x1_in[MV], y1_in[MV], x_out[MV], y_out[MV];

        n0 = j1*nx1p+i1;       n1 = j1*nx1p+i1+1;
        n2 = (j1+1)*nx1p+i1+1; n3 = (j1+1)*nx1p+i1;
        x1_in[0] = input_grid_lon[n0]; y1_in[0] = input_grid_lat[n0];
        x1_in[1] = input_grid_lon[n1]; y1_in[1] = input_grid_lat[n1];
        x1_in[2] = input_grid_lon[n2]; y1_in[2] = input_grid_lat[n2];
        x1_in[3] = input_grid_lon[n3]; y1_in[3] = input_grid_lat[n3];
        input_grid_lat_min = minval_double_acc(4, y1_in);
        input_grid_lat_max = maxval_double_acc(4, y1_in);
        n1_in = fix_lon_acc(x1_in, y1_in, 4, M_PI);
        input_grid_lon_min = minval_double_acc(n1_in, x1_in);
        input_grid_lon_max = maxval_double_acc(n1_in, x1_in);
        input_grid_lon_avg = avgval_double_acc(n1_in, x1_in);
        for(ij=0; ij<nx2*ny2; ij++) {
          int n_in, n_out, i2, j2, n2_in;
          double xarea, dx, output_grid_lon_min, output_grid_lon_max;
          double x2_in[MAX_V], y2_in[MAX_V];

          i2 = ij%nx2;
          j2 = ij/nx2;

          if(output_grid_lat_min_list[ij] >= input_grid_lat_max || output_grid_lat_max_list[ij] <= input_grid_lat_min ) continue;
          /* adjust x2_in according to input_grid_lon_avg*/
          n2_in = n2_list[ij];
          for(l=0; l<n2_in; l++) {
            x2_in[l] = output_grid_lon_list[ij*MAX_V+l];
            y2_in[l] = output_grid_lat_list[ij*MAX_V+l];
          }
          output_grid_lon_min = output_grid_lon_min_list[ij];
          output_grid_lon_max = output_grid_lon_max_list[ij];
          dx = output_grid_lon_avg[ij] - input_grid_lon_avg;
          if(dx < -M_PI ) {
            output_grid_lon_min += TPI;
            output_grid_lon_max += TPI;
            for (l=0; l<n2_in; l++) x2_in[l] += TPI;
          }
          else if (dx >  M_PI) {
            output_grid_lon_min -= TPI;
            output_grid_lon_max -= TPI;
            for (l=0; l<n2_in; l++) x2_in[l] -= TPI;
          }

          /* x2_in should in the same range as x1_in after lon_fix, so no need to
             consider cyclic condition
          */
          if(output_grid_lon_min >= input_grid_lon_max || output_grid_lon_max <= input_grid_lon_min ) continue;
          if (  (n_out = clip_2dx2d_acc( x1_in, y1_in, n1_in, x2_in, y2_in, n2_in, x_out, y_out )) > 0) {
            double min_area;
            int nn;
            xarea = poly_area_acc(x_out, y_out, n_out ) ;
            min_area = min(area_in[j1*nx1+i1], area_out[j2*nx2+i2]);
            if( xarea/min_area > AREA_RATIO_THRESH ) {
              xgrid_area[nxgrid] = xarea;
              xgrid_clon[nxgrid] = poly_ctrlon_acc(x_out, y_out, n_out, input_grid_lon_avg);
              xgrid_clat[nxgrid] = poly_ctrlat_acc(x_out, y_out, n_out );
              i_in[nxgrid]       = i1;
              j_in[nxgrid]       = j1;
              i_out[nxgrid]      = i2;
              j_out[nxgrid]      = j2;
              nxgrid++;
            }
          }
        }
      }
  }

  free(area_in);
  free(area_out);
  free(output_grid_lon_min_list);
  free(output_grid_lon_max_list);
  free(output_grid_lat_min_list);
  free(output_grid_lat_max_list);
  free(output_grid_lon_avg);
  free(n2_list);
  free(output_grid_lon_list);
  free(output_grid_lat_list);

  return nxgrid;

};/* get_xgrid_2Dx2D_order2 */


int create_xgrid_great_circle_acc(const int *nlon_input_cells, const int *nlat_input_cells,
                                  const int *nlon_output_cells, const int *nlat_output_cells,
                                  const double *input_grid_lon, const double *input_grid_lat,
                                  const double *output_grid_lon, const double *output_grid_lat,
                                  const double *skip_input_cells,
                                  int *i_in, int *j_in, int *i_out, int *j_out,
                                  double *xgrid_area, double *xgrid_clon, double *xgrid_clat)
{

  int nx1, nx2, ny1, ny2, nx1p, nx2p, ny1p, ny2p, nxgrid, n1_in, n2_in;
  int n0, n1, n2, n3, i1, j1, i2, j2, l, n;
  double x1_in[MV], y1_in[MV], z1_in[MV];
  double x2_in[MV], y2_in[MV], z2_in[MV];
  double x_out[MV], y_out[MV], z_out[MV];
  double *x1=NULL, *y1=NULL, *z1=NULL;
  double *x2=NULL, *y2=NULL, *z2=NULL;

  double xctrlon, xctrlat;
  double *area1, *area2, min_area;

  nx1 = *nlon_input_cells;
  ny1 = *nlat_input_cells;
  nx2 = *nlon_output_cells;
  ny2 = *nlat_output_cells;
  nxgrid = 0;
  nx1p = nx1 + 1;
  nx2p = nx2 + 1;
  ny1p = ny1 + 1;
  ny2p = ny2 + 1;

  /* first convert lon-lat to cartesian coordinates */
  x1 = (double *)malloc(nx1p*ny1p*sizeof(double));
  y1 = (double *)malloc(nx1p*ny1p*sizeof(double));
  z1 = (double *)malloc(nx1p*ny1p*sizeof(double));
  x2 = (double *)malloc(nx2p*ny2p*sizeof(double));
  y2 = (double *)malloc(nx2p*ny2p*sizeof(double));
  z2 = (double *)malloc(nx2p*ny2p*sizeof(double));

  latlon2xyz_acc(nx1p*ny1p, input_grid_lon, input_grid_lat, x1, y1, z1);
  latlon2xyz_acc(nx2p*ny2p, output_grid_lon, output_grid_lat, x2, y2, z2);

  area1  = (double *)malloc(nx1*ny1*sizeof(double));
  area2 = (double *)malloc(nx2*ny2*sizeof(double));
  get_grid_great_circle_area_acc(nlon_input_cells, nlat_input_cells, input_grid_lon, input_grid_lat, area1);
  get_grid_great_circle_area_acc(nlon_output_cells, nlat_output_cells, output_grid_lon, output_grid_lat, area2);
  n1_in = 4;
  n2_in = 4;

  for(int j1=0; j1<ny1; j1++) for(int i1=0; i1<nx1; i1++) {
      if( skip_input_cells[j1*nx1+i1] > MASK_THRESH ) {
          /* clockwise */
          n0 = j1*nx1p+i1;       n1 = (j1+1)*nx1p+i1;
          n2 = (j1+1)*nx1p+i1+1; n3 = j1*nx1p+i1+1;
          x1_in[0] = x1[n0]; y1_in[0] = y1[n0]; z1_in[0] = z1[n0];
          x1_in[1] = x1[n1]; y1_in[1] = y1[n1]; z1_in[1] = z1[n1];
          x1_in[2] = x1[n2]; y1_in[2] = y1[n2]; z1_in[2] = z1[n2];
          x1_in[3] = x1[n3]; y1_in[3] = y1[n3]; z1_in[3] = z1[n3];

          for(j2=0; j2<ny2; j2++) for(i2=0; i2<nx2; i2++) {
              int n_in, n_out;
              double xarea;

              n0 = j2*nx2p+i2;       n1 = (j2+1)*nx2p+i2;
              n2 = (j2+1)*nx2p+i2+1; n3 = j2*nx2p+i2+1;
              x2_in[0] = x2[n0]; y2_in[0] = y2[n0]; z2_in[0] = z2[n0];
              x2_in[1] = x2[n1]; y2_in[1] = y2[n1]; z2_in[1] = z2[n1];
              x2_in[2] = x2[n2]; y2_in[2] = y2[n2]; z2_in[2] = z2[n2];
              x2_in[3] = x2[n3]; y2_in[3] = y2[n3]; z2_in[3] = z2[n3];

              if (  (n_out = clip_2dx2d_great_circle_acc( x1_in, y1_in, z1_in, n1_in, x2_in, y2_in, z2_in, n2_in,
                                                          x_out, y_out, z_out)) > 0) {
                xarea = great_circle_area_acc( n_out, x_out, y_out, z_out ) ;
                min_area = min(area1[j1*nx1+i1], area2[j2*nx2+i2]);
                if( xarea/min_area > AREA_RATIO_THRESH ) {
#ifdef debug_test_create_xgrid
                  printf("(i2,j2)=(%d,%d), (i1,j1)=(%d,%d), xarea=%g\n", i2, j2, i1, j1, xarea);
#endif
                  xgrid_area[nxgrid] = xarea;
                  xgrid_clon[nxgrid] = 0; /*z1l: will be developed very soon */
                  xgrid_clat[nxgrid] = 0;
                  i_in[nxgrid]       = i1;
                  j_in[nxgrid]       = j1;
                  i_out[nxgrid]      = i2;
                  j_out[nxgrid]      = j2;
                  ++nxgrid;
            }
          }
        }
     }
  }


  free(area1);
  free(area2);

  free(x1);
  free(y1);
  free(z1);
  free(x2);
  free(y2);
  free(z2);

  return nxgrid;

};/* create_xgrid_great_circle */
