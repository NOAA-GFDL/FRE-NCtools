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
#include "globals.h"
#include "mosaic_util.h"
#include "create_xgrid.h"
#include "create_xgrid_util.h"
#include "create_xgrid_acc.h"
#include "constant.h"

#define AREA_RATIO_THRESH (1.e-6)
#define MASK_THRESH       (0.5)
#define EPSLN8            (1.e-8)
#define EPSLN30           (1.0e-30)
#define EPSLN10           (1.0e-10)
#define MAX_V 8

/*******************************************************************************
  prepare_create_xgrid_2dx2d_order2_acc
*******************************************************************************/
int prepare_create_xgrid_2dx2d_order2_acc(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
                                          const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                                          Minmaxavg_lists *out_minmaxavg_lists, const double *mask_in,
                                          int *counts_per_ij1, int *ij2_start, int *ij2_end)
{

#define MAX_V 8
  int nx1, nx2, ny1, ny2, nx1p, nx2p;
  size_t nxgrid;

  nx1 = *nlon_in;
  ny1 = *nlat_in;
  nx2 = *nlon_out;
  ny2 = *nlat_out;
  nx1p = nx1 + 1;
  nx2p = nx2 + 1;

  nxgrid = 0;

#pragma acc data present(lon_out[0:(nx2+1)*(ny2+1)], lat_out[0:(nx2+1)*(ny2+1)])
#pragma acc data present(lon_in[0:(nx1+1)*(ny1+1)], lat_in[0:(nx1+1)*(ny1+1)], mask_in[0:nx1*ny1])
#pragma acc data present(out_minmaxavg_lists->lon_list[0:MAX_V*nx2*ny2], out_minmaxavg_lists->lat_list[0:MAX_V*nx2*ny2])
#pragma acc data present(out_minmaxavg_lists->n_list[0:nx2*ny2], out_minmaxavg_lists->lon_avg[0:nx2*ny2])
#pragma acc data present(out_minmaxavg_lists->lat_min_list[0:nx2*ny2], out_minmaxavg_lists->lat_max_list[0:nx2*ny2])
#pragma acc data present(out_minmaxavg_lists->lon_min_list[0:nx2*ny2], out_minmaxavg_lists->lon_max_list[0:nx2*ny2])
#pragma acc data present(counts_per_ij1[0:nx1*ny1], ij2_start[0:nx1*ny1], ij2_end[0:nx1*ny1])
#pragma acc data copy(nxgrid)
#pragma acc kernels
{
#pragma acc loop independent reduction(+:nxgrid)
  for( int ij1=0 ; ij1 < nx1*ny1 ; ij1++) {

    int i1, j1;
    int icount=0 ;
    int ij2_max=0 , ij2_min=nx2*ny2+1;

    i1 = ij1%nx1;
    j1 = ij1/nx1;

    counts_per_ij1[ij1]=0;

    if( mask_in[ij1] > MASK_THRESH ) {

      int n0, n1, n2, n3, n1_in;
      double lat_in_min, lat_in_max, lon_in_min, lon_in_max, lon_in_avg;
      double x1_in[MV], y1_in[MV];

      n0 = j1*nx1p+i1;       n1 = j1*nx1p+i1+1;
      n2 = (j1+1)*nx1p+i1+1; n3 = (j1+1)*nx1p+i1;
      x1_in[0] = lon_in[n0]; y1_in[0] = lat_in[n0];
      x1_in[1] = lon_in[n1]; y1_in[1] = lat_in[n1];
      x1_in[2] = lon_in[n2]; y1_in[2] = lat_in[n2];
      x1_in[3] = lon_in[n3]; y1_in[3] = lat_in[n3];
      lat_in_min = minval_double(4, y1_in);
      lat_in_max = maxval_double(4, y1_in);
      n1_in = fix_lon(x1_in, y1_in, 4, M_PI);
      lon_in_min = minval_double(n1_in, x1_in);
      lon_in_max = maxval_double(n1_in, x1_in);
      lon_in_avg = avgval_double(n1_in, x1_in);

#pragma acc loop independent reduction(+:nxgrid) reduction(+:icount) reduction(min:ij2_min) reduction(max:ij2_max)
      for(int ij2=0; ij2<nx2*ny2; ij2++) {

        int i2, j2, l;
        double dx, lon_out_min, lon_out_max;
        double x2_in[MAX_V], y2_in[MAX_V],  x_out[MV], y_out[MV];;

        i2 = ij2%nx2;
        j2 = ij2/nx2;

        if(out_minmaxavg_lists->lat_min_list[ij2] >= lat_in_max || out_minmaxavg_lists->lat_max_list[ij2] <= lat_in_min ) continue;

        /* adjust x2_in according to lon_in_avg*/
        lon_out_min = out_minmaxavg_lists->lon_min_list[ij2];
        lon_out_max = out_minmaxavg_lists->lon_max_list[ij2];
        dx = out_minmaxavg_lists->lon_avg[ij2] - lon_in_avg;

        if(dx < -M_PI ) {
          lon_out_min += TPI;
          lon_out_max += TPI;
        }

        else if (dx >  M_PI) {
          lon_out_min -= TPI;
          lon_out_max -= TPI;
        }

        /* x2_in should in the same range as x1_in after lon_fix, so no need to consider cyclic condition */
        if(lon_out_min >= lon_in_max || lon_out_max <= lon_in_min ) continue;

        nxgrid++;
        icount++;
        ij2_min = min(ij2_min, ij2);
        ij2_max = max(ij2_max, ij2);

      } //ij2

      counts_per_ij1[ij1] = icount;
      ij2_start[ij1] = ij2_min ;
      ij2_end[ij1] = ij2_max;

    } // mask
  } //ij1
} //kernel

  return nxgrid;

}


/*******************************************************************************
  create_xgrid_2dx2d_order2 OPENACC version
*******************************************************************************/
int create_xgrid_2dx2d_order2_acc(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
                                  const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                                  Minmaxavg_lists *out_minmaxavg_lists, const double *mask_in, const int approx_nxgrid,
                                  const int *counts_per_ij1, const int *ij2_start, const int *ij2_end,
                                  int *i_in, int *j_in, int *i_out, int *j_out,
                                  double *xgrid_area, double *xgrid_clon, double *xgrid_clat)
{

#define MAX_V 8
  int nx1, nx2, ny1, ny2, nx1p, nx2p;
  double *area_in, *area_out;
  int mxxgrid;

  size_t nxgrid;

  nx1 = *nlon_in;
  ny1 = *nlat_in;
  nx2 = *nlon_out;
  ny2 = *nlat_out;
  nx1p = nx1 + 1;
  nx2p = nx2 + 1;
  mxxgrid = get_maxxgrid();

  area_in  = (double *)malloc(nx1*ny1*sizeof(double));
  area_out = (double *)malloc(nx2*ny2*sizeof(double));
  get_grid_area(nlon_in, nlat_in, lon_in, lat_in, area_in);
  get_grid_area(nlon_out, nlat_out, lon_out, lat_out, area_out);

  nxgrid = 0;

#pragma acc data present(lon_out[0:(nx2+1)*(ny2+1)], lat_out[0:(nx2+1)*(ny2+1)])
#pragma acc data present(lon_in[0:(nx1+1)*(ny1+1)], lat_in[0:(nx1+1)*(ny1+1)], mask_in[0:nx1*ny1])
#pragma acc data present(out_minmaxavg_lists->lon_list[0:MAX_V*nx2*ny2], out_minmaxavg_lists->lat_list[0:MAX_V*nx2*ny2])
#pragma acc data present(out_minmaxavg_lists->n_list[0:nx2*ny2], out_minmaxavg_lists->lon_avg[0:nx2*ny2])
#pragma acc data present(out_minmaxavg_lists->lat_min_list[0:nx2*ny2], out_minmaxavg_lists->lat_max_list[0:nx2*ny2])
#pragma acc data present(out_minmaxavg_lists->lon_min_list[0:nx2*ny2], out_minmaxavg_lists->lon_max_list[0:nx2*ny2])
#pragma acc data present(xgrid_area[0:mxxgrid], xgrid_clon[0:mxxgrid], xgrid_clat[0:mxxgrid], \
                         i_in[0:mxxgrid], j_in[0:mxxgrid], i_out[0:mxxgrid],j_out[0:mxxgrid])
#pragma acc data copyin(area_in[0:nx1*ny1], area_out[0:nx2*ny2])
#pragma acc data copy(nxgrid)
#pragma acc kernels
{
#pragma acc loop independent //reduction(+:nxgrid)
  for( int ij1=0 ; ij1<nx1*ny1 ; ij1++) {

    int i1, j1;

    i1 = ij1%nx1;
    j1 = ij1/nx1;

    if( mask_in[ij1] > MASK_THRESH ) {

      int n0, n1, n2, n3, n1_in;
      double lat_in_min,lat_in_max,lon_in_min,lon_in_max,lon_in_avg;
      double x1_in[MV], y1_in[MV];

      n0 = j1*nx1p+i1;       n1 = j1*nx1p+i1+1;
      n2 = (j1+1)*nx1p+i1+1; n3 = (j1+1)*nx1p+i1;
      x1_in[0] = lon_in[n0]; y1_in[0] = lat_in[n0];
      x1_in[1] = lon_in[n1]; y1_in[1] = lat_in[n1];
      x1_in[2] = lon_in[n2]; y1_in[2] = lat_in[n2];
      x1_in[3] = lon_in[n3]; y1_in[3] = lat_in[n3];
      lat_in_min = minval_double(4, y1_in);
      lat_in_max = maxval_double(4, y1_in);

      n1_in = fix_lon(x1_in, y1_in, 4, M_PI);
      lon_in_min = minval_double(n1_in, x1_in);
      lon_in_max = maxval_double(n1_in, x1_in);
      lon_in_avg = avgval_double(n1_in, x1_in);

#pragma acc loop independent //reduction(+:nxgrid)
      for(int ij2=0; ij2<nx2*ny2; ij2++) {
        int n_out, i2, j2, n2_in, l;
        double xarea, dx, lon_out_min, lon_out_max;
        double x2_in[MAX_V], y2_in[MAX_V],  x_out[MV], y_out[MV];;

        i2 = ij2%nx2;
        j2 = ij2/nx2;

        /* adjust x2_in according to lon_in_avg*/
        n2_in = out_minmaxavg_lists->n_list[ij2];
#pragma acc loop seq
        for(l=0; l<n2_in; l++) {
          x2_in[l] = out_minmaxavg_lists->lon_list[ij2*MAX_V+l];
          y2_in[l] = out_minmaxavg_lists->lat_list[ij2*MAX_V+l];
        }
        lon_out_min = out_minmaxavg_lists->lon_min_list[ij2];
        lon_out_max = out_minmaxavg_lists->lon_max_list[ij2];
        dx = out_minmaxavg_lists->lon_avg[ij2] - lon_in_avg;
        if(dx < -M_PI ) {
          lon_out_min += TPI;
          lon_out_max += TPI;
#pragma acc loop seq
          for (l=0; l<n2_in; l++) x2_in[l] += TPI;
        }
        else if (dx >  M_PI) {
          lon_out_min -= TPI;
          lon_out_max -= TPI;
#pragma acc loop seq
          for (l=0; l<n2_in; l++) x2_in[l] -= TPI;
        }

        n_out = 1;
        if (  (n_out = clip_2dx2d( x1_in, y1_in, n1_in, x2_in, y2_in, n2_in, x_out, y_out )) > 0) {
          double min_area;
          xarea = poly_area (x_out, y_out, n_out ) * mask_in[ij1];
          min_area = min(area_in[ij1], area_out[ij2]);
          printf("HEREEND3 %d %d\n", ij1, ij2);
          if( xarea/min_area > AREA_RATIO_THRESH ) {
            xgrid_area[nxgrid] = xarea;
            xgrid_clon[nxgrid] = poly_ctrlon(x_out, y_out, n_out, lon_in_avg);
            xgrid_clat[nxgrid] = poly_ctrlat (x_out, y_out, n_out );
            i_in[nxgrid]       = i1;
            j_in[nxgrid]       = j1;
            i_out[nxgrid]      = i2;
            j_out[nxgrid]      = j2;
#pragma atomic update
            nxgrid++;
          } //if
        } //if
      } //ij2
    } //mask
  } //ij1
} //kernel

  free(area_in);
  free(area_out);

  return nxgrid;

};/* get_xgrid_2Dx2D_order2 */
