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
#include "mosaic_util.h"
#include "create_xgrid_utils_acc.h"
#include "globals.h"
#include "parameters.h"

/*******************************************************************************
void get_grid_area(const int *nlon, const int *nlat, const double *lon, const double *lat, const double *area)
  return the grid area.
*******************************************************************************/
void get_grid_area_acc(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area)
{
  int nx, ny, nxp, i, j, n_in;
  double x_in[20], y_in[20];

  nx = *nlon;
  ny = *nlat;
  nxp = nx + 1;

  for(j=0; j<ny; j++) for(i=0; i < nx; i++) {
      x_in[0] = lon[j*nxp+i];
      x_in[1] = lon[j*nxp+i+1];
      x_in[2] = lon[(j+1)*nxp+i+1];
      x_in[3] = lon[(j+1)*nxp+i];
      y_in[0] = lat[j*nxp+i];
      y_in[1] = lat[j*nxp+i+1];
      y_in[2] = lat[(j+1)*nxp+i+1];
      y_in[3] = lat[(j+1)*nxp+i];
      n_in = fix_lon(x_in, y_in, 4, M_PI);
      area[j*nx+i] = poly_area(x_in, y_in, n_in);
    }

};  /* get_grid_area */

void get_grid_great_circle_area(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area)
{
  int nx, ny, nxp, nyp, i, j, n_in;
  int n0, n1, n2, n3;
  double x_in[20], y_in[20], z_in[20];
  struct Node *grid=NULL;
  double *x=NULL, *y=NULL, *z=NULL;


  nx = *nlon;
  ny = *nlat;
  nxp = nx + 1;
  nyp = ny + 1;

  x = (double *)malloc(nxp*nyp*sizeof(double));
  y = (double *)malloc(nxp*nyp*sizeof(double));
  z = (double *)malloc(nxp*nyp*sizeof(double));

  latlon2xyz(nxp*nyp, lon, lat, x, y, z);

  for(j=0; j<ny; j++) for(i=0; i < nx; i++) {
      /* clockwise */
      n0 = j*nxp+i;
      n1 = (j+1)*nxp+i;
      n2 = (j+1)*nxp+i+1;
      n3 = j*nxp+i+1;
      rewindList();
      grid = getNext();
      addEnd(grid, x[n0], y[n0], z[n0], 0, 0, 0, -1);
      addEnd(grid, x[n1], y[n1], z[n1], 0, 0, 0, -1);
      addEnd(grid, x[n2], y[n2], z[n2], 0, 0, 0, -1);
      addEnd(grid, x[n3], y[n3], z[n3], 0, 0, 0, -1);
      area[j*nx+i] = gridArea(grid);
    }

  free(x);
  free(y);
  free(z);

};  /* get_grid_great_circle_area */


void get_cell_minmaxavg_latlons( const int nlon, const int nlat, const double *lon, const double *lat,
                                 Minmaxavg_list *minmaxavg_list )
{

  minmaxavg_list->lon_min = (double *)malloc(nlon*nlat*sizeof(double));
  minmaxavg_list->lon_max = (double *)malloc(nlon*nlat*sizeof(double));
  minmaxavg_list->lat_min = (double *)malloc(nlon*nlat*sizeof(double));
  minmaxavg_list->lat_max = (double *)malloc(nlon*nlat*sizeof(double));
  minmaxavg_list->lon_avg = (double *)malloc(nlon*nlat*sizeof(double));
  minmaxavg_list->n_vertices   = (int *)malloc(nlon*nlat*sizeof(int));
  minmaxavg_list->lon_vertices = (double *)malloc(MAX_V*nlon*nlat*sizeof(double));
  minmaxavg_list->lat_vertices = (double *)malloc(MAX_V*nlon*nlat*sizeof(double));

  for(int ij=0; ij<nlon*nlat; ij++){
    int n;
    double x[MV], y[MV];

    get_cell_verticies( ij, nlon, lon, lat, x, y );

    minmaxavg_list->lat_min[ij] = minval_double(4, y);
    minmaxavg_list->lat_max[ij] = maxval_double(4, y);

    n = fix_lon(x, y, 4, M_PI);
    minmaxavg_list->n_vertices[ij] = n;

    //if(n > MAX_V) error_handler("get_cell_minmaxavg_latlons: number of cell vertices is greater than MAX_V");
    minmaxavg_list->lon_min[ij] = minval_double(n, x);
    minmaxavg_list->lon_max[ij] = maxval_double(n, x);
    minmaxavg_list->lon_avg[ij] = avgval_double(n, x);

    for(int l=0; l<n; l++) {
      minmaxavg_list->lon_vertices[ij*MAX_V+l] = x[l];
      minmaxavg_list->lat_vertices[ij*MAX_V+l] = y[l];
    }
  }

}

void get_cell_vertices( const int ij, const nlon, const double *lon, const double *lat, double *x, double *y )
{

  int i, j;
  int n0, n1, n2, n3;

  i = ij%nlon;
  j = ij/nlon;
  n0 = j*(nlon+1)+i;
  n1 = j*(nlon+1)+i+1;
  n2 = (j+1)*(nlon+1)+i+1;
  n3 = (j+1)*(nlon+1)+i;

  x[0] = lon[n0]; y[0] = lat[n0];
  x[1] = lon[n1]; y[1] = lat[n1];
  x[2] = lon[n2]; y[2] = lat[n2];
  x[3] = lon[n3]; y[3] = lat[n3];

}
