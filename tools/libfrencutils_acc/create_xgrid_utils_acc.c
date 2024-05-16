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
                                 Cell *cell )
{

  int ncell=nlon*nlat;
  int npts=(nlon+1)*(nlat+1);

  cell->lon_min  = (double *)malloc(ncell*sizeof(double));
  cell->lon_max  = (double *)malloc(ncell*sizeof(double));
  cell->lat_min  = (double *)malloc(ncell*sizeof(double));
  cell->lat_max  = (double *)malloc(ncell*sizeof(double));
  cell->lon_cent = (double *)malloc(ncell*sizeof(double));
  cell->n_vertices = (int *)malloc(ncell*sizeof(int));
  cell->lon_vertices  = (double **)malloc(ncell*sizeof(double));
  cell->lat_vertices  = (double **)malloc(ncell*sizeof(double));

  for(int icell=0 ; icell<ncell ; icell++) {
    cell->lat_vertices[icell] = (double *)malloc(MAX_V*sizeof(double));
    cell->lon_vertices[icell] = (double *)malloc(MAX_V*sizeof(double));
  }

#pragma acc enter data copyin(cell[:1])
#pragma acc enter data copyin(cell->lon_min[:ncell], cell->lon_max[:ncell], \
                              cell->lat_min[:ncell], cell->lat_max[:ncell], \
                              cell->lon_cent[:ncell], cell->n_vertices[:ncell])
#pragma acc enter data copyin(cell->lon_vertices[:ncell][:MAX_V], \
                              cell->lat_vertices[:ncell][:MAX_V])

#pragma acc data present(cell[:1], lon[:npts], lat[:npts])
#pragma acc parallel loop independent
  for(int icell=0; icell<ncell; icell++){
    int n;
    double x[MV], y[MV];

    get_cell_vertices( icell, nlon, lon, lat, x, y );
    cell->lat_min[icell] = minval_double(4, y);
    cell->lat_max[icell] = maxval_double(4, y);

    n = fix_lon(x, y, 4, M_PI);
    cell->n_vertices[icell] = n;

    if(n > MAX_V) printf("ERROR get_cell_minmaxavg_latlons:  number of cell vertices is greater than MAX_V\n");
    cell->lon_min[icell]  = minval_double(n, x);
    cell->lon_max[icell]  = maxval_double(n, x);
    cell->lon_cent[icell] = avgval_double(n, x);

    for(int ivertex=0 ; ivertex<n ; ivertex++) {
      cell->lon_vertices[icell][ivertex] = x[ivertex];
      cell->lat_vertices[icell][ivertex] = y[ivertex];
    }
  }

}

void get_cell_vertices( const int icell, const int nlon, const double *lon, const double *lat, double *x, double *y )
{

  int i, j;
  int n0, n1, n2, n3;

  i = icell%nlon;
  j = icell/nlon;
  n0 = j*(nlon+1)+i;
  n1 = j*(nlon+1)+i+1;
  n2 = (j+1)*(nlon+1)+i+1;
  n3 = (j+1)*(nlon+1)+i;

  x[0] = lon[n0]; y[0] = lat[n0];
  x[1] = lon[n1]; y[1] = lat[n1];
  x[2] = lon[n2]; y[2] = lat[n2];
  x[3] = lon[n3]; y[3] = lat[n3];

}
