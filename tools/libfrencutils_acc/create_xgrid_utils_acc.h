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
#ifndef CREATE_XGRID_UTILS_ACC_H_
#define CREATE_XGRID_UTILS_ACC_H_

#include "globals.h"

#define MV 50
/* this value is small compare to earth area */

void get_grid_area_acc(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area);

void get_grid_great_circle_area_acc(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area);

double poly_ctrlon_acc(const double x[], const double y[], int n, double clon);

double poly_ctrlat_acc(const double x[], const double y[], int n);

int clip_2dx2d_acc(const double lon1_in[], const double lat1_in[], int n1_in,
                   const double lon2_in[], const double lat2_in[], int n2_in,
                   double lon_out[], double lat_out[]);

int clip_2dx2d_great_circle_acc(const double x1_in[], const double y1_in[], const double z1_in[], int n1_in,
                                const double x2_in[], const double y2_in[], const double z2_in [], int n2_in,
                                double x_out[], double y_out[], double z_out[]);

void get_cell_minmaxavg_latlons_acc( const int nlon, const int nlat, const double *lon, const double *lat, Cell *cell);

#pragma acc routine seq
void get_cell_vertices_acc( const int ij, const int nlon, const double *lon, const double *lat, double *x, double *y );

#endif
