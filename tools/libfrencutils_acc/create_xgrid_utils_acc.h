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
#ifndef CREATE_XGRID_ACC_H_
#define CREATE_XGRID_ACC_H_

#include "globals.h"

#define MV 50
/* this value is small compare to earth area */

void get_grid_area_acc(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area);

void get_grid_great_circle_area(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area);

void get_cell_minmaxavg_latlons( const int nlon, const int nlat, const double *lon, const double *lat,
                                 Minmaxavg_list *minmaxavg_list );

#pragma acc routine seq
void get_cell_vertices( const int ij, const int nlon, const double *lon, const double *lat, double *x, double *y );

#endif
