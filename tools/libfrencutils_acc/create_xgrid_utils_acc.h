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

#include "globals_acc.h"

#define MV 50
/* this value is small compare to earth area */

void get_grid_area_acc(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area);

void get_grid_great_circle_area_acc(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area);

#pragma acc routine seq
void poly_ctrlon_acc(const double *x, const double *y, int n, double clon_in, double *crtlon);

#pragma acc routine seq
void poly_ctrlat_acc(const double *x, const double *y, int n, double *crtlat);

#pragma acc routine seq
int clip_2dx2d_acc(const double lon1_in[], const double lat1_in[], int n1_in,
                   const double lon2_in[], const double lat2_in[], int n2_in,
                   double lon_out[], double lat_out[]);

int clip_2dx2d_great_circle_acc(const double x1_in[], const double y1_in[], const double z1_in[], int n1_in,
                                const double x2_in[], const double y2_in[], const double z2_in [], int n2_in,
                                double x_out[], double y_out[], double z_out[]);

void get_grid_cells_struct_acc( const int nlon, const int nlat, const double *lon, const double *lat,
                                Grid_cells_struct_config *grid_cells);

void free_grid_cell_struct_acc( const int ncells, Grid_cells_struct_config *grid_cells);

#pragma acc routine seq
void get_cell_vertices_acc( const int ij, const int nlon, const double *lon, const double *lat, double *x, double *y );

void create_upbound_nxcells_arrays_on_device_acc(const int n, int **approx_nxcells_per_ij1,
                                                 int **ij2_start, int **ij2_end);

void free_upbound_nxcells_array_from_all_acc( const int n, int *approx_nxcells_per_ij1,
                                              int *ij2_start, int *ij2_end);

void free_output_grid_cell_struct_from_all_acc(const int n, Grid_cells_struct_config *grid_cells);

void copy_data_to_xgrid_on_device_acc(const int nxcells, const int input_ncells, const int upbound_nxcells,
                                      int *xcells_per_ij1, double *xcell_clon, double *xcell_clat,
                                      int *approx_xcells_per_ij1, int *parent_input_indices, int *parent_output_indices,
                                      double *xcell_areas, Interp_per_input_tile *interp_for_input_tile);
#endif
