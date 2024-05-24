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

#include "globals_acc.h"

int create_xgrid_2dx2d_order1_acc(const int *nlon_input_cells, const int *nlat_input_cells,
                                  const int *nlon_output_cells, const int *nlat_output_cells,
                                  const double *input_grid_lon, const double *input_grid_lat,
                                  const double *output_grid_lon, const double *output_grid_lat,
                                  const double *skip_input_cells,
                                  int *i_in, int *j_in, int *i_out, int *j_out, double *xgrid_area);

int create_xgrid_2dx2d_order2_acc(const int *nlon_input_cells, const int *nlat_input_cells,
                                  const int *nlon_output_cells, const int *nlat_output_cells,
                                  const double *input_grid_lon, const double *input_grid_lat,
                                  const double *output_grid_lon, const double *output_grid_lat,
                                  const double *skip_input_cells,
                                  int *i_in, int *j_in, int *i_out, int *j_out, double *xgrid_area,
                                  double *xgrid_clon, double *xgrid_clat);

int create_xgrid_great_circle_acc(const int *nlon_input_cells, const int *nlat_input_cells,
                                  const int *nlon_output_cells, const int *nlat_output_cells,
                                  const double *input_grid_lon, const double *input_grid_lat,
                                  const double *output_grid_lon, const double *output_grid_lat,
                                  const double *skip_input_cells,
                                  int *i_in, int *j_in, int *i_out, int *j_out,
                                  double *xgrid_area, double *xgrid_clon, double *xgrid_clat);

#endif
