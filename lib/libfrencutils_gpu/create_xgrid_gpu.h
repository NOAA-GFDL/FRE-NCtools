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
#ifndef CREATE_XGRID_GPU_H_
#define CREATE_XGRID_GPU_H_

#include "globals_gpu.h"

int get_upbound_nxcells_2dx2d_gpu(const int nlon_input_cells,  const int nlat_input_cells,
                                  const int nlon_output_cells, const int nlat_output_cells,
                                  const int jlat_overlap_starts, const int jlat_overlap_ends,
                                  const double *input_grid_lon, const double *input_grid_lat,
                                  const double *output_grid_lon, const double *output_grid_lat,
                                  const double *skip_input_cells,
                                  const Grid_cells_struct_config *output_grid_cells,
                                  int *approx_nxcells_per_ij1, int *ij2_start, int *ij2_end);

int create_xgrid_2dx2d_order1_gpu(const int nlon_input_cells,  const int nlat_input_cells,
                                  const int nlon_output_cells, const int nlat_output_cells,
                                  const int jlat_overlap_starts, const int jlat_overlap_ends,
                                  const double *input_grid_lon,  const double *input_grid_lat,
                                  const double *output_grid_lon, const double *output_grid_lat,
                                  const int upbound_nxcells, const double *skip_input_cells,
                                  const Grid_cells_struct_config *output_grid_cells,
                                  int *approx_nxcells_per_ij1, int *ij2_start, int *ij2_end,
                                  Interp_per_input_tile *interp_for_input_tile);

int create_xgrid_2dx2d_order2_gpu(const int nlon_input_cells,  const int nlat_input_cells,
                                  const int nlon_output_cells, const int nlat_output_cells,
                                  const int jlat_overlap_starts, const int jlat_overlap_ends,
                                  const double *input_grid_lon,  const double *input_grid_lat,
                                  const double *output_grid_lon, const double *output_grid_lat,
                                  const int upbound_nxcells, const double *skip_input_cells,
                                  const Grid_cells_struct_config *output_grid_cells,
                                  int *approx_nxcells_per_ij1, int *ij2_start, int *ij2_end,
                                  Interp_per_input_tile *interp_for_input_tile, double *readin_input_area);

int create_xgrid_great_circle_gpu(const int *nlon_input_cells, const int *nlat_input_cells,
                                  const int *nlon_output_cells, const int *nlat_output_cells,
                                  const double *input_grid_lon, const double *input_grid_lat,
                                  const double *output_grid_lon, const double *output_grid_lat,
                                  const double *skip_input_cells,
                                  int *i_in, int *j_in, int *i_out, int *j_out,
                                  double *xgrid_area, double *xgrid_clon, double *xgrid_clat);

#endif
