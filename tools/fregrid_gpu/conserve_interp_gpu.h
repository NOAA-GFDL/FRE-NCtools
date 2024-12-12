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
#ifndef CONSERVE_INTERP_GPU_H_
#define CONSERVE_INTERP_GPU_H_

#include "globals_gpu.h"

void setup_conserve_interp_gpu(int ntiles_in, Grid_config *grid_in, int ntiles_out,
			   Grid_config *grid_out, Interp_config_gpu *interp_gpu, unsigned int opcode);

void do_scalar_conserve_interp_gpu(Interp_config_gpu *interp_gpu, int varid, int ntiles_in, const Grid_config *grid_in,
			       int ntiles_out, const Grid_config *grid_out, const Field_config *field_in,
                                   Field_config *field_out, unsigned int opcode);

void read_remap_file_gpu(int ntiles_input_grid, int ntiles_output_grid,
                         Grid_config *output_grid, Grid_config *input_grid,
                         Interp_config_gpu *interp_gpu, unsigned int opcode);

void write_remap_file(const int ntiles_out, const int ntiles_in, Grid_config *output_grid,
                      Grid_config *input_grid, Interp_config_gpu *interp_gpu, unsigned int opcode);

void check_area_conservation(const int ntiles_output_grid, const int ntiles_input_grid, Grid_config *output_grid,
                             Interp_config_gpu *interp_gpu);

void get_input_area_weight(const int weights_exist, const int cell_measures, const int cell_methods,
                           const Field_config *field_in, const Grid_config *input_grid,
                           double *input_area_weight);

void interp_data_order1(const Grid_config *output_grid, const Grid_config *input_grid,
                        Interp_per_input_tile *interp_for_itile, double *input_area_weight, double *fieldin_data,
                        double *fieldout_data, double *out_area, int *out_miss, double missing);

void interp_data_order2( const Grid_config *output_grid, const Grid_config *input_grid,
                         Interp_per_input_tile *interp_for_itile, double *input_area_weight, double *fieldin_data,
                         double *fieldout_data, double *out_area, int *out_miss,
                         int *grad_mask, double *grad_y, double *grad_x, double missing);

void get_bounding_indices_gpu(const int ref_nlon_cells, const int ref_nlat_cells,
                              const int nlon_cells, const int nlat_cells,
                              const double *ref_grid_lat, const double *grid_lat,
                              int *overlap_starts_here_index, int *overlap_ends_here_index);

void create_interp_gpu_itile_arrays_on_device_gpu(const int nxcells, const unsigned int opcode,
                                                  Interp_per_input_tile *interp_per_itile);

#endif
