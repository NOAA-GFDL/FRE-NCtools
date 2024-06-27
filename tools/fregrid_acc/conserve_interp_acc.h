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
#ifndef CONSERVE_INTERP_ACC_H_
#define CONSERVE_INTERP_ACC_H_

#include "globals_acc.h"

void setup_conserve_interp_acc(int ntiles_in, Grid_config *grid_in, int ntiles_out,
			   Grid_config *grid_out, Interp_config_acc *interp_acc, unsigned int opcode);

void do_scalar_conserve_interp_acc(Interp_config_acc *interp_acc, int varid, int ntiles_in, const Grid_config *grid_in,
			       int ntiles_out, const Grid_config *grid_out, const Field_config *field_in,
			       Field_config *field_out, unsigned int opcode, int nz);

void do_scalar_conserve_interp_order2_acc(Interp_config_acc *interp_acc, int varid, int ntiles_in, const Grid_config *grid_in,
			       int ntiles_out, const Grid_config *grid_out, const Field_config *field_in,
			       Field_config *field_out, unsigned int opcode, int nz);

void read_remap_file_acc(int ntiles_input_grid, int ntiles_output_grid,
                         Grid_config *output_grid, Grid_config *input_grid,
                         Interp_config_acc *interp_acc, unsigned int opcode);

void write_remap_file(const int ntiles_out, const int ntiles_in, Grid_config *output_grid,
                      Grid_config *input_grid, Interp_config_acc *interp_acc, unsigned int opcode);

void check_area_conservation(const int ntiles_output_grid, const int ntiles_input_grid, Grid_config *output_grid,
                             Interp_config_acc *interp_acc);

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

#endif
