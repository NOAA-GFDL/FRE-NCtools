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
			   Grid_config *grid_out, Xgrid_config *xgrid, unsigned int opcode);
void do_scalar_conserve_interp_acc(Xgrid_config *xgrid, int varid, int ntiles_in, const Grid_config *grid_in,
			       int ntiles_out, const Grid_config *grid_out, const Field_config *field_in,
			       Field_config *field_out, unsigned int opcode, int nz);
void do_vector_conserve_interp_acc(Xgrid_config *xgrid, int varid, int ntiles_in, const Grid_config *grid_in, int ntiles_out,
                               const Grid_config *grid_out, const Field_config *u_in,  const Field_config *v_in,
                               Field_config *u_out, Field_config *v_out, unsigned int opcode);

void read_remap_file_acc(int ntiles_in, int ntiles_out, Grid_config *output_grid, Grid_config *input_grid,
                         Xgrid_config *xgrid, unsigned int opcode);

void write_remap_file(const int ntiles_out, const int ntiles_in, Grid_config *output_grid,
                      Grid_config *input_grid, Xgrid_config *xgrid, unsigned int opcode);

void check_area_conservation(const int ntiles_output_grid, const int ntiles_input_grid, Grid_config *output_grid,
                             Xgrid_config *xgrid);

#endif
