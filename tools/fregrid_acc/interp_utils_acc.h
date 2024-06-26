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
#ifndef FREGRID_UTILS_ACC_H_
#define FREGRID_UTILS_ACC_H_

#include "globals_acc.h"

void copy_grid_to_device_acc( const int npoints, const double *lat, const double *lon );

void delete_grid_from_device_acc( const int npoints, const double *lat, const double *lon );

void copy_interp_to_device_acc( const int ntiles_in, const int ntiles_out, const Interp_config_acc *interp_acc,
                                const unsigned int opcode );

void get_bounding_indices_acc(const int ref_nlon_cells, const int ref_nlat_cells,
                              const int nlon_cells, const int nlat_cells,
                              const double *ref_grid_lat, const double *grid_lat,
                              int *overlap_starts_here_index, int *overlap_ends_here_index);

void get_input_grid_mask_acc(const int mask_size, double **input_grid_mask);

void free_input_grid_mask_acc(const int mask_size, double **input_grid_mask);

void create_interp_per_intile_arrays_on_device_acc(const int nxcells, const unsigned int opcode,
                                                   Interp_per_input_tile *interp_per_mile);

#endif
