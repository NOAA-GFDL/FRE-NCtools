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
#ifndef FREGRID_UTILS_GPU_H_
#define FREGRID_UTILS_GPU_H_

#include "globals_gpu.h"

void copy_grid_to_device_gpu( const int npoints, const double *lat, const double *lon );

void delete_grid_from_device_gpu( const int npoints, const double *lat, const double *lon );

void copy_interp_to_device_gpu( const int ntiles_in, const int ntiles_out, const Interp_config_gpu *interp_gpu,
                                const unsigned int opcode );

void get_input_grid_mask_gpu(const int mask_size, double **input_grid_mask);

void free_input_grid_mask_gpu(const int mask_size, double **input_grid_mask);

void create_interp_per_intile_arrays_on_device_gpu(const int nxcells, const unsigned int opcode,
                                                   Interp_per_input_tile *interp_per_itile);

#endif
