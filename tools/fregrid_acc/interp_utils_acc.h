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

void get_input_skip_cells(const int mask_size, double **skip_cells);

void copy_xgrid_to_device_acc( const int ntiles_in, const int ntiles_out, const Xgrid_config *xgrid,
                               const unsigned int opcode );

void get_bounding_indices(const int ref_nx, const int ref_ny, const int nx, const int ny,
                          const double *ref_lat, const double *lat, int *jstart, int *jend, int *new_ny);

#endif
