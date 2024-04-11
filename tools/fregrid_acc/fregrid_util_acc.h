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
#ifndef FREGRID_UTIL_ACC_H_
#define FREGRID_UTIL_ACC_H_

#include "globals.h"

void get_output_grid_by_size_acc(int ntiles, Grid_config *grid, double lonbegin, double lonend, double latbegin, double latend,
                             int nlon, int nlat, int finer_steps, int center_y, unsigned int opcode);

void setup_vertical_interp_acc(VGrid_config *vgrid_in, VGrid_config *vgrid_out);

void do_vertical_interp_acc(VGrid_config *vgrid_in, VGrid_config *vgrid_out, Grid_config *grid_out, Field_config *field, int varid);

#endif
