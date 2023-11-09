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
#ifndef CONSERVE_INTERP_UTIL_H_
#define CONSERVE_INTERP_UTIL_H_
#include "globals.h"

void read_remap_file( int ntiles_in, int ntiles_out, Grid_config *grid_out,
                      Interp_config *interp, unsigned int opcode);

void malloc_xgrid_arrays( int nsize, int **i_in, int **j_in, int **i_out, int **j_out,
                          double **xgrid_area, double **xgrid_clon, double **xgrid_clat );

void get_CellStruct(const int tile_in, const int nx_in, const int nxgrid, int *i_in, int *j_in,
                    double *xgrid_area, double *xgrid_clon, double *xgrid_clat,
                    CellStruct *cell_in);

void get_interp( const int opcode, const int nxgrid, Interp_config *interp, const int m, const int n,
                 const int *i_in, const int *j_in, const int *i_out, const int *j_out,
                 const double *xgrid_clon, const double *xgrid_clat, const double *xgrid_area ) ;

void get_jstart_jend( const int nx_out, const int ny_out, const int nx_in, const int ny_in,
                      const double *lat_out, const double *lat_in,
                      int *jstart, int *jend, int *ny_now ) ;


#endif
