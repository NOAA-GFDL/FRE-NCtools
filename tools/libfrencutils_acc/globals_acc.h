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
#ifndef GLOBALS_ACC_H_
#define GLOBALS_ACC_H_
#include "globals.h"
#include "parameters.h"

typedef struct {
  int nxcells;
  int *input_parent_cell_index;
  int *output_parent_cell_index;
  double *xcell_area;
  double *dcentroid_lon;
  double *dcentroid_lat;
} Interp_per_input_tile;

typedef struct {
  size_t nxcells;
  int *input_parent_lon_index;
  int *input_parent_lat_index;
  int *output_parent_lon_index;
  int *output_parent_lat_index;
  int *output_parent;
  int *input_parent_tile;
  double *dcentroid_lon;
  double *dcentroid_lat;
  Interp_per_input_tile *mtile;
  double *xcell_area;
  double *weight;
  int    *index;
  char   remap_file[STRING];
  int    file_exist;
} Interp_config_acc;

typedef struct {
  double *lon_min;
  double *lon_max;
  double *lat_min;
  double *lat_max;
  double *lon_cent;
  double *area;
  int *nvertices;
  double **lon_vertices;
  double **lat_vertices;
  double *recomputed_area;
  double *centroid_lon;
  double *centroid_lat;
  double *skip_cells;

} Grid_cells_struct_config;

#endif
