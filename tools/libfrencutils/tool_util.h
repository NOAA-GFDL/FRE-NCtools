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
/***********************************************************************
                      tool_util.h
    This header file provide some utilities routine that will be used in many tools.

    contact: Zhi.Liang@noaa.gov
***********************************************************************/
#ifndef TOOL_UTIL_H_
#define TOOL_UTIL_H_
#define MAX_GRID_LENGTH (10000)
#define VERSION_1 1
#define VERSION_2 2
#define VERSION_3 3
int round_to_nearest_int(double r);
void get_file_path(const char *file, char *dir);
int get_int_entry(char *line, int *value);
int get_double_entry(char *line, double *value);
double spherical_dist(double x1, double y1, double x2, double y2);
double spherical_area(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4 );
double bipolar_dist(double x1, double y1, double x2, double y2, double bpeq, double bpsp, double bpnp, double rp );
double bipolar_area(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4 );
void tp_trans(double *lon, double *lat, double lon_ref, double lon_start,
              double lam0, double bpeq, double bpsp, double bpnp, double rp );
double* compute_grid_bound(int nb, const double *bnds, const int *npts, int *grid_size, const char *center);
double* compute_grid_bound_legacy(int nb, const double *bnds, const double *dbnds, double stretch, int *grid_size, const char *center);
int get_legacy_grid_size(int nb, const double *bnds, const double *dbnds);
void get_boundary_type( const char *grid_file, int grid_version, int *cyclic_x, int *cyclic_y, int *is_tripolar );

void print_provenance_gv_gca(int fid,  const char * history, char * grid_version, int gca_flag);
void print_provenance_gv(int fid,  const char * history, char * grid_version);
void print_provenance(int fid, const char * history);
#endif
