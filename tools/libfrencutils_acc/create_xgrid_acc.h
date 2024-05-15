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
#ifndef CREATE_XGRID_ACC_H_
#define CREATE_XGRID_ACC_H_

/* this value is small compare to earth area */

int create_xgrid_2dx2d_order1_acc(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
			      const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
			      int *i_in, int *j_in, int *i_out, int *j_out, double *xgrid_area);

int create_xgrid_2dx2d_order2_acc(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
			      const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
			      int *i_in, int *j_in, int *i_out, int *j_out, double *xgrid_area, double *xgrid_clon, double *xgrid_clat);

int create_xgrid_great_circle_acc(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
            const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
            int *i_in, int *j_in, int *i_out, int *j_out, double *xgrid_area, double *xgrid_clon, double *xgrid_clat);


#endif
