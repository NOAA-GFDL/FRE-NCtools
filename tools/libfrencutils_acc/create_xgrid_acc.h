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
#ifndef MAXXGRID
#define MAXXGRID 1e6
#endif

#define MV 50
/* this value is small compare to earth area */

double poly_ctrlon_acc(const double lon[], const double lat[], int n, double clon);

double poly_ctrlat_acc(const double lon[], const double lat[], int n);

double box_ctrlon_acc(double ll_lon, double ll_lat, double ur_lon, double ur_lat, double clon);

double box_ctrlat_acc(double ll_lon, double ll_lat, double ur_lon, double ur_lat);

int get_maxxgrid_acc(void);

void get_grid_area_acc(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area);

void get_grid_great_circle_area_acc(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area);

//void get_grid_area_dimensionless(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area);

void get_grid_area_no_adjust_acc(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area);

int clip_acc(const double lon_in[], const double lat_in[], int n_in, double ll_lon, double ll_lat,
	 double ur_lon, double ur_lat, double lon_out[], double lat_out[]);

void pimod_acc(double x[],int nn);

int clip_2dx2d_acc(const double lon1_in[], const double lat1_in[], int n1_in,
	       const double lon2_in[], const double lat2_in[], int n2_in,
	       double lon_out[], double lat_out[]);

int create_xgrid_1dx2d_order_acc(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out, const double *lon_in,
			      const double *lat_in, const double *lon_out, const double *lat_out,
			      const double *mask_in, int *i_in, int *j_in, int *i_out,
			      int *j_out, double *xgrid_area);

int create_xgrid_1dx2d_order2_acc(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
			      const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
			      const double *mask_in, int *i_in, int *j_in, int *i_out, int *j_out,
			      double *xgrid_area, double *xgrid_clon, double *xgrid_clat);

int create_xgrid_2dx1d_order1_acc(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out, const double *lon_in,
			      const double *lat_in, const double *lon_out, const double *lat_out,
			      const double *mask_in, int *i_in, int *j_in, int *i_out,
            int *j_out, double *xgrid_area);

int create_xgrid_2dx1d_order2_acc(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
			      const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
			      const double *mask_in, int *i_in, int *j_in, int *i_out, int *j_out,
            double *xgrid_area, double *xgrid_clon, double *xgrid_clat);

int create_xgrid_2dx2d_order1_acc(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
			      const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
			      const double *mask_in, int *i_in, int *j_in, int *i_out,
            int *j_out, double *xgrid_area);

int create_xgrid_2dx2d_order2_acc(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
			      const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
			      const double *mask_in, int *i_in, int *j_in, int *i_out, int *j_out,
            double *xgrid_area, double *xgrid_clon, double *xgrid_clat);

int clip_2dx2d_great_circle_acc(const double x1_in[], const double y1_in[], const double z1_in[], int n1_in,
			    const double x2_in[], const double y2_in[], const double z2_in [], int n2_in,
          double x_out[], double y_out[], double z_out[]);

int create_xgrid_great_circle_acc(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
			      const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
			      const double *mask_in, int *i_in, int *j_in, int *i_out, int *j_out,
			      double *xgrid_area, double *xgrid_clon, double *xgrid_clat);
