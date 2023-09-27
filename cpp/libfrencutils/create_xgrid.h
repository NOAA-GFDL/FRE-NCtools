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
#ifndef CREATE_XGRID_H_
#define CREATE_XGRID_H_

#include <vector>
#include <array>
#include <span>

#define MV 50

#include "constant.h"
#include "BBox3D.h"
using std::vector;



/* this value is small compare to earth area */

//inline constexpr unsigned int get_MAXXGRID(){
inline
unsigned int get_MAXXGRID(){
  std::cout << std::endl << "*** Hi from get_MAXXDRID*** " << std::endl;
  return 10000000;
}


double poly_ctrlon(const double lon[], const double lat[], int n, double clon);
double poly_ctrlat(const double lon[], const double lat[], int n);
double box_ctrlon(double ll_lon, double ll_lat, double ur_lon, double ur_lat, double clon);
double box_ctrlat(double ll_lon, double ll_lat, double ur_lon, double ur_lat);
int get_maxxgrid(void);
void get_grid_area(const int nlon, const int nlat, const double *lon, const double *lat, double *area);
void get_grid_great_circle_area(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area);
void get_grid_area_dimensionless(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area);
void get_grid_area_no_adjust(const int nlon, const int nlat, const double *lon, const double *lat, double *area);
int clip(const double lon_in[], const double lat_in[], int n_in, double ll_lon, double ll_lat,
	 double ur_lon, double ur_lat, double lon_out[], double lat_out[]);
void pimod(double x[],int nn);
int clip_2dx2d(const double lon1_in[], const double lat1_in[], int n1_in, 
	       const double lon2_in[], const double lat2_in[], int n2_in, 
	       double lon_out[], double lat_out[]);
int create_xgrid_1dx2d_order1(const int nlon_in, const int nlat_in, const int nlon_out, const int nlat_out, const double *lon_in,
                              const double *lat_in, const double *lon_out, const double *lat_out,
                              const double *mask_in, int *i_in, int *j_in, int *i_out,
                              int *j_out, double *xgrid_area);
int create_xgrid_1dx2d_order2(const int nlon_in, const int nlat_in, const int nlon_out, const int nlat_out,
                              const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                              const double *mask_in, int *i_in, int *j_in, int *i_out, int *j_out,
                              double *xgrid_area, double *xgrid_clon, double *xgrid_clat);
int create_xgrid_2dx1d_order1(const int nlon_in, const int nlat_in, const int nlon_out, const int nlat_out, const double *lon_in,
                              const double *lat_in, const double *lon_out, const double *lat_out,
                              const double *mask_in, int *i_in, int *j_in, int *i_out,
                              int *j_out, double *xgrid_area);
int create_xgrid_2dx1d_order2(const int nlon_in, const int nlat_in, const int nlon_out, const int nlat_out,
                              const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                              const double *mask_in, int *i_in, int *j_in, int *i_out, int *j_out,
                              double *xgrid_area, double *xgrid_clon, double *xgrid_clat);
int create_xgrid_2dx2d_order1(const int nlon_in, const int nlat_in, const int nlon_out, const int nlat_out,
                              const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                              const double *mask_in, int *&i_in, int *&j_in, int *&i_out,
                              int *&j_out, double *&xgrid_area);
int create_xgrid_2dx2d_order2(const int nlon_in, const int nlat_in, const int nlon_out, const int nlat_out,
                              const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                              const double *mask_in, int *&i_in, int *&j_in, int *&i_out, int *&j_out,
                              double *&xgrid_area, double *&xgrid_clon, double *&xgrid_clat);
int clip_2dx2d_great_circle(const double x1_in[], const double y1_in[], const double z1_in[], int n1_in, 
			    const double x2_in[], const double y2_in[], const double z2_in [], int n2_in, 
			    double x_out[], double y_out[], double z_out[]);
int create_xgrid_great_circle(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
			      const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
			      const double *mask_in, int *&i_in, int *&j_in, int *&i_out, int *&j_out,
			      double *&xgrid_area, double *&xgrid_clon, double *&xgrid_clat);


nct::BBox3D getBoxForSphericalPolygon(const double lat_m[], const double lon_m[],
                               const std::array<size_t, 4> &is, bool debugf  = false);

void  create_xgrid_2dx2d_st(const int nlon_in, const int nlat_in, const int nlon_out, const int nlat_out,
                            const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                            const double *mask_in, std::vector<size_t> &i_in, std::vector<size_t> &j_in,
                            std::vector<size_t> &i_out, std::vector<size_t> &j_out, std::vector<double> &xgrid_area,
                            std::vector<double> &xgrid_clon, std::vector<double> &xgrid_clat,
                            int order);

void  create_xgrid_2dx2d_bf(const int nlon_in, const int nlat_in, const int nlon_out, const int nlat_out,
                            const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                            const double *mask_in, std::vector<size_t> &i_in, std::vector<size_t> &j_in,
                            std::vector<size_t> &i_out, std::vector<size_t> &j_out, std::vector<double> &xgrid_area,
                            std::vector<double> &xgrid_clon, std::vector<double> &xgrid_clat,
                            int order);

inline
void latlon2xyz(const double lat, const double lon,  std::array<double,3> &  v){
  v[0] = RADIUS * cos(lat) * cos(lon );
  v[1] = RADIUS * cos(lat) * sin(lon);
  v[2] = RADIUS * sin(lat);
}

inline
size_t pt_idx(const size_t i, const size_t j,  const size_t nx) {
  return ( j * nx + i);
}

/**
 *  Generate four indices into a 1D array of points; such that the data of these
 *  four points represent a counter-clockwise grid cell. The 1D array of points
 *  can be though of as a mesh of 2D (lat-long) grid points.
 *  array.
 * @param i  i lon index of the lower left cell
 * @param j  j or lat index of the lower left cell
 * @param NX  With in number of points in the 2D grid.
 * @return an array of the indices
 */
inline std::array<size_t, 4>
get_cell_idxs_ccw_4(const size_t i, const size_t j, const size_t nx) {
  std::array<size_t, 4> idxs;
  idxs[0] = pt_idx(i, j, nx); //ll
  idxs[1] = pt_idx(i + 1, j , nx); //lr
  idxs[2] = pt_idx(i + 1, j + 1, nx); //ur
  idxs[3] = pt_idx(i, j + 1, nx);//ul
  return idxs;
}

#endif
