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
/*
  Copyright (C) 2011 NOAA Geophysical Fluid Dynamics Lab, Princeton, NJ
*/
#ifndef INTERP_H_
#define INTERP_H_
/*********************************************************************
                     interp.h
   This header files contains defition of some interpolation routine  (1-D or 2-D).
   contact: Zhi.Liang@noaa.gov
*********************************************************************/
void cubic_spline_sp(int size1, int size2, const double *grid1, const double *grid2, const double *data1,
                  double *data2 );
void cubic_spline(int size1, int size2, const double *grid1, const double *grid2, const double *data1,
		  double *data2, double yp1, double ypn  );
void conserve_interp(int nx_src, int ny_src, int nx_dst, int ny_dst, const double *x_src,
		     const double *y_src, const double *x_dst, const double *y_dst,
		     const double *mask_src, const double *data_src, double *data_dst );
void conserve_interp_great_circle(int nx_src, int ny_src, int nx_dst, int ny_dst, const double *x_src,
		     const double *y_src, const double *x_dst, const double *y_dst,
		     const double *mask_src, const double *data_src, double *data_dst );
void linear_vertical_interp(int nx, int ny, int nk1, int nk2, const double *grid1, const double *grid2,
			    double *data1, double *data2);
#endif
