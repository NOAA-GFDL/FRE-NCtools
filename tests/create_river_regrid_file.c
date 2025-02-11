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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <netcdf.h>
#include "constant.h"

#define FILE_NAME "river_regrid_input.nc"
#define NLON 720
#define NLAT 360
#define LON_MIN 0.25
#define LAT_MIN -89.75
#define BASIN_DIVISIONS 4
#define MISSING_VALUE -9999.
#define XY_MIN -10.0
#define XY_MAX 10.0

void check_err(int stat, int line) {
    if (stat != NC_NOERR) {
        printf("Error at line %d: %s\n", line, nc_strerror(stat));
        exit(EXIT_FAILURE);
    }
}

int
main(int argc, char **argv)
{
    int ncid;
    int lon_dimid, lat_dimid;
    int lon_varid, lat_varid;

    check_err(nc_create(FILE_NAME, (NC_CLOBBER | NC_NETCDF4), &ncid), __LINE__);

    // Define the dimensions
    check_err(nc_def_dim(ncid, "lon", NLON, &lon_dimid), __LINE__);
    check_err(nc_def_dim(ncid, "lat", NLAT, &lat_dimid), __LINE__);

    // Define the dimension variables
    check_err(nc_def_var(ncid, "lon", NC_FLOAT, 1, &lon_dimid, &lon_varid), __LINE__);
    check_err(nc_put_att_text(ncid, lon_varid, "units", 12, "degrees_east"), __LINE__);
    check_err(nc_put_att_text(ncid, lon_varid, "long_name", 34, "Nominal Longitude of T-cell center"), __LINE__);
    check_err(nc_put_att_text(ncid, lon_varid, "cartesian_axis", 1, "X"), __LINE__);
    check_err(nc_def_var(ncid, "lat", NC_FLOAT, 1, &lat_dimid, &lat_varid), __LINE__);
    check_err(nc_put_att_text(ncid, lat_varid, "units", 13, "degrees_north"), __LINE__);
    check_err(nc_put_att_text(ncid, lat_varid, "long_name", 35, "Nominal Latitude of T-cell center"), __LINE__);
    check_err(nc_put_att_text(ncid, lat_varid, "cartesian_axis", 1, "Y"), __LINE__);

    // Define the other variables
    int basin_varid;
    check_err(nc_def_var(ncid, "basin", NC_INT, 2, (int[]){lat_dimid, lon_dimid}, &basin_varid), __LINE__);
    check_err(nc_put_att_text(ncid, basin_varid, "units", 4, "none"), __LINE__);
    check_err(nc_put_att_text(ncid, basin_varid, "long_name", 14, "River basin id"), __LINE__);
    check_err(nc_put_att_double(ncid, basin_varid, "missing_value", NC_DOUBLE, 1, (double[]){MISSING_VALUE}), __LINE__);
    int cellarea_varid;
    check_err(nc_def_var(ncid, "cellarea", NC_DOUBLE, 2, (int[]){lat_dimid, lon_dimid}, &cellarea_varid), __LINE__);
    check_err(nc_put_att_text(ncid, cellarea_varid, "units", 3, "m^2"), __LINE__);
    check_err(nc_put_att_text(ncid, cellarea_varid, "long_name", 9, "Cell area"), __LINE__);
    check_err(nc_put_att_double(ncid, cellarea_varid, "missing_value", NC_DOUBLE, 1, (double[]){MISSING_VALUE}), __LINE__);
    int suba_varid;
    check_err(nc_def_var(ncid, "subA", NC_DOUBLE, 2, (int[]){lat_dimid, lon_dimid}, &suba_varid), __LINE__);
    check_err(nc_put_att_text(ncid, suba_varid, "units", 3, "m^2"), __LINE__);
    check_err(nc_put_att_text(ncid, suba_varid, "long_name", 13, "subbasin area"), __LINE__);
    check_err(nc_put_att_double(ncid, suba_varid, "missing_value", NC_DOUBLE, 1, (double[]){MISSING_VALUE}), __LINE__);
    int celllength_varid;
    check_err(nc_def_var(ncid, "celllength", NC_DOUBLE, 2, (int[]){lat_dimid, lon_dimid}, &celllength_varid), __LINE__);
    check_err(nc_put_att_text(ncid, celllength_varid, "units", 1, "m"), __LINE__);
    check_err(nc_put_att_text(ncid, celllength_varid, "long_name", 11, "Cell length"), __LINE__);
    check_err(nc_put_att_double(ncid, celllength_varid, "missing_value", NC_DOUBLE, 1, (double[]){MISSING_VALUE}), __LINE__);
    int tocell_varid;
    check_err(nc_def_var(ncid, "tocell", NC_DOUBLE, 2, (int[]){lat_dimid, lon_dimid}, &tocell_varid), __LINE__);
    check_err(nc_put_att_text(ncid, tocell_varid, "units", 4, "none"), __LINE__);
    check_err(nc_put_att_text(ncid, tocell_varid, "long_name", 28, "direction to downstream cell"), __LINE__);
    check_err(nc_put_att_double(ncid, tocell_varid, "missing_value", NC_DOUBLE, 1, (double[]){MISSING_VALUE}), __LINE__);
    int travel_varid;
    check_err(nc_def_var(ncid, "travel", NC_DOUBLE, 2, (int[]){lat_dimid, lon_dimid}, &travel_varid), __LINE__);
    check_err(nc_put_att_text(ncid, travel_varid, "units", 4, "none"), __LINE__);
    check_err(nc_put_att_text(ncid, travel_varid, "long_name", 42, "cells left to travel before reaching ocean"), __LINE__);
    check_err(nc_put_att_double(ncid, travel_varid, "missing_value", NC_DOUBLE, 1, (double[]){MISSING_VALUE}), __LINE__);
    int land_frac_varid;
    check_err(nc_def_var(ncid, "land_frac", NC_DOUBLE, 2, (int[]){lat_dimid, lon_dimid}, &land_frac_varid), __LINE__);
    check_err(nc_put_att_text(ncid, land_frac_varid, "units", 4, "none"), __LINE__);
    check_err(nc_put_att_text(ncid, land_frac_varid, "long_name", 23, "land/sea mask(land = 1)"), __LINE__);
    check_err(nc_put_att_double(ncid, land_frac_varid, "missing_value", NC_DOUBLE, 1, (double[]){MISSING_VALUE}), __LINE__);

    // End define mode
    check_err(nc_enddef(ncid), __LINE__);

    // Calculate the latitude and longitude values
    float *lon = malloc(NLON * sizeof(float));
    float *lat = malloc(NLAT * sizeof(float));
    for (int i = 0; i < NLON; i++)
        lon[i] = (float)i * 0.5 + LON_MIN;
    for (int i = 0; i < NLAT; i++)
        lat[i] = (float)i * 0.5 + LAT_MIN;
    // Write the latitude and longitude values
    check_err(nc_put_var_float(ncid, lon_varid, lon), __LINE__);
    check_err(nc_put_var_float(ncid, lat_varid, lat), __LINE__);

    // Calculate the basin values.  Using a step function for simplicity
    double *basin = malloc(NLON * NLAT * sizeof(double));
    for (int i = 0; i < NLAT; i++)
        for (int j = 0; j < NLON; j++)
        {
            if (i == 0 || j == 0 || i == NLAT - 1 || j == NLON - 1)
                basin[i * NLON + j] = (double)MISSING_VALUE;
            else
                basin[i * NLON + j] = (double)(1 + i / (NLAT / BASIN_DIVISIONS) * NLON / (NLON / BASIN_DIVISIONS)
                    + j / (NLON / BASIN_DIVISIONS));
        }
    // Write the basin values
    check_err(nc_put_var_double(ncid, basin_varid, basin), __LINE__);
    free(basin);

    // Calculate the cell area values
    double *cellarea = malloc(NLON * NLAT * sizeof(double));
    // Since the lat/lon use a spacing of 0.5 degrees, we'll use half of that
    // to calculate the area of each cell
    for (int i = 0; i < NLAT; i++)
        for (int j = 0; j < NLON; j++)
            cellarea[i * NLON + j] = RADIUS * RADIUS * (sin((double)(lat[i] + 0.25) * D2R) - sin((double)(lat[i] - 0.25) * D2R))
                * 0.5 * D2R;
    check_err(nc_put_var_double(ncid, cellarea_varid, cellarea), __LINE__);
    // Using the cellarea as the subbasin area, except for the boundary cells
    for (int i = 0; i < NLAT; i++)
        for (int j = 0; j < NLON; j++)
            cellarea[i * NLON + j] = (i == 0 || j == 0 || i == NLAT - 1 || j == NLON - 1) ? MISSING_VALUE : cellarea[i * NLON + j];
    check_err(nc_put_var_double(ncid, suba_varid, cellarea), __LINE__);
    free(cellarea);

    double *celllength = malloc(NLON * NLAT * sizeof(double));
    // Using the great circle distance formula
    for (int i = 0; i < NLAT; i++)
        for (int j = 0; j < NLON; j++)
            celllength[i * NLON + j] = RADIUS * acos(cos((double)(lat[i] + 0.25) * D2R) * cos((double)(lat[i] - 0.25) * D2R) * cos((double)(0.5 * D2R))
                + sin((double)(lat[i] + 0.25) * D2R) * sin((double)(lat[i] - 0.25) * D2R));
    // Write the cell area values, since
    check_err(nc_put_var_double(ncid, celllength_varid, celllength), __LINE__);
    free(celllength);

    double *tocell = malloc(NLON * NLAT * sizeof(double));
    for (int i = 0; i < NLAT; i++)
    {
        double yval = XY_MIN + (XY_MAX - XY_MIN) * (double)i / (double)(NLAT);
        for (int j = 0; j < NLON; j++)
        {
            double xval= XY_MIN + (XY_MAX - XY_MIN) * (double)j / (double)(NLON);
            if (i == 0 || j == 0 || i == NLAT - 1 || j == NLON - 1)
                tocell[i * NLON + j] = MISSING_VALUE;
            else
            {
                double exp = (int)(8.0 * sin(xval*xval+yval*yval) * sin(xval*xval+yval*yval));
                tocell[i * NLON + j] = pow(2.0, exp);
            }
        }
    }
    check_err(nc_put_var_double(ncid, tocell_varid, tocell), __LINE__);
    free(tocell);

    double *travel = malloc(NLON * NLAT * sizeof(double));
    for (int i = 0; i < NLAT; i++)
        for (int j = 0; j < NLON; j++)
        {
            if (i == 0 || j == 0 || i == NLAT - 1 || j == NLON - 1)
                travel[i * NLON + j] = MISSING_VALUE;
            else
                travel[i * NLON + j] = fmin(fmin(i, NLAT - i - 1), fmin(j, NLON - j - 1));
        }
    check_err(nc_put_var_double(ncid, travel_varid, travel), __LINE__);
    free(travel);

    double *land_frac = malloc(NLON * NLAT * sizeof(double));
    for (int i = 0; i < NLAT; i++)
        for (int j = 0; j < NLON; j++)
            land_frac[i * NLON + j] = (i == 0 || j == 0 || i == NLAT - 1 || j == NLON - 1) ? 0.0 : 1.0;
    check_err(nc_put_var_double(ncid, land_frac_varid, land_frac), __LINE__);
    free(land_frac);

    // Close the file
    check_err(nc_close(ncid), __LINE__);
    printf("Successfully created %s\n", FILE_NAME);
    free(lon);
    free(lat);
    return EXIT_SUCCESS;
}