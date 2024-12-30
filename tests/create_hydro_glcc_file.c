#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <netcdf.h>
#include "constant.h"

#define FILE_NAME "hydro_glcc.nc"
#define NLON 2880
#define NLAT 1440
#define DLON 360.0 / NLON
#define DLAT 180.0 / NLAT
#define LON_MIN -180.0 + DLON / 2.0
#define LAT_MIN -90.0 + DLAT / 2.0
#define MISSING_VALUE -999
#define WET_HALO 5
#define LAT_WBODIES 20
#define LON_WBODIES 40

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
    int waterbod_varid;
    check_err(nc_def_var(ncid, "WaterBod", NC_INT, 2, (int[]){lat_dimid, lon_dimid}, &waterbod_varid), __LINE__);
    check_err(nc_put_att_text(ncid, waterbod_varid, "units", 4, "none"), __LINE__);
    check_err(nc_put_att_text(ncid, waterbod_varid, "long_name", 12, "Water Bodies"), __LINE__);
    check_err(nc_put_att_int(ncid, waterbod_varid, "missing_value", NC_INT, 1, (int[]){MISSING_VALUE}), __LINE__);
    check_err(nc_put_att_double(ncid, waterbod_varid, "scale_factor", NC_DOUBLE, 1, (double[]){0.00444444455206394}), __LINE__);
    int pwetland_varid;
    check_err(nc_def_var(ncid, "PWetland", NC_INT, 2, (int[]){lat_dimid, lon_dimid}, &pwetland_varid), __LINE__);
    check_err(nc_put_att_text(ncid, pwetland_varid, "units", 4, "none"), __LINE__);
    check_err(nc_put_att_text(ncid, pwetland_varid, "long_name", 12, "Permanent Wetlands"), __LINE__);
    check_err(nc_put_att_int(ncid, pwetland_varid, "missing_value", NC_INT, 1, (int[]){MISSING_VALUE}), __LINE__);
    check_err(nc_put_att_double(ncid, pwetland_varid, "scale_factor", NC_DOUBLE, 1, (double[]){0.00444444455206394}), __LINE__);

    // End define mode
    check_err(nc_enddef(ncid), __LINE__);

    // Calculate the latitude and longitude values
    float *lon = malloc(NLON * sizeof(float));
    float *lat = malloc(NLAT * sizeof(float));
    for (int i = 0; i < NLON; i++)
        lon[i] = (float)i * DLON + LON_MIN;
    for (int i = 0; i < NLAT; i++)
        lat[i] = (float)i * DLAT + LAT_MIN;
    // Write the latitude and longitude values
    check_err(nc_put_var_float(ncid, lon_varid, lon), __LINE__);
    check_err(nc_put_var_float(ncid, lat_varid, lat), __LINE__);

    // Variable to hold water bodies and wetlands data
    int *data = malloc(NLON * NLAT * sizeof(int));
    // Wetlands data
    for (int i = 0; i < NLAT; i++)
        for (int j = 0; j < NLON; j++)
            if (((i > 0 && i < WET_HALO) || (i > NLAT - WET_HALO - 1 && i < NLAT - 1))
                || ((j > 0 && j < WET_HALO) || (j > NLON - WET_HALO - 1 && j < NLON - 1)))
                data[i * NLON + j] = 1;
            else
                data[i * NLON + j] = 0;
    check_err(nc_put_var_int(ncid, pwetland_varid, data), __LINE__);

    // Water bodies data
    for (int i = 0; i < NLAT; i++)
    {
        int i_write = 0;
        if ((i+1)%LAT_WBODIES == 0 && (i+1)%LON_WBODIES != 0)
            i_write = 1;
        for (int j = 0; j < NLON; j++)
        {
            int j_write = 0;
            if ((j+1)%LON_WBODIES == 0)
                j_write = 1;
            // Default value
            data[i * NLON + j] = 0;
            // The outermost cells are ocean
            if (i == 0 || i == NLAT -1)
                data[i * NLON + j] = 1;
            if (j==0 || j == NLON - 1)
                data[i * NLON + j] = 1;
            // Cells that flow east-west
            if (i_write == 1 && (j < i+1 && j < NLAT-i-1))
                data[i * NLON + j] = 1;
            if (i_write == 1 && (j > NLON-i-2 && j > NLON-NLAT+i))
                data[i * NLON + j] = 1;
            // Cells that flow north-south
            if (j_write == 1 && (i < j+1 && i < NLON-j-1))
                data[i * NLON + j] = 1;
            if (j_write == 1 && (i > NLAT-j-2 && i > NLAT-NLON+j))
                data[i * NLON + j] = 1;
        }
    }
    check_err(nc_put_var_int(ncid, waterbod_varid, data), __LINE__);
    free(data);

    // Close the file
    check_err(nc_close(ncid), __LINE__);
    printf("Successfully created %s\n", FILE_NAME);
    free(lon);
    free(lat);
    return EXIT_SUCCESS;
}