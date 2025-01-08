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
#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define FILE_NAME "ocean_temp_salt.nc"
#define NX 360
#define NY 180
#define NZ 33

#define TIME 1
#define TEMP_MIN -2
#define TEMP_MAX 32
#define NX_MIN -10
#define NX_MAX 10
#define NY_MIN -10
#define NY_MAX 10

void check_err(int stat, int line) {
    if (stat != NC_NOERR) {
        printf("Error at line %d: %s\n", line, nc_strerror(stat));
        exit(1);
    }
}

int main(int argc, char **argv) {
    int ncid, xdim, ydim, zdim, tdim, bnds_dim;
    int xid, yid, zid, tid, bnds_id, temp_id, salt_id;
    float missing_value = (float)NC_FILL_FLOAT;
    // Create the file
    check_err(nc_create(FILE_NAME, NC_CLOBBER | NC_64BIT_OFFSET, &ncid), __LINE__);

    // Define dimensions
    check_err(nc_def_dim(ncid, "xaxis", NX, &xdim), __LINE__);
    check_err(nc_def_dim(ncid, "yaxis", NY, &ydim), __LINE__);
    check_err(nc_def_dim(ncid, "zaxis", NZ, &zdim), __LINE__);
    check_err(nc_def_dim(ncid, "bnds", 2, &bnds_dim), __LINE__);
    check_err(nc_def_dim(ncid, "Time", NC_UNLIMITED, &tdim), __LINE__);

    // Define variables
    check_err(nc_def_var(ncid, "xaxis", NC_DOUBLE, 1, &xdim, &xid), __LINE__);
    check_err(nc_put_att_text(ncid, xid, "units", 4, "none"), __LINE__);
    check_err(nc_put_att_text(ncid, xid, "long_name", 7, "xaxis_1"), __LINE__);
    check_err(nc_put_att_text(ncid, xid, "cartesian_axis", 1, "X"), __LINE__);

    check_err(nc_def_var(ncid, "yaxis", NC_DOUBLE, 1, &ydim, &yid), __LINE__);
    check_err(nc_put_att_text(ncid, yid, "units", 4, "none"), __LINE__);
    check_err(nc_put_att_text(ncid, yid, "long_name", 7, "yaxis_1"), __LINE__);
    check_err(nc_put_att_text(ncid, yid, "cartesian_axis", 1, "Y"), __LINE__);

    check_err(nc_def_var(ncid, "zaxis", NC_DOUBLE, 1, &zdim, &zid), __LINE__);
    check_err(nc_put_att_text(ncid, zid, "units", 1, "m"), __LINE__);
    check_err(nc_put_att_text(ncid, zid, "long_name", 7, "zaxis_1"), __LINE__);
    check_err(nc_put_att_text(ncid, zid, "cartesian_axis", 1, "Z"), __LINE__);
    check_err(nc_put_att_text(ncid, zid, "positive", 4, "down"), __LINE__);
    check_err(nc_put_att_text(ncid, zid, "bounds", 10, "zaxis_bnds"), __LINE__);

    check_err(nc_def_var(ncid, "zaxis_bnds", NC_DOUBLE, 2, (int[]){zdim, bnds_dim}, &bnds_id), __LINE__);

    check_err(nc_def_var(ncid, "Time", NC_DOUBLE, 1, &tdim, &tid), __LINE__);
    check_err(nc_put_att_text(ncid, tid, "units", 10, "time level"), __LINE__);
    check_err(nc_put_att_text(ncid, tid, "long_name", 4, "Time"), __LINE__);
    check_err(nc_put_att_text(ncid, tid, "cartesian_axis", 1, "T"), __LINE__);

    check_err(nc_def_var(ncid, "temp", NC_DOUBLE, 4, (int[]){tdim, zdim, ydim, xdim}, &temp_id), __LINE__);
    check_err(nc_put_att_text(ncid, temp_id, "units", 4, "none"), __LINE__);
    check_err(nc_put_att_text(ncid, temp_id, "long_name", 4, "temp"), __LINE__);
    check_err(nc_put_att_float(ncid, temp_id, "missing_value", NC_FLOAT, 1, &missing_value), __LINE__);

    check_err(nc_def_var(ncid, "salt", NC_DOUBLE, 4, (int[]){tdim, zdim, ydim, xdim}, &salt_id), __LINE__);
    check_err(nc_put_att_text(ncid, salt_id, "units", 4, "none"), __LINE__);
    check_err(nc_put_att_text(ncid, salt_id, "long_name", 4, "salt"), __LINE__);
    check_err(nc_put_att_float(ncid, salt_id, "missing_value", NC_FLOAT, 1, &missing_value), __LINE__);

    // Add Global attribute
    check_err(nc_put_att_text(ncid, NC_GLOBAL, "description", 15, "Test input file"), __LINE__);

    // End define mode
    check_err(nc_enddef(ncid), __LINE__);

    // Allocate memory for sample data and write to file
    double xdata[NX], ydata[NY];

    double zdata[NZ] = {   0,   10,   20,   30,   50,   75,  100,  125,  150,  200,  250,
                         300,  400,  500,  600,  700,  800,  900, 1000, 1100, 1200, 1300,
                        1400, 1500, 1750, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500 };

    for (int i = 0; i < NX; i++)
        xdata[i] = (float)i+0.5;
    for (int i = 0; i < NY; i++)
        ydata[i] = -89.5 +(float)(i)*180.0/(float)NY;
    double zdata_bnds[NZ][2];
    for (int i = 0; i < NZ; i++) {
        double up_bnd_diff, lo_bnd_diff;
        if (i == NZ - 1) {
            lo_bnd_diff = (zdata[i] - zdata[i-1])/2;
            up_bnd_diff = lo_bnd_diff;
        } else if (i == 0) {
            up_bnd_diff = (zdata[i+1] - zdata[i])/2;
            lo_bnd_diff = up_bnd_diff;
        } else {
            up_bnd_diff = (zdata[i+1] - zdata[i])/2;
            lo_bnd_diff = (zdata[i] - zdata[i-1])/2;
        }
        double up_bnd = (zdata[i+1] - zdata[i])/2;
        zdata_bnds[i][0] = zdata[i] - lo_bnd_diff;
        zdata_bnds[i][1] = zdata[i] + up_bnd_diff;
    }
    check_err(nc_put_var_double(ncid, xid, &xdata[0]), __LINE__);
    check_err(nc_put_var_double(ncid, yid, &ydata[0]), __LINE__);
    check_err(nc_put_var_double(ncid, zid, &zdata[0]), __LINE__);
    check_err(nc_put_var_double(ncid, bnds_id, &zdata_bnds[0][0]), __LINE__);

    double *temp_data = malloc(NZ*NY*NX*sizeof(double));
    double *salt_data = malloc(NZ*NY*NX*sizeof(double));
    double time_data;
    int i;
    for (int rec = 0; rec < TIME; rec++)
    {
        for (int z = 0; z < NZ; z++)
            for (int y = 0; y < NY; y++) {
                float yval = (float)NY_MIN + (float)(NY_MAX - NY_MIN) * y / (float)(NY - 1);
                for (int x = 0; x < NX; x++)
                {
                    float xval = (float)NX_MIN + (float)(NX_MAX - NX_MIN) * x / (float)(NX - 1);
                    i = z*NY*NX + y*NX + x;
                    temp_data[i] = (TEMP_MAX-TEMP_MIN) *
                                   ((float)x/NX + (float)y/NY + (float)z/NZ + (float)rec/TIME);
                    salt_data[i] = 20.0 * (xval * xval + yval * yval) - 1.0 - zdata[z]/4.0;
                    if (salt_data[i] < 0.0) salt_data[i] = (float)NC_FILL_FLOAT;
                }
            }
        time_data = (float)rec;
        check_err(nc_put_vara_double(ncid, tid,
                                    (size_t[1]){rec},
                                    (size_t[1]){1},
                                    &time_data), __LINE__);
        check_err(nc_put_vara_double(ncid, temp_id,
                                    (size_t[4]){rec, 0, 0, 0},
                                    (size_t[4]){1, NZ, NY, NX},
                                    &temp_data[0]), __LINE__);
        check_err(nc_put_vara_double(ncid, salt_id,
                                    (size_t[4]){rec, 0, 0, 0},
                                    (size_t[4]){1, NZ, NY, NX},
                                    &salt_data[0]), __LINE__);
    }
    // Close the file
    check_err(nc_close(ncid), __LINE__);
    free(temp_data);
    free(salt_data);
    printf("Netcdf file %s created successfully!\n", FILE_NAME);
    return 0;
}
