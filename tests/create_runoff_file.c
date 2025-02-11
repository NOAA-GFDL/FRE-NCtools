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
#include <math.h>

#define FILE_NAME "runoff_input.nc"
#define NX 360
#define NY 180

#define NX_MIN -10.0
#define NX_MAX 10.0
#define NY_MIN -10.0
#define NY_MAX -10.0
#define TIME 1

void check_err(int stat, int line) {
    if (stat != NC_NOERR) {
        printf("Error at line %d: %s\n", line, nc_strerror(stat));
        exit(1);
    }
}

int main(int argc, char **argv) {
    int ncid, xdim, ydim, tdim;
    int xid, yid, tid, runoff_id;

    // Create the file
    check_err(nc_create(FILE_NAME, NC_CLOBBER | NC_64BIT_OFFSET, &ncid), __LINE__);

    // Define dimensions
    check_err(nc_def_dim(ncid, "xaxis", NX, &xdim), __LINE__);
    check_err(nc_def_dim(ncid, "yaxis", NY, &ydim), __LINE__);
    check_err(nc_def_dim(ncid, "Time", NC_UNLIMITED, &tdim), __LINE__);

    // Define variables
    check_err(nc_def_var(ncid, "xaxis", NC_DOUBLE, 1, &xdim, &xid), __LINE__);
    check_err(nc_put_att_text(ncid, xid, "units", 12, "degrees east"), __LINE__);
    check_err(nc_put_att_text(ncid, xid, "long_name", 9, "longitude"), __LINE__);
    check_err(nc_put_att_text(ncid, xid, "cartesian_axis", 1, "X"), __LINE__);

    check_err(nc_def_var(ncid, "yaxis", NC_DOUBLE, 1, &ydim, &yid), __LINE__);
    check_err(nc_put_att_text(ncid, yid, "units", 13, "degrees north"), __LINE__);
    check_err(nc_put_att_text(ncid, yid, "long_name", 7, "latitude"), __LINE__);
    check_err(nc_put_att_text(ncid, yid, "cartesian_axis", 1, "Y"), __LINE__);

    check_err(nc_def_var(ncid, "Time", NC_DOUBLE, 1, &tdim, &tid), __LINE__);
    check_err(nc_put_att_text(ncid, tid, "units", 10, "time level"), __LINE__);
    check_err(nc_put_att_text(ncid, tid, "long_name", 4, "Time"), __LINE__);
    check_err(nc_put_att_text(ncid, tid, "cartesian_axis", 1, "T"), __LINE__);

    check_err(nc_def_var(ncid, "runoff", NC_DOUBLE, 3, (int[]){tdim, ydim, xdim}, &runoff_id), __LINE__);
    check_err(nc_put_att_text(ncid, runoff_id, "units", 8, "kg/2/m^2"), __LINE__);
    check_err(nc_put_att_text(ncid, runoff_id, "long_name", 12, "river runoff"), __LINE__);

    // Add Global attribute
    check_err(nc_put_att_text(ncid, NC_GLOBAL, "description", 15, "Test input file"), __LINE__);

    // End define mode
    check_err(nc_enddef(ncid), __LINE__);

    // Allocate memory for sample data and write to file
    double xdata[NX], ydata[NY];

    for (int i = 0; i < NX; i++)
        xdata[i] = (float)i - 179.5;
    for (int i = 0; i < NY; i++)
        ydata[i] = (float)i - 89.5;

    check_err(nc_put_var_double(ncid, xid, &xdata[0]), __LINE__);
    check_err(nc_put_var_double(ncid, yid, &ydata[0]), __LINE__);

    double *runoff_data = malloc(NY*NX*sizeof(double));

    double time_data;
    int i;
    for (int rec = 0; rec < TIME; rec++)
    {
        for (int y = 0; y < NY; y++) {
            double yval = (double)NY_MIN + (double)(NY_MAX - NY_MIN) * y / (double)NY;
            for (int x = 0; x < NX; x++)
            {
                double xval = (double)NX_MIN + (double)(NX_MAX - NX_MIN) * x / (double)NX;
                i = y*NX + x;
                runoff_data[i] = sin(sqrt(xval*xval+yval*yval))/sqrt(xval*xval+yval*yval);
            }
        }
        time_data = (float)rec;
        check_err(nc_put_vara_double(ncid, tid,
                                    (size_t[1]){rec},
                                    (size_t[1]){1},
                                    &time_data), __LINE__);
        check_err(nc_put_vara_double(ncid, runoff_id,
                                    (size_t[3]){rec, 0, 0},
                                    (size_t[4]){1, NY, NX},
                                    &runoff_data[0]), __LINE__);
    }
    // Close the file
    check_err(nc_close(ncid), __LINE__);
    free(runoff_data);
    printf("Netcdf file %s created successfully!\n", FILE_NAME);
    return 0;
}
