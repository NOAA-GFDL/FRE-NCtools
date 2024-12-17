#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define FILE_NAME "ocean_temp_salt.nc"
#define NX 360
#define NY 200
#define NZ 50

#define TEMP_MIN -2
#define TEMP_MAX 32
#define SALT_MIN 0
#define SALT_MAX 41
#define TIME 1

void check_err(int stat, int line) {
    if (stat != NC_NOERR) {
        printf("Error at line %d: %s\n", line, nc_strerror(stat));
        exit(1);
    }
}

int main(int argc, char **argv) {
    int ncid, xdim, ydim, zdim, tdim;
    int xid, yid, zid, tid, temp_id, salt_id;

    // Create the file
    check_err(nc_create(FILE_NAME, NC_CLOBBER | NC_64BIT_OFFSET, &ncid), __LINE__);

    // Define dimensions
    check_err(nc_def_dim(ncid, "xaxis", NX, &xdim), __LINE__);
    check_err(nc_def_dim(ncid, "yaxis", NY, &ydim), __LINE__);
    check_err(nc_def_dim(ncid, "zaxis", NZ, &zdim), __LINE__);
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
    check_err(nc_put_att_text(ncid, zid, "units", 4, "none"), __LINE__);
    check_err(nc_put_att_text(ncid, zid, "long_name", 7, "zaxis_1"), __LINE__);
    check_err(nc_put_att_text(ncid, zid, "cartesian_axis", 1, "Z"), __LINE__);

    check_err(nc_def_var(ncid, "Time", NC_DOUBLE, 1, &tdim, &tid), __LINE__);
    check_err(nc_put_att_text(ncid, tid, "units", 10, "time level"), __LINE__);
    check_err(nc_put_att_text(ncid, tid, "long_name", 4, "Time"), __LINE__);
    check_err(nc_put_att_text(ncid, tid, "cartesian_axis", 1, "T"), __LINE__);

    check_err(nc_def_var(ncid, "temp", NC_DOUBLE, 4, (int[]){tdim, zdim, ydim, xdim}, &temp_id), __LINE__);
    check_err(nc_put_att_text(ncid, temp_id, "units", 4, "none"), __LINE__);
    check_err(nc_put_att_text(ncid, temp_id, "long_name", 4, "temp"), __LINE__);

    check_err(nc_def_var(ncid, "salt", NC_DOUBLE, 4, (int[]){tdim, zdim, ydim, xdim}, &salt_id), __LINE__);
    check_err(nc_put_att_text(ncid, salt_id, "units", 4, "none"), __LINE__);
    check_err(nc_put_att_text(ncid, salt_id, "long_name", 4, "salt"), __LINE__);

    // Add Global attribute
    check_err(nc_put_att_text(ncid, NC_GLOBAL, "description", 16, "Test input file"), __LINE__);

    // End define mode
    check_err(nc_enddef(ncid), __LINE__);

    // Allocate memory for sample data and write to file
    double xdata[NX], ydata[NY], zdata[NZ];

    for (int i = 0; i < NX; i++)
        xdata[i] = (i+1);
    for (int i = 0; i < NY; i++)
        ydata[i] = (i+1);
    for (int i = 0; i < NZ; i++)
        zdata[i] = (i+1);

    check_err(nc_put_var_double(ncid, xid, &xdata[0]), __LINE__);
    check_err(nc_put_var_double(ncid, yid, &ydata[0]), __LINE__);
    check_err(nc_put_var_double(ncid, zid, &zdata[0]), __LINE__);

    double *temp_data = malloc(NZ*NY*NX*sizeof(double));
    double *salt_data = malloc(NZ*NY*NX*sizeof(double));
    double time_data;
    int i;
    for (int rec = 0; rec < TIME; rec++)
    {
        for (int z = 0; z < NZ; z++)
            for (int y = 0; y < NY; y++)
                for (int x = 0; x < NX; x++)
                {
                    i = z*NY*NX + y*NX + x;
                    temp_data[i] = (TEMP_MAX-TEMP_MIN) *
                                   ((float)x/NX + (float)y/NY + (float)z/NZ + (float)rec/TIME);
                    salt_data[i] = (SALT_MAX-SALT_MIN)/NZ *
                                   ((float)x/NX + (float)y/NY + (float)z/NZ + (float)rec/TIME);
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
