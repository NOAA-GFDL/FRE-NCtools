#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define FILE_NAME_TEMPLATE "test_input_%d.nc"
#define NUM_FILES 5
#define NX 5
#define NY 5
#define NBND 2
#define TIME 12
#define MAX_FNAME 256

void check_err(int stat, int line) {
    if (stat != NC_NOERR) {
        printf("Error at line %d: %s\n", line, nc_strerror(stat));
        exit(1);
    }
}

int main(int argc, char **argv) {
    int ncid, time_dim, xt_dim, yt_dim, bnds_dim;
    int times_id, xt_id, yt_id, time_bnds_id, var1_id, var2_id;
    char file_name[MAX_FNAME];

    int nc_format = NC_CLOBBER;

    size_t start_var[3] = {0, 0, 0};
    size_t count_var[3] = {1, NY, NX};

    // Allocate memory for sample data and write to file
    float xt_data[NX], yt_data[NY], bnds_data[NBND];
    float var1_data[NY][NX], var2_data[NY][NX];

    // Init the static data
    for (int i = 0; i < NX; i++) xt_data[i] = i;
    for (int i = 0; i < NY; i++) yt_data[i] = i;

    for (int i=0; i<NUM_FILES; i++) {
        sprintf(file_name, FILE_NAME_TEMPLATE, i+1);

        // Create the file
        check_err(nc_create(file_name, NC_CLOBBER, &ncid), __LINE__);

        // Define dimensions
        check_err(nc_def_dim(ncid, "time", NC_UNLIMITED, &time_dim), __LINE__);
        check_err(nc_def_dim(ncid, "bnds", NBND, &bnds_dim), __LINE__);
        check_err(nc_def_dim(ncid, "lat", NX, &xt_dim), __LINE__);
        check_err(nc_def_dim(ncid, "lon", NY, &yt_dim), __LINE__);

        // Define variables
        check_err(nc_def_var(ncid, "time", NC_DOUBLE, 1, &time_dim, &times_id), __LINE__);
        check_err(nc_put_att_text(ncid, times_id, "bounds", 10, "time_bnds"), __LINE__);

        check_err(nc_def_var(ncid, "lat", NC_FLOAT, 1, &yt_dim, &xt_id), __LINE__);
        check_err(nc_def_var(ncid, "lon", NC_FLOAT, 1, &xt_dim, &yt_id), __LINE__);
        check_err(nc_def_var(ncid, "time_bnds", NC_FLOAT, 2, (int[]){time_dim, bnds_dim}, &time_bnds_id), __LINE__);

        check_err(nc_def_var(ncid, "var1", NC_FLOAT, 3, (int[]){time_dim, yt_dim, xt_dim}, &var1_id), __LINE__);
        check_err(nc_put_att_text(ncid, var1_id, "type", 11, "times data"), __LINE__);

        check_err(nc_def_var(ncid, "var2", NC_FLOAT, 3, (int[]){time_dim, yt_dim, xt_dim}, &var2_id), __LINE__);
        check_err(nc_put_att_text(ncid, var2_id, "type", 12, "random data"), __LINE__);

        // End define mode
        check_err(nc_enddef(ncid), __LINE__);

        // Set the axis data
        check_err(nc_put_var_float(ncid, xt_id, &xt_data[0]), __LINE__);
        check_err(nc_put_var_float(ncid, yt_id, &yt_data[0]), __LINE__);

        // Allocate memory for sample data and write to file
        double times_data = 0;
        double bnds_data[NBND] = {0, 0};

        for (int rec=0; rec < TIME; rec++)
        {
            times_data = (double)(i*TIME+rec+1);
            bnds_data[0] = times_data - 0.5;
            bnds_data[1] = times_data + 0.5;

            for (int y = 0; y < NY; y++)
                for (int x = 0; x < NX; x++)
                    var1_data[y][x] = times_data;

            for (int y = 0; y < NY; y++)
                for (int x = 0; x < NX; x++)
                    var2_data[y][x] = (rand() % 100);

            check_err(nc_put_vara_double(ncid, times_id,
                                        (size_t[1]){rec},
                                        (size_t[1]){1},
                                        &times_data),
                      __LINE__);
            check_err(nc_put_vara_double(ncid, time_bnds_id,
                                         (size_t[2]){rec, 0},
                                         (size_t[2]){1, NBND},
                                         &bnds_data[0]),
                      __LINE__);
            check_err(nc_put_vara_float(ncid, var1_id,
                                        (size_t[3]){rec, 0, 0},
                                        (size_t[3]){1, NY, NX},
                                        &var1_data[0][0]),
                      __LINE__);
            check_err(nc_put_vara_float(ncid, var2_id,
                                        (size_t[3]){rec, 0, 0},
                                        (size_t[3]){1, NY, NX},
                                        &var2_data[0][0]),
                      __LINE__);
        }
        // Close the file
        check_err(nc_close(ncid), __LINE__);

        printf("NetCDF file \"%s\" created successfully!\n", file_name);
    }
    return 0;
}
