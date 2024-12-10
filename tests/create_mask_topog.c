#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define FILE_NAME "mask_topog.nc"
#define NX 360
#define NY 200

#define NX_MIN -6
#define NX_MAX 6
#define NY_MIN -6
#define NY_MAX 6

void check_err(int stat, int line) {
    if (stat != NC_NOERR) {
        printf("Error at line %d: %s\n", line, nc_strerror(stat));
        exit(1);
    }
}

int main(int argc, char **argv) {
    int ncid, xt_dim, yt_dim, nt_dim;
    int depth_id;

    // Create the file
    check_err(nc_create(FILE_NAME, NC_CLOBBER, &ncid), __LINE__);

    // Define dimensions
    check_err(nc_def_dim(ncid, "grid_xt", NX, &xt_dim), __LINE__);
    check_err(nc_def_dim(ncid, "grid_yt", NY, &yt_dim), __LINE__);
    check_err(nc_def_dim(ncid, "ntiles", 1, &nt_dim), __LINE__);

    // Define variables
    check_err(nc_def_var(ncid, "depth", NC_DOUBLE, 2, (int[]){yt_dim, xt_dim}, &depth_id), __LINE__);

    // End define mode
    check_err(nc_enddef(ncid), __LINE__);

    // Allocate memory for sample data and write to file
    double depth_data[NY][NX];

    float xval, yval;
    for (int x = 0; x < NX; x++) {
        xval = (NX_MIN + (double)(NX_MAX - NX_MIN) * x / (NX - 1));
        for (int y = 0; y < NY; y++) {
            yval = (NY_MIN + (double)(NY_MAX - NY_MIN) * y / (NY - 1));
            depth_data[y][x] = 76.0 * (xval * xval + yval * yval) - 100.0;
            if (depth_data[y][x] < 0.0) depth_data[y][x] = 0.0;
        }
    }

    check_err(nc_put_var_double(ncid, depth_id, &depth_data[0][0]), __LINE__);

    // Close the file
    check_err(nc_close(ncid), __LINE__);
}
