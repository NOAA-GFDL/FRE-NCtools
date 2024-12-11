#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define FILE_NAME "ocean_mask_grid.nc"
#define NX 76
#define NY 82

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
    int ncid, xt_dim, yt_dim;
    int wet_id;

    // Create the file
    check_err(nc_create(FILE_NAME, NC_CLOBBER, &ncid), __LINE__);

    // Define dimensions
    check_err(nc_def_dim(ncid, "grid_xt", NX, &xt_dim), __LINE__);
    check_err(nc_def_dim(ncid, "grid_yt", NY, &yt_dim), __LINE__);

    // Define variables
    check_err(nc_def_var(ncid, "wet", NC_DOUBLE, 2, (int[]){yt_dim, xt_dim}, &wet_id), __LINE__);

    check_err(nc_put_att_text(ncid, NC_GLOBAL, "y_boundary_type", 11, "solid_walls"), __LINE__);
    check_err(nc_put_att_text(ncid, NC_GLOBAL, "x_boundary_type", 11, "solid_walls"), __LINE__);

    // End define mode
    check_err(nc_enddef(ncid), __LINE__);

    // Allocate memory for sample data and write to file
    double wet_data[NY][NX];

    float xval, yval;
    for (int x = 0; x < NX; x++) {
        xval = (NX_MIN + (double)(NX_MAX - NX_MIN) * x / (NX - 1));
        for (int y = 0; y < NY; y++) {
            yval = (NY_MIN + (double)(NY_MAX - NY_MIN) * y / (NY - 1));
            wet_data[y][x] = 76.0 * (xval * xval + yval * yval) - 100.0;
            if (wet_data[y][x] < 0.0) wet_data[y][x] = 0.0;
        }
    }

    check_err(nc_put_var_double(ncid, wet_id, &wet_data[0][0]), __LINE__);

    // Close the file
    check_err(nc_close(ncid), __LINE__);
}
