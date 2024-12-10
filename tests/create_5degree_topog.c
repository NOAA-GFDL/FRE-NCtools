#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#define FILE_NAME "5degree_topog.nc"
#define NX 719
#define NY 359
#define NX_MIN 0
#define NY_MIN -90
#define X_min -10.0
#define X_max 10.0
#define Y_min -5.0
#define Y_max 5.0
#define TOPO_MIN -5800
#define TOPO_MAX 3450

void check_err(int stat, int line) {
    if (stat != NC_NOERR) {
        printf("Error at line %d: %s\n", line, nc_strerror(stat));
        exit(1);
    }
}

int main(int argc, char **argv) {
    int ncid, xt_dim, yt_dim;
    int xt_id, yt_id, topo_id;

    // Create the file
    check_err(nc_create(FILE_NAME, NC_CLOBBER, &ncid), __LINE__);

    // Define dimensions
    check_err(nc_def_dim(ncid, "grid_xt", NX, &xt_dim), __LINE__);
    check_err(nc_def_dim(ncid, "grid_yt", NY, &yt_dim), __LINE__);

    // Define variables
    check_err(nc_def_var(ncid, "grid_xt", NC_FLOAT, 1, &xt_dim, &xt_id), __LINE__);
        check_err(nc_put_att_text(ncid, xt_id, "units", 13, "degrees east"), __LINE__);
    check_err(nc_put_att_text(ncid, xt_id, "point_spacing", 5, "even"), __LINE__);

    check_err(nc_def_var(ncid, "grid_yt", NC_FLOAT, 1, &yt_dim, &yt_id), __LINE__);
    check_err(nc_put_att_text(ncid, yt_id, "units", 14, "degrees north"), __LINE__);
    check_err(nc_put_att_text(ncid, yt_id, "point_spacing", 5, "even"), __LINE__);

    check_err(nc_def_var(ncid, "TOPO", NC_FLOAT, 2, (int[]){yt_dim, xt_dim}, &topo_id), __LINE__);
    check_err(nc_put_att_text(ncid, topo_id, "units", 1, "m"), __LINE__);
    check_err(nc_put_att_text(ncid, topo_id, "long_name", 9, "Topography"), __LINE__);
    check_err(nc_put_att_text(ncid, topo_id, "long_name_mod", 30, "regrid: G_P5DEGREE on X@AAV, on Y@AAV"), __LINE__);
    check_err(nc_put_att_float(ncid, topo_id, "missing_value", NC_FLOAT, 1, (float[]){0.f}), __LINE__);
    check_err(nc_put_att_float(ncid, topo_id, "_FillValue", NC_FLOAT, 1, (float[]){0.f}), __LINE__);

    // End define mode
    check_err(nc_enddef(ncid), __LINE__);

    // Allocate memory for sample data and write to file
    float xt_data[NX], yt_data[NY], pt_delta;

    pt_delta = 360.0 / (NX + 1);
    for (int i = 0; i < NX; i++)
        xt_data[i] = (i+1) * pt_delta + NX_MIN;

    pt_delta = 180.0 / (NY+1);
    for (int i = 0; i < NY; i++)
       yt_data[i] = (i+1) * pt_delta + NY_MIN;

    check_err(nc_put_var_float(ncid, xt_id, &xt_data[0]), __LINE__);
    check_err(nc_put_var_float(ncid, yt_id, &yt_data[0]), __LINE__);

    float topo_data[NY][NX];
    float xval, yval;
    for (int i = 0; i < NY; i++)
    {
        xval = (X_min + (double)(X_max - X_min) * i / (NX));
        for (int j = 0; j < NX; j++)
        {
            yval = (Y_min + (double)(Y_max - Y_min) * j / (NY));
            topo_data[i][j] = (float)(TOPO_MAX - TOPO_MIN)/2.0 * sin(xval + yval) + (float)(TOPO_MAX + TOPO_MIN);
        }
    }
    check_err(nc_put_var_float(ncid, topo_id, &topo_data[0][0]), __LINE__);

    // Close the file
    check_err(nc_close(ncid), __LINE__);
}
