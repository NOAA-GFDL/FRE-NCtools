#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <libgen.h>
#include <math.h>

#include <netcdf.h>

#define OUT_FILE_BASE "test_output"
#define TIME 10

void check_err(int stat, int line) {
    if (stat != NC_NOERR) {
        printf("Error at line %d: %s\n", line, nc_strerror(stat));
        exit(1);
    }
}

int main(int argc, char **argv) {
    int ncid_in;
    int ncid_out;

    // The user should pass in the input file name
    // without the .tile?.nc suffix
    if (argc != 2) {
        fprintf(stderr,
                "Usage: %s <input_file_name>\n",
                basename(argv[0]));
        exit(1);
    }
    // Where to store the filename
    size_t in_fn_length = snprintf(NULL, 0, "%s.tile1.nc", argv[1]);
    char *in_tile_file_name = malloc(in_fn_length + 1);
    size_t out_fn_length = snprintf(NULL, 0, "%s.tile1.nc", OUT_FILE_BASE);
    char *out_tile_file_name = malloc(out_fn_length + 1);
    // Check that all tile files exist
    for (int i = 1; i <= 6; i++) {
        snprintf(in_tile_file_name, in_fn_length + 1, "%s.tile%d.nc", argv[1], i);
        if (access(in_tile_file_name, F_OK) == -1) {
            fprintf(stderr,
                    "Error: %s does not exist\n",
                    in_tile_file_name);
            exit(1);
        }
    }
    for (int i = 1; i <= 6; i++) {
        snprintf(in_tile_file_name, in_fn_length + 1, "%s.tile%d.nc", argv[1], i);
        // Open the input file
        check_err(nc_open(in_tile_file_name, NC_NOWRITE, &ncid_in), __LINE__);
        // Get the dimensions from the input file.
        int nx_dimid, ny_dimid, nxp_dimid, nyp_dimid;
        size_t nx_in, ny_in, nxp_in, nyp_in;
        check_err(nc_inq_dimid(ncid_in, "nx", &nx_dimid), __LINE__);
        check_err(nc_inq_dimid(ncid_in, "ny", &ny_dimid), __LINE__);
        check_err(nc_inq_dimid(ncid_in, "nxp", &nxp_dimid), __LINE__);
        check_err(nc_inq_dimid(ncid_in, "nyp", &nyp_dimid), __LINE__);
        check_err(nc_inq_dimlen(ncid_in, nx_dimid, &nx_in), __LINE__);
        check_err(nc_inq_dimlen(ncid_in, ny_dimid, &ny_in), __LINE__);
        check_err(nc_inq_dimlen(ncid_in, nxp_dimid, &nxp_in), __LINE__);
        check_err(nc_inq_dimlen(ncid_in, nyp_dimid, &nyp_in), __LINE__);
        // Read in the x and y variables
        double *x_in = malloc(nxp_in * nyp_in * sizeof(double));
        double *y_in = malloc(nxp_in * nyp_in * sizeof(double));
        int x_varid, y_varid;
        check_err(nc_inq_varid(ncid_in, "x", &x_varid), __LINE__);
        check_err(nc_inq_varid(ncid_in, "y", &y_varid), __LINE__);
        check_err(nc_get_var_double(ncid_in, x_varid, x_in), __LINE__);
        check_err(nc_get_var_double(ncid_in, y_varid, y_in), __LINE__);
        // Close the input file
        check_err(nc_close(ncid_in), __LINE__);

        // Create the output file
        snprintf(out_tile_file_name, out_fn_length + 1, "%s.tile%d.nc", OUT_FILE_BASE, i);
        check_err(nc_create(out_tile_file_name, NC_CLOBBER, &ncid_out), __LINE__);
        // Define the dimensions
        int xt_dimid, yt_dimid;
        check_err(nc_def_dim(ncid_out, "grid_xt", nx_in/2, &xt_dimid), __LINE__);
        check_err(nc_def_dim(ncid_out, "grid_yt", ny_in/2, &yt_dimid), __LINE__);
        // Define the grid variables
        int xt_varid, yt_varid;
        check_err(nc_def_var(ncid_out, "grid_xt", NC_FLOAT, 1, &xt_dimid, &xt_varid), __LINE__);
        check_err(nc_put_att_text(ncid_out, xt_varid, "units", 13, "degrees east"), __LINE__);
        check_err(nc_put_att_text(ncid_out, xt_varid, "long_name", 16, "T-cell latitude"), __LINE__);
        check_err(nc_put_att_text(ncid_out, xt_varid, "cartesian_axis", 2, "X"), __LINE__);
        check_err(nc_def_var(ncid_out, "grid_yt", NC_FLOAT, 1, &yt_dimid, &yt_varid), __LINE__);
        check_err(nc_put_att_text(ncid_out, yt_varid, "units", 14, "degrees north"), __LINE__);
        check_err(nc_put_att_text(ncid_out, yt_varid, "long_name", 17, "T-cell longitude"), __LINE__);
        check_err(nc_put_att_text(ncid_out, yt_varid, "cartesian_axis", 2, "Y"), __LINE__);
        // Define a static variable
        int static_varid;
        check_err(nc_def_var(ncid_out, "static", NC_FLOAT, 2, (int[]){yt_dimid, xt_dimid}, &static_varid), __LINE__);

        // Define the time dimension
        int time_dimid;
        check_err(nc_def_dim(ncid_out, "time", TIME, &time_dimid), __LINE__);
        // Define the time variable
        int time_varid;
        check_err(nc_def_var(ncid_out, "time", NC_DOUBLE, 1, &time_dimid, &time_varid), __LINE__);
        check_err(nc_put_att_text(ncid_out, time_varid, "units", 31, "days since 0001-01-01 00:00:00"), __LINE__);
        check_err(nc_put_att_text(ncid_out, time_varid, "long_name", 5, "Time"), __LINE__);
        check_err(nc_put_att_text(ncid_out, time_varid, "calendar", 7, "noleap"), __LINE__);
        check_err(nc_put_att_text(ncid_out, time_varid, "calendar_type", 7, "noleap"), __LINE__);
        // Define a time-dependent variable
        int tvar_varid;
        check_err(nc_def_var(ncid_out, "tvar", NC_FLOAT, 3, (int[]){time_dimid, yt_dimid, xt_dimid}, &tvar_varid), __LINE__);

        // End define mode
        check_err(nc_enddef(ncid_out), __LINE__);

        // Write the static variables
        float xt[nx_in/2], yt[nx_in/2];
        for (int i = 0; i < nx_in/2; i++)
            xt[i] = (float)(i+1);
        for (int i = 0; i < ny_in/2; i++)
            yt[i] = (float)(i+1);
        check_err(nc_put_var_float(ncid_out, xt_varid, xt), __LINE__);
        check_err(nc_put_var_float(ncid_out, yt_varid, yt), __LINE__);

        float *var_data = malloc(nx_in/2 * ny_in/2 * sizeof(float));
        float xval, yval;
        for (int i = 0; i < ny_in/2; i++) {
            yval = y_in[i*2];
            for (int j = 0; j < nx_in/2; j++) {
                xval = x_in[j*2];
                var_data[i*ny_in/2 + j] = 10.0 * sin((xval + yval) * M_PI / 180.0);
            }
        }
        check_err(nc_put_var_float(ncid_out, static_varid, var_data), __LINE__);

        // Write the time-dependent variables
        for (int rec=0; rec < TIME; rec++) {
            for (int i = 0; i < ny_in/2; i++) {
                yval = y_in[i*2];
                for (int j = 0; j < nx_in/2; j++) {
                    xval = x_in[j*2];
                    var_data[i*ny_in/2 + j] = 10.0 * sin(rec + (xval + yval) * M_PI / 180.0);
                }
            }
            check_err(nc_put_vara_float(ncid_out, tvar_varid,
                                        (size_t[3]){rec, 0, 0},
                                        (size_t[3]){1, ny_in/2, nx_in/2},
                                        var_data), __LINE__);
        }
        // Close the output file
        check_err(nc_close(ncid_out), __LINE__);
        free(var_data);
    }
    free (in_tile_file_name);
    free (out_tile_file_name);
    return 0;
}
