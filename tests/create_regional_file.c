#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>

#define FILE_NAME "rregion_input_file"

#define DEFAULT_GRID_SIZE 48
#define DEFAULT_TILE_NUMBER 3
//#define NY 10
//#define PFULL 33
//#define PHALF 34
#define TIME 10
//#define PS_MIN 75223
//#define PS_MAX 101325
//#define TEMP_MIN 250
//#define TEMP_MAX 330

void check_err(int stat, int line) {
    if (stat != NC_NOERR) {
        printf("Error at line %d: %s\n", line, nc_strerror(stat));
        exit(1);
    }
}

int main(int argc, char **argv) {
    int grid_size = DEFAULT_GRID_SIZE;
    int tile_number = DEFAULT_TILE_NUMBER;
    int ncid, time_dim, xt_dim, yt_dim;
    int time_id, xt_id, yt_id;
    int stat;
    int c;

    while (( c = getopt (argc, argv, "g:t:") ) != -1)
    {
        switch (c)
        {
            case 'g':
                grid_size = atoi(optarg);
                if (grid_size < 1)
                {
                    fprintf(stderr, "Gridsize must be positive\n");
                    exit(EXIT_FAILURE);
                }
                break;
            case 't':
                tile_number = atoi(optarg);
                if (tile_number > 6 || tile_number < 1)
                {
                    fprintf(stderr, "Tile number must be in the range [1,6]\n");
                    exit(EXIT_FAILURE);
                }
                break;
            case '?':
            default:
                exit(EXIT_FAILURE);
                break;
        }
    }
    // Seed the random number generator
    srand(time(NULL));

    size_t fn_length = snprintf(NULL, 0, "%s.tile1.nc", FILE_NAME);
    char *file_name = malloc(fn_length + 1);
    snprintf(file_name, fn_length + 1, "%s.tile%d.nc", FILE_NAME, tile_number);
    // Create the file
    check_err(nc_create(file_name, NC_CLOBBER, &ncid), __LINE__);

    // Set the regional dimensions
    int nx = 1 + (int)rand() % grid_size;
    int ny = 1 + (int)rand() % grid_size;
    int nx_start = 1 + (int)rand() % (grid_size - nx + 1);
    int ny_start = 1 + (int)rand() % (grid_size - ny + 1);

    // Define dimensions
    check_err(nc_def_dim(ncid, "time", NC_UNLIMITED, &time_dim), __LINE__);
    check_err(nc_def_dim(ncid, "grid_xt_sub01", nx, &xt_dim), __LINE__);
    check_err(nc_def_dim(ncid, "grid_yt_sub01", ny, &yt_dim), __LINE__);

    // Define variables
    check_err(nc_def_var(ncid, "time", NC_DOUBLE, 1, &time_dim, &time_id), __LINE__);
    check_err(nc_put_att_text(ncid, time_id, "units", 30, "days since 0001-01-01 00:00:00"), __LINE__);
    check_err(nc_put_att_text(ncid, time_id, "long_name", 4, "Time"), __LINE__);
    check_err(nc_put_att_text(ncid, time_id, "calendar", 6, "noleap"), __LINE__);
    check_err(nc_put_att_text(ncid, time_id, "calendar_type", 6, "noleap"), __LINE__);
    check_err(nc_put_att_text(ncid, time_id, "cartesian_axis", 1, "T"), __LINE__);

    check_err(nc_def_var(ncid, "grid_xt_sub01", NC_FLOAT, 1, &xt_dim, &xt_id), __LINE__);
    check_err(nc_put_att_text(ncid, xt_id, "units", 12, "degrees east"), __LINE__);
    check_err(nc_put_att_text(ncid, xt_id, "long_name", 15, "T-cell latitude"), __LINE__);
    check_err(nc_put_att_text(ncid, xt_id, "cartesian_axis", 1, "X"), __LINE__);

    check_err(nc_def_var(ncid, "grid_yt_sub01", NC_FLOAT, 1, &yt_dim, &yt_id), __LINE__);
    check_err(nc_put_att_text(ncid, yt_id, "units", 13, "degrees north"), __LINE__);
    check_err(nc_put_att_text(ncid, yt_id, "long_name", 16, "T-cell longitude"), __LINE__);
    check_err(nc_put_att_text(ncid, yt_id, "cartesian_axis", 1, "Y"), __LINE__);

    // Assign global attributes
    check_err(nc_put_att_text(ncid, NC_GLOBAL, "description", 54, "Sample netCDF file for testing run_timepressure_interp"), __LINE__);

    // End define mode
    check_err(nc_enddef(ncid), __LINE__);

    // Allocate memory for sample data and write to file
    double times_data = 0;
    float *xt_data = malloc(nx * sizeof(int));
    float *yt_data = malloc(ny * sizeof(int));

    for (int i = 0; i < nx; i++) xt_data[i] = 1 + nx_start + i;
    for (int i = 0; i < ny; i++) yt_data[i] = 1 + ny_start + i;

    // Writing static data
    stat = nc_put_var_float(ncid, xt_id, &xt_data[0]);
    check_err(stat, __LINE__);
    stat = nc_put_var_float(ncid, yt_id, &yt_data[0]);
    check_err(stat, __LINE__);

    // Write time-variable data
    for (int rec=0; rec < TIME; rec++)
    {
        times_data = (double)rec;

        stat = nc_put_vara_double(ncid, time_id,
                                  (size_t[1]){rec},
                                  (size_t[1]){1},
                                  &times_data);
    }

    // Close the file
    stat = nc_close(ncid);
    check_err(stat, __LINE__);

    printf("NetCDF file \"%s\" created successfully!\n", file_name);
    free(file_name);
    return EXIT_SUCCESS;
}
