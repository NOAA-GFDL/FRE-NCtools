#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <libgen.h>

#include <netcdf.h>

#define OUT_FILE_BASE "land_remap"
#define LAND_TILE 5
#define TILE_INDEX 5
#define TILE_INDEX_MAX LAND_TILE*LAT*LON

void check_err(int stat, int line)
{
    if (stat != NC_NOERR)
    {
        printf("Error at line %d: %s\n", line, nc_strerror(stat));
        exit(EXIT_FAILURE);
    }
}

void generate_unique_random_numbers(int *array, int size, int range) {
    int count = 0;
    while (count < size) {
        int num = rand() % range;
        int unique = 1;
        for (int i = 0; i < count; i++) {
            if (array[i] == num) {
                unique = 0;
                break;
            }
        }
        if (unique) {
            array[count] = num;
            count++;
        }
    }
}

int main(int argc, char **argv)
{
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

    // Loop over the tile files
    for (int i = 1; i <= 6; i++)
    {
        // Open the input tile file
        snprintf(in_tile_file_name, in_fn_length + 1, "%s.tile%d.nc", argv[1], i);
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
        size_t out_fn_length = snprintf(NULL, 0, "%s%zu.tile%d.nc", OUT_FILE_BASE, nx_in/2, i);
        char *out_file_name = malloc(out_fn_length + 1);
        snprintf(out_file_name, out_fn_length + 1, "%s%zu.tile%d.nc", OUT_FILE_BASE, nx_in/2, i);
        check_err(nc_create(out_file_name, NC_CLOBBER, &ncid_out), __LINE__);
        // Define the dimensions
        int lon_dimid, lat_dimid, tile_dimid, tile_index_dimid;
        check_err(nc_def_dim(ncid_out, "lon", nx_in/2, &lon_dimid), __LINE__);
        check_err(nc_def_dim(ncid_out, "lat", ny_in/2, &lat_dimid), __LINE__);
        check_err(nc_def_dim(ncid_out, "tile", LAND_TILE, &tile_dimid), __LINE__);
        check_err(nc_def_dim(ncid_out, "tile_index", TILE_INDEX, &tile_index_dimid), __LINE__);
        // Define dimension variables
        int lon_varid, lat_varid, tile_varid, tile_index_varid;
        check_err(nc_def_var(ncid_out, "lon", NC_DOUBLE, 1, &lon_dimid, &lon_varid), __LINE__);
        check_err(nc_put_att_text(ncid_out, lon_varid, "units", 12, "degrees east"), __LINE__);
        check_err(nc_put_att_text(ncid_out, lon_varid, "long_name", 9, "longitude"), __LINE__);
        check_err(nc_put_att_text(ncid_out, lon_varid, "cartesian_axis", 1, "X"), __LINE__);
        check_err(nc_def_var(ncid_out, "lat", NC_DOUBLE, 1, &lat_dimid, &lat_varid), __LINE__);
        check_err(nc_put_att_text(ncid_out, lon_varid, "units", 13, "degrees north"), __LINE__);
        check_err(nc_put_att_text(ncid_out, lon_varid, "long_name", 8, "latitude"), __LINE__);
        check_err(nc_put_att_text(ncid_out, lon_varid, "cartesian_axis", 1, "Y"), __LINE__);
        check_err(nc_def_var(ncid_out, "tile_index", NC_INT, 1, &tile_index_dimid, &tile_varid), __LINE__);
        check_err(nc_put_att_text(ncid_out, tile_varid, "long_name", 27, "compressed land point index"), __LINE__);
        check_err(nc_put_att_text(ncid_out, tile_varid, "compress", 12, "tile lat lon"), __LINE__);
        check_err(nc_put_att_float(ncid_out, tile_varid, "valid_min", NC_FLOAT, 1, (float[]){0.0}), __LINE__);
        // Define other variables -- Using land restart-like variables
        int frac_varid, soil_varid, glac_varid, lake_varid, vegn_varid;
        check_err(nc_def_var(ncid_out, "frac", NC_DOUBLE, 1, &tile_index_dimid, &frac_varid), __LINE__);
        check_err(nc_put_att_text(ncid_out, frac_varid, "long_name", 23,"fractional area of tile"), __LINE__);
        check_err(nc_def_var(ncid_out, "soil", NC_INT, 1, &tile_index_dimid, &soil_varid), __LINE__);
        check_err(nc_put_att_text(ncid_out, soil_varid, "long_name", 17, "tag of soil tiles"), __LINE__);
        check_err(nc_def_var(ncid_out, "glac", NC_INT, 1, &tile_index_dimid, &glac_varid), __LINE__);
        check_err(nc_put_att_text(ncid_out, glac_varid, "long_name", 17, "tag of glac tiles"), __LINE__);
        check_err(nc_def_var(ncid_out, "lake", NC_INT, 1, &tile_index_dimid, &lake_varid), __LINE__);
        check_err(nc_put_att_text(ncid_out, lake_varid, "long_name", 17, "tag of lake tiles"), __LINE__);
        check_err(nc_def_var(ncid_out, "vegn", NC_INT, 1, &tile_index_dimid, &vegn_varid), __LINE__);
        check_err(nc_put_att_text(ncid_out, vegn_varid, "long_name", 17, "tag of vegn tiles"), __LINE__);

        // End define mode of out file
        check_err(nc_enddef(ncid_out), __LINE__);

        // Write the static variable data
        double lon[nx_in/2], lat[ny_in/2];
        for (int xi = 0; xi < nx_in/2; xi++)
        {
            // We need to get the lon from the in file's x[1][1::2] (numpy format)
            lon[xi] = x_in[nxp_in + 2 * xi + 1];
        }
        for (int yi = 0; yi < ny_in/2; yi++)
        {
            // We need to get the lat from the in file's y[1::2][1] (numpy format)
            lat[yi] = y_in[(2*yi+1)*nxp_in + 1];
        }
        check_err(nc_put_var_double(ncid_out, lon_varid, lon), __LINE__);
        check_err(nc_put_var_double(ncid_out, lat_varid, lat), __LINE__);

        // Generate and write the tile index data
        int tile_index_data[TILE_INDEX];
        generate_unique_random_numbers(tile_index_data, TILE_INDEX, LAND_TILE*nx_in*ny_in/4);
        check_err(nc_put_var_int(ncid_out, tile_varid, tile_index_data), __LINE__);
        // Generate and write the frac data (using a random number generator)
        double frac_data[TILE_INDEX];
        for (int ti = 0; ti < TILE_INDEX; ti++)
            frac_data[ti] = (double)rand() / RAND_MAX;
        check_err(nc_put_var_double(ncid_out, frac_varid, frac_data), __LINE__);
        // Generate and write the soil data
        int soil_data[TILE_INDEX];
        for (int ti = 0; ti < TILE_INDEX; ti++)
            soil_data[ti] = (int)(rand() % LAND_TILE);
        check_err(nc_put_var_int(ncid_out, soil_varid, soil_data), __LINE__);
        // Generate and write the glac data
        int glac_data[TILE_INDEX];
        for (int ti = 0; ti < TILE_INDEX; ti++)
            if ((float)rand() / RAND_MAX < 0.5)
                glac_data[ti] = NC_FILL_INT;
            else
                glac_data[ti] = 1;
        check_err(nc_put_var_int(ncid_out, glac_varid, glac_data), __LINE__);
        // Generate and write the lake data
        int lake_data[TILE_INDEX];
        for (int ti = 0; ti < TILE_INDEX; ti++)
            if ((float)rand() / RAND_MAX < 0.5)
                lake_data[ti] = NC_FILL_INT;
            else
                lake_data[ti] = 1;
        check_err(nc_put_var_int(ncid_out, lake_varid, lake_data), __LINE__);
        // Generate and write the vegn data
        int vegn_data[TILE_INDEX];
        for (int ti = 0; ti < TILE_INDEX; ti++)
            if ((float)rand() / RAND_MAX < 0.5)
                vegn_data[ti] = NC_FILL_INT;
            else
                vegn_data[ti] = 1;
        check_err(nc_put_var_int(ncid_out, vegn_varid, vegn_data), __LINE__);

        // Close the output file
        check_err(nc_close(ncid_out), __LINE__);
        free(out_file_name);
        free(x_in);
        free(y_in);
    }
    free(in_tile_file_name);
    return EXIT_SUCCESS;
}
