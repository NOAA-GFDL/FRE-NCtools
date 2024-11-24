#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define FILE_NAME "test_input.nc"
#define NX 10
#define NY 10
#define PFULL 33
#define PHALF 34
#define TIME 10
#define PS_MIN 75223
#define PS_MAX 101325
#define TEMP_MIN 250
#define TEMP_MAX 330

void check_err(int stat, int line) {
    if (stat != NC_NOERR) {
        printf("Error at line %d: %s\n", line, nc_strerror(stat));
        exit(1);
    }
}

int main(int argc, char **argv) {
    int ncid, time_dim, xt_dim, yt_dim, pfull_dim, phalf_dim;
    int times_id, xt_id, yt_id, phalf_id, pfull_id, pk_id, bk_id;
    int ps_id, temp_id, dummy_id;

    int nc_format = NC_CLOBBER;

    int stat;
    int c;

    while (( c = getopt (argc, argv, "abc:")) != -1)
    {
        switch (c)
        {
            case '4':
                nc_format = NC_CLOBBER | NC_NETCDF4;
                break;
            case '6':
                nc_format = NC_CLOBBER | NC_64BIT_OFFSET;
                break;
        }
    }

    size_t start_ps[3] = {0, 0, 0};
    size_t start_temp[4] = {0, 0, 0, 0};
    size_t count_ps[3] = {1, NY, NX};
    size_t count_temp[4] = {1, PFULL, NX, NY};

    // Create the file
    stat = nc_create(FILE_NAME, NC_CLOBBER, &ncid);
    check_err(stat, __LINE__);

    // Define dimensions
    stat = nc_def_dim(ncid, "time", NC_UNLIMITED, &time_dim);
    check_err(stat, __LINE__);
    stat = nc_def_dim(ncid, "grid_xt", NX, &xt_dim);
    check_err(stat, __LINE__);
    stat = nc_def_dim(ncid, "grid_yt", NY, &yt_dim);
    check_err(stat, __LINE__);
    stat = nc_def_dim(ncid, "pfull", PFULL, &pfull_dim);
    check_err(stat, __LINE__);
    stat = nc_def_dim(ncid, "phalf", PHALF, &phalf_dim);
    check_err(stat, __LINE__);

    // Define variables
    stat = nc_def_var(ncid, "time", NC_DOUBLE, 1, &time_dim, &times_id);
    check_err(stat, __LINE__);
    stat = nc_def_var(ncid, "grid_xt", NC_FLOAT, 1, &xt_dim, &xt_id);
    check_err(stat, __LINE__);
    stat = nc_def_var(ncid, "grid_yt", NC_FLOAT, 1, &yt_dim, &yt_id);
    check_err(stat, __LINE__);
    stat = nc_def_var(ncid, "phalf", NC_FLOAT, 1, &phalf_dim, &phalf_id);
    check_err(stat, __LINE__);
    stat = nc_def_var(ncid, "pfull", NC_FLOAT, 1, &pfull_dim, &pfull_id);
    check_err(stat, __LINE__);
    stat = nc_def_var(ncid, "pk", NC_FLOAT, 1, &phalf_dim, &pk_id);
    check_err(stat, __LINE__);
    stat = nc_def_var(ncid, "bk", NC_FLOAT, 1, &phalf_dim, &bk_id);
    check_err(stat, __LINE__);
    stat = nc_def_var(ncid, "ps", NC_FLOAT, 3, (int[]){time_dim, xt_dim, yt_dim}, &ps_id);
    check_err(stat, __LINE__);
    stat = nc_def_var(ncid, "temp", NC_FLOAT, 4, (int[]){time_dim, pfull_dim, xt_dim, yt_dim}, &temp_id);
    check_err(stat, __LINE__);
    stat = nc_def_var(ncid, "dummy", NC_FLOAT, 3, (int[]){time_dim, yt_dim, xt_dim}, &dummy_id);

    // Assign global attributes
    stat = nc_put_att_text(ncid, NC_GLOBAL, "description", 55, "Sample netCDF file for testing run_timepressure_interp");

    // Assign attributes to variables (for brevity, a few key attributes shown)
    stat = nc_put_att_text(ncid, times_id, "units", 30, "days since 0001-01-01 00:00:00");
    check_err(stat, __LINE__);
    stat = nc_put_att_text(ncid, xt_id, "units", 13, "degrees east");
    check_err(stat, __LINE__);
    stat = nc_put_att_text(ncid, yt_id, "units", 14, "degrees north");
    check_err(stat, __LINE__);

    // End define mode
    stat = nc_enddef(ncid);
    check_err(stat, __LINE__);

    // Allocate memory for sample data and write to file
    double times_data = 0;
    float xt_data[NX], yt_data[NY], temp_data[PFULL][NX][NY];
    float ps_data[NX][NY];
    float dummy_data[NX][NY];

    for (int i = 0; i < NX; i++) xt_data[i] = i;
    for (int i = 0; i < NY; i++) yt_data[i] = i;

    float pfull_data[PFULL] = {
        2.16404256133345, 5.84530754043837, 10.7450801553507,
        17.1065372595328, 25.1138051289766, 35.2211968195201,
        48.1379036870309, 64.560183540236, 85.1144822187467,
        110.419626793719, 141.092610422464, 177.729387604925,
        220.89239709789, 271.066623628003, 328.516336882098,
        392.785272564752, 461.94726249906, 532.465906646758,
        600.430867020261, 663.107382784874, 719.307118080772,
        768.814283964814, 811.846868750618, 848.836020639326,
        880.346139043136, 906.995722281721, 929.359836869601,
        947.9591141953, 963.258551283114, 975.628859671795,
        985.399151129319, 992.786752544069, 997.948596287477};
    float phalf_data[PHALF] = {
        1.0, 4.0, 8.1860211, 13.7888653, 20.9179519, 29.8364084,
        41.217896, 55.7922148, 74.2019063, 97.0478639,
        124.9666476, 158.5495525, 198.3969594, 245.0272213,
        298.8885756, 360.0401788, 427.4580251, 498.2435727,
        568.2205347, 633.8360465, 693.2663294, 745.9919863,
        792.0973733, 831.9219455, 865.9778137, 894.8725252,
        919.2279199, 939.5659324, 956.4021322, 970.1476613,
        981.1306646, 989.68, 995.9, 1000.0};
    float bk_data[PHALF] = {
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00513, 0.01969,
        0.04299, 0.07477, 0.11508, 0.16408, 0.22198, 0.28865,
        0.36281, 0.44112, 0.51882, 0.59185, 0.6581, 0.71694,
        0.76843, 0.81293, 0.851, 0.88331, 0.91055, 0.93331,
        0.95214, 0.9675, 0.97968, 0.98908, 0.99575, 1};
    float pk_data[PHALF] = {
        100.0, 400.0, 818.6021, 1378.886, 2091.795, 2983.641,
        4121.79, 5579.222, 6907.19, 7735.787, 8197.665, 8377.955,
        8331.696, 8094.722, 7690.857, 7139.018, 6464.803, 5712.357,
        4940.054, 4198.604, 3516.633, 2905.199, 2366.737, 1899.195,
        1497.781, 1156.253, 867.792, 625.5933, 426.2132, 264.7661,
        145.0665, 60.0, 15.0, 0.0};

    // Writing static data
    stat = nc_put_var_float(ncid, xt_id, &xt_data[0]);
    check_err(stat, __LINE__);
    stat = nc_put_var_float(ncid, yt_id, &yt_data[0]);
    check_err(stat, __LINE__);
    stat = nc_put_var_float(ncid, phalf_id, &phalf_data[0]);
    check_err(stat, __LINE__);
    stat=nc_put_var_float(ncid, pfull_id, &pfull_data[0]);
    check_err(stat, __LINE__);
    stat = nc_put_var_float(ncid, pk_id, &pk_data[0]);
    check_err(stat, __LINE__);
    stat = nc_put_var_float(ncid, bk_id, &bk_data[0]);
    check_err(stat, __LINE__);

    // Write time-variable data
    for (int rec=0; rec < TIME; rec++)
    {
        // Generate random data as in the Python example
        for (int i = 0; i < PFULL; i++)
            for (int j = 0; j < NX; j++)
                for (int k = 0; k < NY; k++)
                    temp_data[i][j][k] = (rand() % (TEMP_MAX - TEMP_MIN)) + TEMP_MIN;

        for (int i = 0; i < NX; i++)
            for (int j = 0; j < NY; j++)
                ps_data[i][j] = (rand() % (PS_MAX-PS_MIN)) + PS_MIN;

        for (int i = 0; i < NX; i++)
            for (int j = 0; j < NY; j++)
                dummy_data[i][j] = i * NX + j;

        times_data = (double)rec;

        stat = nc_put_vara_double(ncid, times_id,
                                  (size_t[1]){rec},
                                  (size_t[1]){1},
                                  &times_data);
        check_err(stat, __LINE__);
        stat = nc_put_vara_float(ncid, temp_id,
                                 (size_t[4]){rec, 0, 0, 0},
                                 (size_t[4]){1, PFULL, NX, NY},
                                 &temp_data[0][0][0]);
        check_err(stat, __LINE__);
        stat = nc_put_vara_float(ncid, ps_id,
                             (size_t[3]){rec, 0, 0},
                             (size_t[3]){1, NX, NY},
                             &ps_data[0][0]);
        check_err(stat, __LINE__);
        stat = nc_put_vara_float(ncid, dummy_id,
                                 (size_t[3]){rec, 0, 0},
                                 (size_t[3]){1, NX, NY},
                                 &dummy_data[0][0]);
        check_err(stat, __LINE__);
    }

    // Close the file
    stat = nc_close(ncid);
    check_err(stat, __LINE__);

    printf("NetCDF file created successfully!\n");
    return 0;
}
