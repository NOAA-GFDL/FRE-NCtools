#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <stdbool.h>
#include "read_mosaic.h"
#include "mosaic_util.h"
#include "tool_util.h"
#include "constant.h"
#include "mpp.h"
#include "mpp_domain.h"
#include "mpp_io.h"
#define  TILE_INDEX_NAME "tile_index"
#define  COHORT_INDEX_NAME "cohort_index"
#define  COHORT_NAME "cohort"
#define  NSPECIES_NAME "nspecies"
#define  TEXTLEN_NAME "textlen"
#define  SC_COHORT_NAME "soilCCohort"
#define  LC_COHORT_NAME "litterCCohort"
#define  SPECIES_NAMES "species_names"
#define  MDF_DAY_NAME "mdf_day"
#define  TILE_NAME "tile"
#define  LON_NAME  "lon"
#define  LAT_NAME  "lat"
#define  FRAC_NAME "frac"
#define  SOIL_NAME "soil"
#define  VEGN_NAME "vegn"
#define  GLAC_NAME "glac"
#define  LAKE_NAME "lake"
#define  LEVEL_NAME "zfull"
#define  N_ACCUM_NAME "n_accum"
#define  NMN_ACM_NAME "nmn_acm"
#define  D2R (M_PI/180.)
#define  R2D (180./M_PI)

char *usage[] = {
  "",
  "                                                                                ",
  "                                                                                ",
  "                    Usage of remap_land                                         ",
  "                                                                                ",
  "   remap_land --src_mosaic src_mosaic --src_restart src_restart                 ",
  "              --dst_mosaic dst_msoaic --dst_restart dst_restart                 ",
  "              --dst_cold_restart dst_cold_restart                               ",
  "              --land_src_restart land_src_restart                               ",
  "              --land_dst_cold_restart land_dst_cold_restart                     ",
  "              [--remap_file remap_file] [--print_memory]                        ",
  "                                                                                ",
          " remap_land remap land restart file from one mosaic grid to another mosaic grid ",
  " remap_land takes the following flag  s,                                          ",
  "                                                                                  ",
  " REQUIRED:                                                                        ",
  "                                                                                  ",
  " --src_mosaic src_mosaic      specify the source mosaic information. This file    ",
  "                              contains list of tile files which specify the grid  ",
  "                              information for each tile.                          ",
  "                                                                                  ",
  " --dst_mosaic dst_mosaic      specify the destination mosaic information. This    ",
  "                              file contains list of tile files which specify the  ",
  "                              grid information for each tile.                     ",
  "                                                                                  ",
  " --src_restart src_restart    specify the source restart file.                    ",
  "                                                                                  ",
  " --dst_restart dst_restart    specify the restart file to be generated on         ",
  "                              destination grid.                                   ",
  "                                                                                  ",
  " --dst_cold_restart file      specify the cold restart file destination grid.     ",
  "                              This is the input file. The dst_cold_restart_file   ",
  "                              could be obtained by running the experiment         ",
  "                              for 1 day with --dst_mosaic using cold restart.     ",
  "                                                                                  ",
  " --file_type file_type        spefify file type. Its value could be 'land',       ",
  "                              'cana', 'snow', 'glac', 'lake', 'soil' or 'vegn'.   ",
  "                              when file_type is 'cana' or 'snow', need to         ",
  "                              specify --land_src_restart and --land_cold_restart  ",
  "                                                                                  ",
  " OPTIONAL FLAGS                                                                   ",
  "                                                                                  ",
  "                                                                                  ",
  " --land_src_restart file      specify the source file of land.res.nc. It is       ",
  "                              required when the restart file is snow.res or       ",
  "                              cana.res                                            ",
  "                                                                                  ",
  " --land_dst_cold_restart file specify the destination file of land.res.nc. It is  ",
  "                              required when the restart file is snow.res or       ",
  "                              cana.res                                            ",
  "                                                                                  ",
  " --remap_file remap_file    specify the file name that saves remapping            ",
  "                            information. If remap_file is specified and the       ",
  "                            file does not exist, remapping information will be    ",
  "                            calculated and stored in remap_file. If remap_file    ",
  "                            is specified and the file exists, remapping           ",
  "                            information will be read from remap_file.             ",
  "                                                                                  ",
  " --print_memory             debug memory usage when it is set                     ",
  "                                                                                  ",
  "                                                                                  ",
  " Example: Remap land restart files from C192 mosaic onto C384 mosaic.             ",
  "                                                                                  ",
  " 1. Run model with C384 grid using cold restart for land. This generates          ",
  "    dst_cold_restart files.                                                       ",
  " 2. Create remapping file and remap land.res file.                                ",
  "    >aprun -n 64 remap_land_parallel --file_type land --src_mosaic C192_mosaic.nc ",
  "          --dst_mosaic c384_mosaic.nc --src_restart src_restart/land.res          ",
  "          --dst_restart dst_restart/land.res --dst_cold_restart dst_cold_restart/land.res ",
  "          --remap_file remap_file_C192_to_C384 --print_memory                     ",
  " 3. Remap soil, snow, cana, glac, lake, vegn1, vegn2                              ",
  "      foreach type ( soil snow cana glac lake vegn1 vegn2 )                       ",
  "         if( $type == 'vegn1' || $type == 'vegn2' ) then                          ",
  "             set filetype = 'vegn'                                                ",
  "         else                                                                     ",
  "             set filetype = $type                                                 ",
  "         endif                                                                    ",
  "         >remap_land --file_type $filetype --src_mosaic C192_mosaic.nc            ",
  "             --dst_mosaic c384_mosaic.nc --src_restart src_restart/$type.res      ",
  "             --dst_restart dst_restart/$type.res --dst_cold_restart dst_cold_restart/$type.res    ",
  "             --land_src_restart src_restart/land.res --land_cold_restart dst_cold_restart/land.res ",
  "             --remap_file remap_file_C192_to_C384 --print_memory                                   ",
  "      end                                                                                          ",
  " 4. Remap static_vegn if needed                                                                    ",
  "    >remap_land --file_type vegn --src_mosaic C192_mosaic.nc                                       ",
  "       --dst_mosaic c384_mosaic.nc --src_restart src_restart/static_vegn                           ",
  "       --dst_restart dst_restart/static_vegn --dst_cold_restart dst_cold_restart/vegn1.res         ",
  "       --land_src_restart src_restart/land.res --land_cold_restart dst_cold_restart/land.res       ",
  "       --remap_file remap_file_C192_to_C384 --print_memory                                         ",
  "                                                                                                   ",
  NULL };

char grid_version[] = "0.2";
char tagname[] = "$Name: fre-nctools-bronx-10 $";

double distance(double lon1, double lat1, double lon2, double lat2);
void get_actual_file_name(int nface, int face, const char *file_orig, char *file);
void get_land_tile_info(int fid, const char *name1, const char *name2, int nidx, const int *idx_in, const double *frac_in,
                        int nx, int ny, int ntile, int isc, int iec, int *count, double *frac, int *tag1, int *tag2, int *idx, int all_tile);
void full_search_nearest(int nface_src, int npts_src, const double *lon_src, const double *lat_src,
                         const int *mask_src, int npts_dst, const double *lon_dst, const double *lat_dst,
                         const int *mask_dst, int *idx_map, int *face_map);
void compress_int_data(int ntile, int npts, int nidx, int nidx_global, const int *land_count,
                       const int *data, int *data_global, int all_tile );
void compress_double_data(int ntile, int npts, int nidx, int nidx_global, const int *land_count,
                          const double *data, double *data_global, int all_tile );
void read_cohort_index(int fid, const char * field_name, int* coho_idx);
int compute_dst_coho_idx(int ns_dst, int nt_dst, int nc_dst, int ns_src, int nt_src, int nc_src, int nface_dst,
                         int *idx_map_soil, int *face_map_soil,  int* soil_count_cold, int* ncoho_idx_src, int *coho_idx_src,
                         int *coho_idx_dst, int * coho_idx_d_to_s,  int *coho_idx_face_d_to_s);



const int LANDTYPE = 1;
const int SOILTYPE = 2;
const int GLACTYPE = 3;
const int LAKETYPE = 4;
const int VEGNTYPE = 5;
const int CANATYPE = 6;
const int SNOWTYPE = 7;
const int FIRETYPE = 8;

int main(int argc, char *argv[]) {
  char *src_mosaic = NULL;
  char *dst_mosaic = NULL;
  char *src_restart_file = NULL;
  char *dst_restart_file = NULL;
  char *dst_cold_restart = NULL;
  char *land_src_restart = NULL;
  char *land_cold_restart = NULL;
  char *remap_file = NULL;
  char *file_type = NULL;

  int face_dst;
  int nface_src = 0, nface_dst = 0;
  int nx_src, ny_src, nx_dst, ny_dst;
  int ntile_src, ntile_cold, ntile_dst;
  int nidx_tot_src, ncohort;
  int *fid_src = NULL;
  int *cohort_data = NULL, *tile_axis_data = NULL;
  int *nidx_src = NULL, *nidx_land_src = NULL;

  double *sc_cohort_data = NULL, *lc_cohort_data = NULL;
  int n_sc_cohort, n_lc_cohort = 0;

  double *mdf_day_data = NULL; //Forfire data
  int n_mdf_day = 0;

  double *nspecies_data = NULL;
  int n_nspecies = 0;
  double *textlen_data = NULL;
  int n_textlen = 0;
  int *spc_idx_src = NULL;  // TODO: put such indecies in structure.

  int *idx_soil_src = NULL, *idx_glac_src = NULL, *idx_lake_src = NULL;
  int *soil_count_src = NULL, *glac_count_src = NULL, *lake_count_src = NULL;
  int *soil_tag_src = NULL, *glac_tag_src = NULL, *lake_tag_src = NULL, *vegn_tag_src = NULL;
  double *soil_frac_src = NULL, *glac_frac_src = NULL, *lake_frac_src = NULL;

  int *ncoho_idx_src = NULL;
  int *coho_idx_src = NULL;
  int ncoho_idx_tot_src = 0;

  int filetype;
  char history[1280];
  double *x_src = NULL, *y_src = NULL;
  int time_exist, zaxis_exist, ntime;
  int *has_taxis = NULL, *var_type = NULL, *ndim_src = NULL, *nz_src = NULL;
  int *zld_pos_src = NULL;
  int *tiled_pos_src = NULL;
  int *tileidx_pos_src = NULL;
  double *time_data = NULL;
  int* has_coho_idx = NULL;
  int has_glac = 0, has_lake = 0;
  int src_has_tile = 0, cold_has_tile = 0;
  int src_has_cohort = 0, cold_has_cohort = 0;
  //int src_has_sc_cohort = 0, src_has_lc_cohort = 0;
  //int src_has_nspecies = 0;
  //int src_has_textlen = 0;
  int nz;
  double *z_axis_data = NULL;
  int print_memory = 0;
  int use_all_tile;
  char timename[32] = "NONE";
  int npes, nvar_src;
  int option_index, c;
  int errflg = (argc == 1);
  /*
   * process command line
   */

  static struct option long_options[] = {
                                         {"src_mosaic", required_argument, NULL, 'a'},
                                         {"dst_mosaic", required_argument, NULL, 'b'},
                                         {"src_restart", required_argument, NULL, 'c'},
                                         {"dst_restart", required_argument, NULL, 'o'},
                                         {"dst_cold_restart", required_argument, NULL, 'd'},
                                         {"file_type", required_argument, NULL, 't'},
                                         {"land_src_restart", required_argument, NULL, 'l'},
                                         {"land_cold_restart", required_argument, NULL, 'm'},
                                         {"remap_file", required_argument, NULL, 'r'},
                                         {"print_memory", no_argument, NULL, 'p'},
                                         {"help", no_argument, NULL, 'h'},
                                         {0, 0, 0, 0},
  };

  /* start parallel */

  mpp_init(&argc, &argv);

  mpp_domain_init();

  while ((c = getopt_long(argc, argv, "", long_options, &option_index)) != -1) switch (c) {
    case 'a':
      src_mosaic = optarg;
      break;
    case 'b':
      dst_mosaic = optarg;
      break;
    case 'c':
      src_restart_file = optarg;
      break;
    case 'o':
      dst_restart_file = optarg;
      break;
    case 'd':
      dst_cold_restart = optarg;
      break;
    case 'l':
      land_src_restart = optarg;
      break;
    case 'm':
      land_cold_restart = optarg;
      break;
    case 't':
      file_type = optarg;
      break;
    case 'r':
      remap_file = optarg;
      break;
    case 'p':
      print_memory = 1;
      break;
    case '?':
      errflg++;
      break;
    }

  if (!src_mosaic) errflg++;
  if (!dst_mosaic) errflg++;
  if (!src_restart_file) errflg++;
  if (!dst_restart_file) errflg++;
  if (!dst_cold_restart) errflg++;
  if (!file_type) errflg++;

  if (errflg) {
    char **u = usage;
    if (mpp_pe() == mpp_root_pe()) {
      while (*u) {
        fprintf(stderr, "%s\n", *u);
        u++;
      }
      if (!src_mosaic) mpp_error("remap_land: src_mosaic is not specified");
      if (!dst_mosaic) mpp_error("remap_land: dst_mosaic is not specified");
      if (!src_restart_file) mpp_error("remap_land: src_restart_file is not specified");
      if (!dst_restart_file) mpp_error("remap_land: dst_restart_file is not specified");
      if (!dst_cold_restart) mpp_error("remap_land: dst_cold_restart is not specified");
      if (!file_type) mpp_error("remap_land: file_type is not specified");
    }
    mpp_error("remap_land: check the command line arguments");
  }

  if (print_memory) print_mem_usage("at the begining of remap_land");

  /* write out arguments */
  if (mpp_pe() == mpp_root_pe()) {
    printf("src_mosaic       is %s\n", src_mosaic);
    printf("dst_mosaic       is %s\n", dst_mosaic);
    printf("src_restart_file is %s\n", src_restart_file);
    printf("dst_restart_file is %s\n", dst_restart_file);
    printf("dst_cold_restart is %s\n", dst_cold_restart);
    printf("file_type        is %s\n", file_type);
    if (land_src_restart) {
      printf("land_src_restart is %s\n", land_src_restart);
      printf("land_cold_restart is %s\n", land_cold_restart);
    } else {
      printf("land_src_restart is not specified\n");
      printf("land_cold_restart is not specified\n");
    }
  }

  /* file type must be the land, cana, snow, soil, lake, glac, vegn, fire */
  if (!strcmp(file_type, "land"))
    filetype = LANDTYPE;
  else if (!strcmp(file_type, "soil"))
    filetype = SOILTYPE;
  else if (!strcmp(file_type, "cana"))
    filetype = CANATYPE;
  else if (!strcmp(file_type, "snow"))
    filetype = SNOWTYPE;
  else if (!strcmp(file_type, "lake"))
    filetype = LAKETYPE;
  else if (!strcmp(file_type, "glac"))
    filetype = GLACTYPE;
  else if (!strcmp(file_type, "vegn"))
    filetype = VEGNTYPE;
  else if (!strcmp(file_type, "fire"))
    filetype = FIRETYPE;
  else
    mpp_error("remap_land: invalid option in --file_type");

  /* when file_type is cana or snow, land_src_restart and land_cold_restart must
   * be specified */
  if (filetype == LANDTYPE) {
    if (land_src_restart)
      mpp_error( "remap_land: land_src_restart must not be specified when file_type is 'land'");
    if (land_cold_restart)
      mpp_error( "remap_land: land_cold_restart must not be specified  when file_type is 'land'");
    land_src_restart = src_restart_file;
    land_cold_restart = dst_cold_restart;
  } else {
    if (!land_src_restart)
      mpp_error("remap_land: land_src_restart must be specified when file_type is not 'land'");
    if (!land_cold_restart)
      mpp_error( "remap_land: land_cold_restart must be specified when file_type is not 'land'");
  }

  npes = mpp_npes();

  /*----------------------------------------------------------------------------
    get source and destination grid size
    --------------------------------------------------------------------------*/
  {
    int *nx, *ny;
    nface_src = read_mosaic_ntiles(src_mosaic);
    nface_dst = read_mosaic_ntiles(dst_mosaic);
    nx = (int *)malloc(nface_src * sizeof(int));
    ny = (int *)malloc(nface_src * sizeof(int));
    read_mosaic_grid_sizes(src_mosaic, nx, ny);
    /* nx, ny of source should have the same value on each face */
    for (int n = 1; n < nface_src; n++) {
      if (nx[n] != nx[0] || ny[n] != ny[0])
        mpp_error("remap_land: all the faces of source grid should have the same number of grid points");
    }
    nx_src = nx[0];
    ny_src = ny[0];
    free(nx);
    free(ny);
    nx = (int *)malloc(nface_dst * sizeof(int));
    ny = (int *)malloc(nface_dst * sizeof(int));
    read_mosaic_grid_sizes(dst_mosaic, nx, ny);
    /* nx, ny of source should have the same value on each face */
    for (int n = 1; n < nface_dst; n++) {
      if (nx[n] != nx[0] || ny[n] != ny[0])
        mpp_error("remap_land: all the faces of destination grid should have the same number of grid points");
    }
    nx_dst = nx[0];
    ny_dst = ny[0];
    free(nx);
    free(ny);
  }

  /*-----------------------------------------------------------------------------
    read the source tile_index size and source grid
    ----------------------------------------------------------------------------*/
  {

    fid_src = (int *)malloc(nface_src * sizeof(int));
    nidx_src = (int *)malloc(nface_src * sizeof(int));
    nidx_land_src = (int *)malloc(nface_src * sizeof(int));
    nidx_tot_src = 0;
    for (int n = 0; n < nface_src; n++) {
      int nlon, nlat;
      char file[512];

      get_actual_file_name(nface_src, n, src_restart_file, file);
      fid_src[n] = mpp_open(file, MPP_READ);
      nlon = mpp_get_dimlen(fid_src[n], LON_NAME);
      nlat = mpp_get_dimlen(fid_src[n], LAT_NAME);
      if (nx_src != nlon)
        mpp_error(
                  "remap_land: mismatch on the longitude dimension size "
                  "between source mosaic grid file and src_restart");
      if (ny_src != nlat)
        mpp_error(
                  "remap_land: mismatch on the latitude dimension size "
                  "between source mosaic grid file and src_restart");
      nidx_src[n] = mpp_get_dimlen(fid_src[n], TILE_INDEX_NAME);

      if (filetype == LANDTYPE)
        nidx_land_src[n] = nidx_src[n];
      else {
        get_actual_file_name(nface_src, n, land_src_restart, file);
        int fid_land = mpp_open(file, MPP_READ);

        nidx_land_src[n] = mpp_get_dimlen(fid_land, TILE_INDEX_NAME);
        mpp_close(fid_land);
      }

      nidx_tot_src += nidx_src[n];
    }

    x_src = (double *)malloc(nface_src * nx_src * ny_src * sizeof(double));
    y_src = (double *)malloc(nface_src * nx_src * ny_src * sizeof(double));

    for (int n = 0; n < nface_src; n++) {
      read_mosaic_grid_data(src_mosaic, "x", nx_src, ny_src, x_src + n * nx_src * ny_src, n, 1, 1);
      read_mosaic_grid_data(src_mosaic, "y", nx_src, ny_src, y_src + n * nx_src * ny_src, n, 1, 1);
    }
    /* convert to radians */
    for (int n = 0; n < nface_src * nx_src * ny_src; n++) {
      x_src[n] *= D2R;
      y_src[n] *= D2R;
    }
  }

  /*-----------------------------------------------------------------------------
    Get the time information, time axis only exists for static vegetation.
    ---------------------------------------------------------------------------*/
  time_exist = mpp_get_record_name(fid_src[0], timename);
  ntime = 1;
  if (time_exist) {
    int vid;

    ntime = mpp_get_dimlen(fid_src[0], timename);
    time_data = (double *)malloc(ntime * sizeof(double));
    vid = mpp_get_varid(fid_src[0], timename);
    mpp_get_var_value(fid_src[0], vid, time_data);
  }
  /* only static vegetation can have more than 1 time level */
  if (ntime != 1 && filetype != VEGNTYPE)
    mpp_error("remap_land: ntime should be 1 when file_type is not vegn");

  /* check if z-axis exist or not */
  zaxis_exist = mpp_dim_exist(fid_src[0], LEVEL_NAME);
  nz = 1;
  if (zaxis_exist) {
    int vid;

    nz = mpp_get_dimlen(fid_src[0], LEVEL_NAME);
    z_axis_data = (double *)malloc(nz * sizeof(double));
    vid = mpp_get_varid(fid_src[0], LEVEL_NAME);
    mpp_get_var_value(fid_src[0], vid, z_axis_data);
  }

  /*-----------------------------------------------------------------------------
    loop through each variable of the source data to see get dimension of each variable
    ----------------------------------------------------------------------------*/
  // TODO: Next version : put such indecies below in structure or class.
  nvar_src = mpp_get_nvars(fid_src[0]);
  has_taxis = (int *)malloc(nvar_src * sizeof(int));
  ndim_src = (int *)malloc(nvar_src * sizeof(int));
  nz_src = (int *)malloc(nvar_src * sizeof(int));
  zld_pos_src = (int *)malloc(nvar_src * sizeof(int));  //zlevel dimension position (-1 if none) //ML
  tiled_pos_src = (int *)malloc(nvar_src * sizeof(int)); //tile dimesion position (-1 if no tile)
  tileidx_pos_src = (int *)malloc(nvar_src * sizeof(int)); //tile index position (-i if none)
  var_type = (int *)malloc(nvar_src * sizeof(int));
  spc_idx_src = (int *)malloc(nvar_src * sizeof(int));
  has_coho_idx = (int *)malloc(nvar_src * sizeof(int));


  printf("**RL starting determination of dim of each var, nvar_scr=%d\n", nvar_src);

  for (int l = 0; l < nvar_src; l++) {
    char varname[64];
    char vdname[64];  // var dimname
    int vid;
    //int klev;

    mpp_get_varname(fid_src[0], l, varname);
    vid = mpp_get_varid(fid_src[0], varname);

    var_type[l] = mpp_get_var_type(fid_src[0], vid);
    ndim_src[l] = mpp_get_var_ndim(fid_src[0], vid);

    has_taxis[l] = 0;
    nz_src[l] = -1;
    zld_pos_src[l] = -1;
    tiled_pos_src[l] = -1;
    tileidx_pos_src[l] = -1;
    has_coho_idx[l] = 0;

    if (ndim_src[l] > 3) {
      mpp_error("remap_land: more than 3 dimensions for the field in src_restart");
    }

    if (var_type[l] == MPP_INT || var_type[l] == MPP_DOUBLE) {

      for (int m = 0; m < ndim_src[l]; m++) {
        mpp_get_var_dimname(fid_src[0], vid, m, vdname);
        if (!strcmp(vdname, timename)) {
          has_taxis[l] = 1;
        } else if (!strcmp(vdname, LEVEL_NAME)) {
          zld_pos_src[l] = m;
          nz_src[l] = mpp_get_dimlen(fid_src[0], vdname);
        } else if (!strcmp(vdname, TILE_INDEX_NAME)) {
          tileidx_pos_src[l] = m;
        } else if (!strcmp(vdname, TILE_NAME)) {
          tiled_pos_src[l] = m;
        } else if (!strcmp(vdname, COHORT_INDEX_NAME)) {
          has_coho_idx[l] = 1;
        }
        printf("vname=%s  dname=%s vtype=%d ndim=%d\n", varname, vdname, var_type[l], ndim_src[l]);
      }

      //TODO: Finalize removal of these checks originally from lm4.0
      //  (nz_src[l] == -1)  (ndim_src[l] == 3 && !has_taxis[l]) {
      // TODO: What is the consequence of not having a time dimension.
    } else if (var_type[l] == MPP_CHAR) {
      for (int m = 0; m < ndim_src[l]; m++) {
        mpp_get_var_dimname(fid_src[0], vid, m, vdname);
        if (!strcmp(vdname, timename)) {
          has_taxis[l] = 1;
        } else if (!strcmp(vdname, LEVEL_NAME)) {
          zld_pos_src[l] = m;
          nz_src[l] = mpp_get_dimlen(fid_src[0], vdname);
        }else if (!strcmp(vdname, TILE_INDEX_NAME)) {
          tileidx_pos_src[l] = m;
        } else if (!strcmp(vdname, TILE_NAME)) {
          tiled_pos_src[l] = m;
        } else if (!strcmp(vdname, NSPECIES_NAME)) {
          spc_idx_src[l] = m;
        }
        printf("vname=%s  dname=%s vtype=%d ndim=%d\n", varname, vdname, var_type[l], ndim_src[l]);
      }
    } else {
      mpp_error("remap_land: field type must be MPP_INT or MPP_DOUBLE or MPP_CHAR");
    }
    printf("vname=%s l=%d vid=%d vtype=%d ndim=%d taxis=%d zld =%d nz=%d tiled_pos=%d tidx_pos=%d\n",
           varname, l, vid, var_type[l], ndim_src[l], has_taxis[l], zld_pos_src[l], nz_src[l],
           tiled_pos_src[l], tileidx_pos_src[l]);
  }

  /*------------------------------------------------------------------------------
    get the cohort data, currently it only contains in vegn data
    -----------------------------------------------------------------------------*/
  {
    ncohort = 0;
    src_has_cohort = 0;
    if (filetype == VEGNTYPE) {
      if (!mpp_dim_exist(fid_src[0], COHORT_NAME)) {
        mpp_error("remap_land: dimension cohort should exist when file_type is vegn");
      }

      ncohort = mpp_get_dimlen(fid_src[0], COHORT_NAME);

      if (mpp_var_exist(fid_src[0], COHORT_NAME)) {
        src_has_cohort = 1;
        cohort_data = (int *)malloc(ncohort * sizeof(int));
        int vid = mpp_get_varid(fid_src[0], COHORT_NAME);
        mpp_get_var_value(fid_src[0], vid, cohort_data);
        for (int m = 1; m < nface_src; m++) {
          int dimsize = mpp_get_dimlen(fid_src[m], COHORT_NAME);
          if (dimsize != ncohort) mpp_error("remap_land: the dimension size of cohort is different between faces");
          int * tmp = (int *)malloc(ncohort * sizeof(int));
          vid = mpp_get_varid(fid_src[m], COHORT_NAME);
          mpp_get_var_value(fid_src[m], vid, tmp);
          for (int i = 0; i < ncohort; i++) {
            if (cohort_data[i] != tmp[i]) {
              mpp_error("remap_land: cohort value is different between faces");
            }
          }
          free(tmp);
        }
      }
    }
  }

  /*------------------------------------------------------------------------------
    get the sc_cohort data, currently it only contains in vegn data
    -----------------------------------------------------------------------------*/
  {
    n_sc_cohort = 0;
    if ((filetype == VEGNTYPE) || (filetype == SOILTYPE)) {
      if (mpp_var_exist(fid_src[0], SC_COHORT_NAME)) {
        //src_has_sc_cohort = 1;
        n_sc_cohort = mpp_get_dimlen(fid_src[0], SC_COHORT_NAME);
        sc_cohort_data = (double *)malloc(n_sc_cohort * sizeof(double));
        int vid = mpp_get_varid(fid_src[0], SC_COHORT_NAME);
        mpp_get_var_value(fid_src[0], vid, sc_cohort_data);
        for (int m = 1; m < nface_src; m++) {
          int dimsize = mpp_get_dimlen(fid_src[m], SC_COHORT_NAME);
          if (dimsize != n_sc_cohort)
            mpp_error("remap_land: the dimension size of sc_cohort is different between faces");
          double *tmp = (double *)malloc(n_sc_cohort * sizeof(double));
          vid = mpp_get_varid(fid_src[m], SC_COHORT_NAME);
          mpp_get_var_value(fid_src[m], vid, tmp);
          for (int i = 0; i < n_sc_cohort; i++) {
            if (sc_cohort_data[i] != tmp[i]) {
              mpp_error("remap_land: sc_cohort value is different between faces");
            }
          }
          free(tmp);
        }
      }
    }
  }

  /*------------------------------------------------------------------------------
    get the lc_cohort data, currently it only contains in vegn data
    -----------------------------------------------------------------------------*/

  {
    n_lc_cohort = 0;
    //src_has_lc_cohort = 0;
    if ((filetype == VEGNTYPE) || (filetype == SOILTYPE)) {
      if (mpp_var_exist(fid_src[0], LC_COHORT_NAME)) {
        n_lc_cohort = mpp_get_dimlen(fid_src[0], LC_COHORT_NAME);
        printf("*RLE n_lc_cohort= %d\n", n_lc_cohort);
        lc_cohort_data = (double *)malloc(n_lc_cohort * sizeof(double));
        int vid = mpp_get_varid(fid_src[0], LC_COHORT_NAME);
        mpp_get_var_value(fid_src[0], vid, lc_cohort_data);
        for (int m = 1; m < nface_src; m++) {
          int dimsize = mpp_get_dimlen(fid_src[m], LC_COHORT_NAME);
          if (dimsize != n_lc_cohort)
            mpp_error("remap_land: the dimension size of lc_cohort is different between faces");
          double *tmp = (double *)malloc(n_lc_cohort * sizeof(double));
          vid = mpp_get_varid(fid_src[m], LC_COHORT_NAME);
          mpp_get_var_value(fid_src[m], vid, tmp);
          for (int i = 0; i < n_lc_cohort; i++) {
            if (lc_cohort_data[i] != tmp[i]) {
              mpp_error("remap_land: lc_cohort value is different between faces");
            }
          }
          free(tmp);
        }
      }
    }
  }

  /*------------------------------------------------------------------------------
    get the nspecies data
    ------------------------------------------------------------------------*/
  {
    n_nspecies = 0;
    if ((filetype == VEGNTYPE) || (filetype == SOILTYPE)) {

      if (mpp_var_exist(fid_src[0], NSPECIES_NAME)) {
        n_nspecies = mpp_get_dimlen(fid_src[0], NSPECIES_NAME);
        nspecies_data = (double *)malloc(n_nspecies * sizeof(double));
        int vid = mpp_get_varid(fid_src[0], NSPECIES_NAME);
        mpp_get_var_value(fid_src[0], vid, nspecies_data);
        for (int m = 1; m < nface_src; m++) {
          int dimsize = mpp_get_dimlen(fid_src[m], NSPECIES_NAME);
          if (dimsize != n_nspecies){
            mpp_error("remap_land: the dimension size of nspecies is different between faces");
          }
          double *tmp = (double *)malloc(n_nspecies * sizeof(double));
          vid = mpp_get_varid(fid_src[m], NSPECIES_NAME);
          mpp_get_var_value(fid_src[m], vid, tmp);
          for (int i = 0; i < n_nspecies; i++) {
            if (nspecies_data[i] != tmp[i]) {
              free(tmp);
              tmp = NULL;
              mpp_error("remap_land: nspecies value is different between faces");
            }
          }
          free(tmp);
          tmp = NULL;
        }
      }
    }
  }

  /*-----------------------------------------------------------------------------
    get the textlen data
    -------------------------------------------------------------------------------*/
  {
    n_textlen = 0;
    if ((filetype == VEGNTYPE) || (filetype == SOILTYPE)) {
      if (mpp_var_exist(fid_src[0], TEXTLEN_NAME)) {
        n_textlen = mpp_get_dimlen(fid_src[0], TEXTLEN_NAME);
        textlen_data = (double *)malloc(n_textlen * sizeof(double));
        int vid = mpp_get_varid(fid_src[0], TEXTLEN_NAME);
        mpp_get_var_value(fid_src[0], vid, textlen_data);
        for (int m = 1; m < nface_src; m++) {
          int dimsize = mpp_get_dimlen(fid_src[m], TEXTLEN_NAME);
          if (dimsize != n_textlen) mpp_error("remap_land: the dimension size of textlen is different between faces");
          double *tmp = (double *)malloc(n_textlen * sizeof(double));
          vid = mpp_get_varid(fid_src[m], TEXTLEN_NAME);
          mpp_get_var_value(fid_src[m], vid, tmp);
          for (int i = 0; i < n_textlen; i++) {
            if (textlen_data[i] != tmp[i]) {
              free(tmp);
              tmp = NULL;
              mpp_error("remap_land: textlen value is different between faces");
            }
          }
          free(tmp);
          tmp = NULL;
        }
      }
    }
  }

 /*------------------------------------------------------------------------------
    get the mdf_dat data
    -----------------------------------------------------------------------------*/

  {
    n_mdf_day = 0;
    if (filetype == FIRETYPE) {
      if (mpp_var_exist(fid_src[0], MDF_DAY_NAME)) {
        n_mdf_day = mpp_get_dimlen(fid_src[0], MDF_DAY_NAME);
        printf("*RLE n_mdf_day= %d\n", n_mdf_day);
        mdf_day_data = (double *)malloc(n_mdf_day * sizeof(double));
        int vid = mpp_get_varid(fid_src[0], MDF_DAY_NAME);
        mpp_get_var_value(fid_src[0], vid, mdf_day_data);
        for (int m = 1; m < nface_src; m++) {
          int dimsize = mpp_get_dimlen(fid_src[m], MDF_DAY_NAME);
          if (dimsize != n_mdf_day)
            mpp_error("remap_land: the dimension size of mdf_day is different between faces");
          double *tmp = (double *)malloc(n_mdf_day * sizeof(double));
          vid = mpp_get_varid(fid_src[m], MDF_DAY_NAME);
          mpp_get_var_value(fid_src[m], vid, tmp);
          for (int i = 0; i < n_mdf_day; i++) {
            if (mdf_day_data[i] != tmp[i]) {
              mpp_error("remap_land: mdf_day value is different between faces");
            }
          }
          free(tmp);
        }
      }
    }
  }

/*-------------------------------
  read the cohort index. Unlike cohort field, the index does not have to have the same elements in all faces.
 -------------------------------------------*/
  //mem allocated *
   printf("Calling read_cohort_index for filetype=%d\n", filetype);

  if (filetype == VEGNTYPE) {
    ncoho_idx_src = (int *)malloc(nface_src * sizeof(int));

    //Get the source total number of cohort indices in all (usually six) faces:
    ncoho_idx_tot_src = 0;
    for (int n = 0; n < nface_src; n++) {
      if (!mpp_dim_exist(fid_src[n], COHORT_INDEX_NAME)) {
        mpp_error("remap_land: dimension cohort_index should exist when file_type is vegn");
      }
      ncoho_idx_src[n] = mpp_get_dimlen(fid_src[n], COHORT_INDEX_NAME);
      ncoho_idx_tot_src += ncoho_idx_src[n];
    }

    //Get the values of all the src cohort indices in all (usually six) faces:
    coho_idx_src = (int *)malloc(ncoho_idx_tot_src * sizeof(int));
    int pos_coho = 0;
    for (int n = 0; n < nface_src; n++) {
      int vid = mpp_get_varid(fid_src[n], COHORT_INDEX_NAME);
      mpp_get_var_value(fid_src[n], vid, coho_idx_src + pos_coho);
      pos_coho += ncoho_idx_src[n];
    }
  }


  /*-----------------------------------------------------------------------------
    Check if the cold_restart has lake or glac. soil is required to be existed
    ---------------------------------------------------------------------------*/
  {
    int max_nidx, vid;
    int *fid = NULL;
    char file[512];
    int *nidx = NULL;
    int *tmp = NULL;

    cold_has_tile = 0;
    cold_has_cohort = 0;
    has_lake = 0;
    has_glac = 0;
    max_nidx = 0;
    fid = (int *)malloc(nface_dst * sizeof(int));
    nidx = (int *)malloc(nface_dst * sizeof(int));

    /* get ntile for dst_cold_restart */
    {
      int fid_cold;

      get_actual_file_name(nface_dst, 0, dst_cold_restart, file);
      fid_cold = mpp_open(file, MPP_READ);
      cold_has_tile = mpp_var_exist(fid_cold, TILE_NAME);
      if (filetype == VEGNTYPE) {
        cold_has_cohort = mpp_var_exist(fid_cold, COHORT_NAME);
      }
      mpp_close(fid_cold);
    }

    for (face_dst = 0; face_dst < nface_dst; face_dst++) {
      get_actual_file_name(nface_dst, face_dst, land_cold_restart, file);
      fid[face_dst] = mpp_open(file, MPP_READ);
      nidx[face_dst] = mpp_get_dimlen(fid[face_dst], TILE_INDEX_NAME);
      if (nidx[face_dst] > max_nidx) max_nidx = nidx[face_dst];
    }

    ntile_cold = mpp_get_dimlen(fid[0], TILE_NAME);
    for (face_dst = 1; face_dst < nface_dst; face_dst++) {
      int ntile;
      ntile = mpp_get_dimlen(fid[face_dst], TILE_NAME);
      if (ntile != ntile_cold) mpp_error("remap_land: mismatch of tile dimension between different faces");
    }

    tmp = (int *)malloc(max_nidx * sizeof(int));
    for (face_dst = 0; face_dst < nface_dst; face_dst++) {
      vid = mpp_get_varid(fid[face_dst], GLAC_NAME);
      mpp_get_var_value(fid[face_dst], vid, tmp);
      for (int i = 0; i < nidx[face_dst]; i++) {
        if (tmp[i] > 0) {
          has_glac = 1;
          goto GLAC_CHECK;
        }
      }
    }
  GLAC_CHECK:
    if (mpp_pe() == mpp_root_pe()) {
      if (has_glac)
        printf("remap_land: there is glac in cold restart file\n");
      else
        printf("remap_land: there is no glac in cold restart file\n");
    }

    for (face_dst = 0; face_dst < nface_dst; face_dst++) {
      vid = mpp_get_varid(fid[face_dst], LAKE_NAME);
      mpp_get_var_value(fid[face_dst], vid, tmp);
      for (int i = 0; i < nidx[face_dst]; i++) {
        if (tmp[i] > 0) {
          has_lake = 1;
          goto LAKE_CHECK;
        }
      }
    }
  LAKE_CHECK:
    if (mpp_pe() == mpp_root_pe()) {
      if (has_lake)
        printf("remap_land: there is lake in cold restart file\n");
      else
        printf("remap_land: there is no lake in cold restart file\n");
    }
    for (face_dst = 0; face_dst < nface_dst; face_dst++) mpp_close(fid[face_dst]);

    free(tmp);
    free(nidx);
    free(fid);
  }

  /*-----------------------------------------------------------------------------
    get the tile data
    ----------------------------------------------------------------------------*/
  printf("remap_land : getting tile data\n");

  if (!mpp_dim_exist(fid_src[0], TILE_NAME)) mpp_error("remap_land: dimension tile should exist");

  ntile_src = mpp_get_dimlen(fid_src[0], TILE_NAME);
  ntile_dst = ntile_src;
  /*  if(!has_glac) ntile_dst--; */
  /* if(!has_lake) ntile_dst--;  */

  src_has_tile = 0;
  if (mpp_var_exist(fid_src[0], TILE_NAME)) {
    src_has_tile = 1;
    tile_axis_data = (int *)malloc(ntile_src * sizeof(int));
    int vid = mpp_get_varid(fid_src[0], TILE_NAME);
    mpp_get_var_value(fid_src[0], vid, tile_axis_data);
    for (int m = 1; m < nface_src; m++) {
      int dimsize = mpp_get_dimlen(fid_src[m], TILE_NAME);
      if (dimsize != ntile_src) mpp_error("remap_land: the dimension size of tile is different between faces");
      int* tmp = (int *)malloc(ntile_src * sizeof(int));
      vid = mpp_get_varid(fid_src[m], TILE_NAME);
      mpp_get_var_value(fid_src[m], vid, tmp);
      for (int i = 0; i < ntile_src; i++) {
        if (tile_axis_data[i] != tmp[i]) mpp_error("remap_land: tile value is different between faces");
      }
      free(tmp);
    }
  }

  /*-----------------------------------------------------------------------
    get source tile information.
    get the tile information for soil, lake and glac
    The following will use large memroy for high resolution source restart file,
    May come back to rewrite to decrease memory usage.
    ---------------------------------------------------------------------*/
  {
    soil_count_src = (int *)malloc(nface_src * nx_src * ny_src * sizeof(int));
    soil_frac_src = (double *)malloc(nface_src * nx_src * ny_src * ntile_src * sizeof(double));
    soil_tag_src = (int *)malloc(nface_src * nx_src * ny_src * ntile_src * sizeof(int));
    vegn_tag_src = (int *)malloc(nface_src * nx_src * ny_src * ntile_src * sizeof(int));
    idx_soil_src = (int *)malloc(nface_src * nx_src * ny_src * ntile_src * sizeof(int));

    if (has_glac) {
      glac_count_src = (int *)malloc(nface_src * nx_src * ny_src * sizeof(int));//ML and 4 below
      glac_frac_src = (double *)malloc(nface_src * nx_src * ny_src * sizeof(double)); /* at most one tile for glac */
      glac_tag_src = (int *)malloc(nface_src * nx_src * ny_src * sizeof(int));        /* at most one tile for glac */
      idx_glac_src = (int *)malloc(nface_src * nx_src * ny_src * ntile_src * sizeof(int));
    }
    if (has_lake) {
      lake_count_src = (int *)malloc(nface_src * nx_src * ny_src * sizeof(int));//ML and 4 below
      lake_frac_src = (double *)malloc(nface_src * nx_src * ny_src * sizeof(double)); /* at most one tile for glac */
      lake_tag_src = (int *)malloc(nface_src * nx_src * ny_src * sizeof(int));        /* at most one tile for glac */
      idx_lake_src = (int *)malloc(nface_src * nx_src * ny_src * ntile_src * sizeof(int));
    }



    {
      double *frac_land_src = NULL;
      int *idx_src = NULL;
      int *idx_land_src = NULL;
      char file[512];

      int max_nidx, pos1, pos2;
      max_nidx = 0;

      for (int n = 0; n < nface_src; n++) {
        if (nidx_land_src[n] > max_nidx) max_nidx = nidx_land_src[n];
      }
      frac_land_src = (double *)malloc(max_nidx * sizeof(double));
      idx_land_src = (int *)malloc(max_nidx * sizeof(int));

      pos1 = 0;
      pos2 = 0;
      use_all_tile = 0;
      if (filetype == LANDTYPE || filetype == CANATYPE || filetype == SNOWTYPE) use_all_tile = 1;

      for (int n = 0; n < nface_src; n++) {
        int vid, fid;

        if (filetype == LANDTYPE) {
          fid = fid_src[n];
        } else {
          get_actual_file_name(nface_src, n, land_src_restart, file);
          fid = mpp_open(file, MPP_READ);
        }
        vid = mpp_get_varid(fid, FRAC_NAME);
        mpp_get_var_value(fid, vid, frac_land_src);
        vid = mpp_get_varid(fid, TILE_INDEX_NAME);
        mpp_get_var_value(fid, vid, idx_land_src);

        /* soil, vegn */
        get_land_tile_info(fid, SOIL_NAME, VEGN_NAME, nidx_land_src[n], idx_land_src, frac_land_src, nx_src, ny_src,
                           ntile_src, 0, nx_src * ny_src - 1, soil_count_src + pos1, soil_frac_src + pos2,
                           soil_tag_src + pos2, vegn_tag_src + pos2, idx_soil_src + pos2, use_all_tile);
        /* glac */
        if (has_glac)
          get_land_tile_info(fid, GLAC_NAME, NULL, nidx_land_src[n], idx_land_src, frac_land_src, nx_src, ny_src, 1, 0,
                             nx_src * ny_src - 1, glac_count_src + pos1, glac_frac_src + pos1, glac_tag_src + pos1,
                             NULL, idx_glac_src + pos1, use_all_tile);
        /* lake */
        if (has_lake)
          get_land_tile_info(fid, LAKE_NAME, NULL, nidx_land_src[n], idx_land_src, frac_land_src, nx_src, ny_src, 1, 0,
                             nx_src * ny_src - 1, lake_count_src + pos1, lake_frac_src + pos1, lake_tag_src + pos1,
                             NULL, idx_lake_src + pos1, use_all_tile);
        pos1 += nx_src * ny_src;
        pos2 += nx_src * ny_src * ntile_src;

        if (filetype != LANDTYPE) mpp_close(fid);

        /* make sure the tile_index consistency between src_restart_file and
         * land_src_restart */
        if (filetype != LANDTYPE) {
          if (filetype == CANATYPE || filetype == SNOWTYPE) {
            if (nidx_land_src[n] != nidx_src[n])
              mpp_error("remap_land: size of tile_index mismatch between "
                        "       src_restart_file and land_src_restart for 'cana' or 'snow'");
          } else {
            if (nidx_land_src[n] < nidx_src[n])
              mpp_error("remap_land: size of tile_index mismatch between "
                        "src_restart_file and land_src_restart for 'soil', 'vegn', 'glac' or 'lake'");
          }
          idx_src = (int *)malloc(max_nidx * sizeof(int));
          vid = mpp_get_varid(fid_src[n], TILE_INDEX_NAME);
          mpp_get_var_value(fid_src[n], vid, idx_src);
          if (filetype == CANATYPE || filetype == SNOWTYPE) {
            int i;
            for (i = 0; i < nidx_land_src[n]; i++)
              if (idx_src[i] != idx_land_src[i])
                mpp_error("remap_land: mismatch of tile_index between src_restart_file "
                          "and land_src_restart for 'soil', 'vegn', 'glac' or 'lake'");
          } else if (filetype == SOILTYPE || filetype == VEGNTYPE) {
            int i, j, k, l, p, m, idx, p2;
            for (m = 0; m < ntile_src; m++) {
              for (l = 0; l < nx_src * ny_src; l++) {
                p = n * nx_src * ny_src + l;
                idx = idx_soil_src[ntile_src * p + m];
                if (idx == MPP_FILL_INT) continue;
                if (idx > nidx_src[n]) mpp_error("remap_land: idx is out of bound for soil consistency check");
                idx = idx_src[idx];
                i = idx % nx_src;
                k = idx / nx_src;
                j = k % ny_src;
                k = k / ny_src;
                p2 = n * nx_src * ny_src + j * nx_src + i;
                if (p != p2) mpp_error("remap_land: mismatch of tile_index for src soil check");
                if (soil_tag_src[ntile_src * p + k] == MPP_FILL_INT) {
                  // TODO: Investigate
                  // mpp_error("remap_land: soil_tag_src is not defined for src soil check");
                }
              }
            }
          }
        }
      }
      free(frac_land_src);
      free(idx_land_src);
      if (filetype != LANDTYPE) free(idx_src);
    }
  }

  /* define history attribute */
  {
    strcpy(history, argv[0]);
    for (int n = 1; n < argc; n++) {
      strcat(history, " ");
      strcat(history, argv[n]);
    }
  }

  if (print_memory) print_mem_usage("After initialization of source information");

  /*------------------------------------------------------------------------------------------
    loop through each face of destination grid, first read the grid, then read
    the tile_index, then find the remapping index, then setup metadata for the
    destination file, last do the remapping and write out the data to
    dst_restart_file
    ----------------------------------------------------------------------------------------*/
  {
    double *data_dst = NULL, *data_src = NULL;
    int *idata_dst = NULL, *idata_src = NULL;
    int *start_pos = NULL;
    double *x_tmp = NULL, *y_tmp = NULL;
    double *lon_axis_dst = NULL, *lat_axis_dst = NULL;
    double *x_dst = NULL, *y_dst = NULL;

    int *idx_dst = NULL;
    int *glac_tag_dst = NULL, *lake_tag_dst = NULL, *soil_tag_dst = NULL, *vegn_tag_dst = NULL;
    int *idx_map_soil = NULL, *face_map_soil = NULL;
    int *idx_map_glac = NULL, *face_map_glac = NULL;
    int *idx_map_lake = NULL, *face_map_lake = NULL;
    int *land_idx_map = NULL, *land_face_map = NULL;
    double *glac_frac_cold = NULL, *lake_frac_cold = NULL, *soil_frac_cold = NULL, *tmp_frac_cold = NULL;
    int *glac_count_cold = NULL, *lake_count_cold = NULL, *soil_count_cold = NULL;
    int *glac_tag_cold = NULL, *lake_tag_cold = NULL, *soil_tag_cold = NULL, *vegn_tag_cold = NULL;
    int *land_count_dst = NULL;
    double *land_frac_dst = NULL;
    //Cohort index related data uses more space than above (unless gegenerate)
    int *coho_idx_dst = NULL; //The valid destination cohort indecies;
    int *coho_idx_d_to_s = NULL; //Mapping of the dst coho index to the src coho index.
    int *coho_idx_face_d_to_s = NULL; //face of the mapped src coho index.
    int *coho_idata_src = NULL;
    double *coho_rdata_src = NULL;
    int *coho_idata_dst = NULL;
    double *coho_rdata_dst = NULL;

    int isc_dst, iec_dst, jsc_dst, jec_dst, nxc_dst;
    int n, pos;
    int layout[2];
    domain2D Dom_dst;

    lon_axis_dst = (double *)malloc(nx_dst * sizeof(double));
    lat_axis_dst = (double *)malloc(ny_dst * sizeof(double));
    x_tmp = (double *)malloc(nx_dst * ny_dst * sizeof(double));
    y_tmp = (double *)malloc(nx_dst * ny_dst * sizeof(double));
    start_pos = (int *)malloc(nface_src * sizeof(int));

    idata_src = (int *)malloc(nidx_tot_src * sizeof(int));
    data_src = (double *)malloc(nidx_tot_src * sizeof(double));

    /* setup domain */
    layout[0] = npes;
    layout[1] = 1;
    mpp_define_domain2d(nx_dst * ny_dst, 1, layout, 0, 0, &Dom_dst);
    mpp_get_compute_domain2d(Dom_dst, &isc_dst, &iec_dst, &jsc_dst, &jec_dst);
    /* for this domain, jsc_dst should equal jec_dst */
    if (jsc_dst != jec_dst) {
      mpp_error("remap_land: This is a 1-D domain decomposition, jsc_dst must "
                "equal to jec_dst");
    }
    nxc_dst = iec_dst - isc_dst + 1;

    data_dst = (double *)malloc(ntile_dst * nxc_dst * sizeof(double));
    idata_dst = (int *)malloc(ntile_dst * nxc_dst * sizeof(int));
    x_dst = (double *)malloc(nxc_dst * sizeof(double));
    y_dst = (double *)malloc(nxc_dst * sizeof(double));

    land_idx_map = (int *)malloc(ntile_dst * nxc_dst * sizeof(int));
    land_face_map = (int *)malloc(ntile_dst * nxc_dst * sizeof(int));
    land_count_dst = (int *)malloc(nxc_dst * sizeof(int));
    idx_dst = (int *)malloc(ntile_dst * nxc_dst * sizeof(int));

    coho_idata_src = (int *)malloc(ncoho_idx_tot_src * sizeof(int));
    coho_rdata_src = (double *)malloc(ncoho_idx_tot_src * sizeof(int));
    coho_idata_dst = (int *)malloc(ncohort * ntile_dst * nxc_dst * sizeof(int));
    coho_rdata_dst = (double *)malloc(ncohort * ntile_dst * nxc_dst * sizeof(double));
    coho_idx_dst = (int *)malloc(ncohort * ntile_dst * nxc_dst * sizeof(int));
    coho_idx_d_to_s = (int *)malloc(ncohort * ntile_dst * nxc_dst * sizeof(int));
    coho_idx_face_d_to_s  = (int *)malloc(ncohort * ntile_dst * nxc_dst * sizeof(int));

    soil_tag_dst = (int *)malloc(ntile_dst * nxc_dst * sizeof(int));
    vegn_tag_dst = (int *)malloc(ntile_dst * nxc_dst * sizeof(int));
    land_frac_dst = (double *)malloc(ntile_dst * nxc_dst * sizeof(double));
    soil_count_cold = (int *)malloc(nxc_dst * sizeof(int));
    soil_frac_cold = (double *)malloc(nxc_dst * sizeof(double));
    tmp_frac_cold = (double *)malloc(ntile_cold * nxc_dst * sizeof(double));
    soil_tag_cold = (int *)malloc(ntile_cold * nxc_dst * sizeof(int));
    vegn_tag_cold = (int *)malloc(ntile_cold * nxc_dst * sizeof(int));
    idx_map_soil = (int *)malloc(nxc_dst * sizeof(int));
    face_map_soil = (int *)malloc(nxc_dst * sizeof(int));

    if (has_glac) {
      glac_tag_dst = (int *)malloc(ntile_dst * nxc_dst * sizeof(int));
      glac_count_cold = (int *)malloc(nxc_dst * sizeof(int));
      glac_frac_cold = (double *)malloc(nxc_dst * sizeof(double));
      glac_tag_cold = (int *)malloc(nxc_dst * sizeof(int));
      idx_map_glac = (int *)malloc(nxc_dst * sizeof(int));
      face_map_glac = (int *)malloc(nxc_dst * sizeof(int));
    }

    if (has_lake) {
      lake_tag_dst = (int *)malloc(ntile_dst * nxc_dst * sizeof(int));
      lake_count_cold = (int *)malloc(nxc_dst * sizeof(int));
      lake_frac_cold = (double *)malloc(nxc_dst * sizeof(double));
      lake_tag_cold = (int *)malloc(nxc_dst * sizeof(int));
      idx_map_lake = (int *)malloc(nxc_dst * sizeof(int));
      face_map_lake = (int *)malloc(nxc_dst * sizeof(int));
    }

    if (print_memory) print_mem_usage("before the loop of face_dst");
    for (face_dst = 0; face_dst < nface_dst; face_dst++) {
      int nlon, nlat, l, i, k, t; //p,j,ll
      int fid_dst, vid_dst, vid_src, fid_cold;  //vid_cold
      int nidx_dst, nidx_dst_global, ncoho_idx_dst, ncoho_idx_dst_global;
      int nidx_cold, vid;
      int fid_land_cold;
      double *frac_cold = NULL;
      int *idx_cold = NULL;

      double *rdata_global = NULL;
      int *idata_global = NULL;
      int *coho_idata_global = NULL;

      size_t start[4], nread[4], nwrite[4];
      char land_cold[512], file_dst[512], file_cold[512];

      get_actual_file_name(nface_dst, face_dst, dst_restart_file, file_dst);
      get_actual_file_name(nface_dst, face_dst, dst_cold_restart, file_cold);

      fid_cold = mpp_open(file_cold, MPP_READ);
      if (filetype == LANDTYPE) {
        nidx_cold = mpp_get_dimlen(fid_cold, TILE_INDEX_NAME);
        fid_land_cold = fid_cold;
      } else {
        get_actual_file_name(nface_dst, face_dst, land_cold_restart, land_cold);
        fid_land_cold = mpp_open(land_cold, MPP_READ);
        nidx_cold = mpp_get_dimlen(fid_land_cold, TILE_INDEX_NAME);
      }

      idx_cold = (int *)malloc(nidx_cold * sizeof(int));//ML
      vid = mpp_get_varid(fid_land_cold, TILE_INDEX_NAME);
      mpp_get_var_value(fid_land_cold, vid, idx_cold);

      /* get destination grid */
      read_mosaic_grid_data(dst_mosaic, "x", nx_dst, ny_dst, x_tmp, face_dst, 1, 1);
      read_mosaic_grid_data(dst_mosaic, "y", nx_dst, ny_dst, y_tmp, face_dst, 1, 1);

      for (i = 0; i < nx_dst; i++) lon_axis_dst[i] = x_tmp[i];
      for (i = 0; i < ny_dst; i++) lat_axis_dst[i] = y_tmp[i * nx_dst];

      fid_dst = mpp_open(file_dst, MPP_WRITE);

      nlon = mpp_get_dimlen(fid_cold, LON_NAME);
      nlat = mpp_get_dimlen(fid_cold, LAT_NAME);
      if (nx_dst != nlon)
        mpp_error("remap_land: mismatch on the longitude dimension size "
                  "between destination mosaic grid file and dst_cold_restart");
      if (ny_dst != nlat)
        mpp_error("remap_land: mismatch on the latitude dimension size "
                  "between destination mosaic grid file and dst_cold_restart");

      if (filetype != LANDTYPE) {
        nlon = mpp_get_dimlen(fid_land_cold, LON_NAME);
        nlat = mpp_get_dimlen(fid_land_cold, LAT_NAME);
        if (nx_dst != nlon)
          mpp_error("remap_land: mismatch on the longitude dimension size "
                    "between destination mosaic grid file and land_cold_restart");
        if (ny_dst != nlat)
          mpp_error("remap_land: mismatch on the latitude dimension size "
                    "between destination mosaic grid file and land_cold_restart");
      }

      /* get destination grid */
      for (i = isc_dst; i <= iec_dst; i++) {
        x_dst[i - isc_dst] = x_tmp[i] * D2R;
        y_dst[i - isc_dst] = y_tmp[i] * D2R;
      }

      for (i = 0; i < nxc_dst; i++) land_count_dst[i] = 0;
      for (i = 0; i < ntile_dst * nxc_dst; i++) land_idx_map[i] = -1;

      {
        frac_cold = (double *)malloc(nidx_cold * sizeof(double));
        vid = mpp_get_varid(fid_land_cold, FRAC_NAME);
        mpp_get_var_value(fid_land_cold, vid, frac_cold);

        for (i = 0; i < nxc_dst; i++) {
          soil_frac_cold[i] = 0.0;
        }

        for (i = 0; i < ntile_dst * nxc_dst; i++) {
          soil_tag_dst[i] = MPP_FILL_INT;
          vegn_tag_dst[i] = MPP_FILL_INT;
          land_frac_dst[i] = MPP_FILL_DOUBLE;
        }

        if (has_glac) {
          for (i = 0; i < ntile_dst * nxc_dst; i++) glac_tag_dst[i] = MPP_FILL_INT;
        }

        if (has_lake) {
          for (i = 0; i < ntile_dst * nxc_dst; i++) lake_tag_dst[i] = MPP_FILL_INT;
        }

        /* soil, vegn */
        get_land_tile_info(fid_land_cold, SOIL_NAME, VEGN_NAME, nidx_cold, idx_cold, frac_cold, nx_dst, ny_dst,
                           ntile_cold, isc_dst, iec_dst, soil_count_cold, tmp_frac_cold, soil_tag_cold, vegn_tag_cold,
                           NULL, 0);

        /* glac */
        if (has_glac)
          get_land_tile_info(fid_land_cold, GLAC_NAME, NULL, nidx_cold, idx_cold, frac_cold, nx_dst, ny_dst, 1, isc_dst,
                             iec_dst, glac_count_cold, glac_frac_cold, glac_tag_cold, NULL, NULL, 0);

        /* lake */
        if (has_lake)
          get_land_tile_info(fid_land_cold, LAKE_NAME, NULL, nidx_cold, idx_cold, frac_cold, nx_dst, ny_dst, 1, isc_dst,
                             iec_dst, lake_count_cold, lake_frac_cold, lake_tag_cold, NULL, NULL, 0);

        for (i = 0; i < nxc_dst; i++)
          for (l = 0; l < soil_count_cold[i]; l++) soil_frac_cold[i] += tmp_frac_cold[ntile_cold * i + l];
      }



      /* find nearest source grid for each destination grid for soil, lake, glac
       */
      /*--------------------------------------------------------------------------------
        if remap_file exists, read the remap file, otherwise
        Find the remap index
        -------------------------------------------------------------------------------*/
      {
        int remap_file_exist;
        int write_remap_file;
        char file[512];
        remap_file_exist = 0;
        write_remap_file = 0;
        if (remap_file) {
          get_actual_file_name(nface_dst, face_dst, remap_file, file);
          remap_file_exist = mpp_file_exist(file);
          if (!remap_file_exist) write_remap_file = 1;
        }
        /* The following mpp_sync is necessary. It is possible that the root pe
           is in the process of writing remap_file and other processors think
           the remap_file exist
        */
        mpp_sync();
        if (remap_file_exist) { /* read from remap file */
          int fid;
          // int nidx

          if (mpp_pe() == mpp_root_pe()) printf("Read remap information from remap_file \n");
          for (l = 0; l < 4; l++) {
            start[l] = 0;
            nread[l] = 1;
          }
          start[0] = isc_dst;
          nread[0] = nxc_dst;
          fid = mpp_open(file, MPP_READ);
          nlon = mpp_get_dimlen(fid, LON_NAME);
          nlat = mpp_get_dimlen(fid, LAT_NAME);
          if (nlon != nx_dst)
            mpp_error("remap_land: mismatch of dimension length of LON between "
                      "destination grid and remap_file");
          if (nlat != ny_dst)
            mpp_error("remap_land: mismatch of dimension length of LAT between "
                      "destination grid and remap_file");

          vid = mpp_get_varid(fid, "remap_face_soil");
          mpp_get_var_value_block(fid, vid, start, nread, face_map_soil);
          vid = mpp_get_varid(fid, "remap_index_soil");
          mpp_get_var_value_block(fid, vid, start, nread, idx_map_soil);

          if (has_glac) {
            vid = mpp_get_varid(fid, "remap_face_glac");
            mpp_get_var_value_block(fid, vid, start, nread, face_map_glac);
            vid = mpp_get_varid(fid, "remap_index_glac");
            mpp_get_var_value_block(fid, vid, start, nread, idx_map_glac);
          }

          if (has_lake) {
            vid = mpp_get_varid(fid, "remap_face_lake");
            mpp_get_var_value_block(fid, vid, start, nread, face_map_lake);
            vid = mpp_get_varid(fid, "remap_index_lake");
            mpp_get_var_value_block(fid, vid, start, nread, idx_map_lake);
          }
          mpp_close(fid);
        } else {
          /* compute remap info for soil */
          full_search_nearest(nface_src, nx_src * ny_src, x_src, y_src, soil_count_src, nxc_dst, x_dst, y_dst,
                              soil_count_cold, idx_map_soil, face_map_soil);
          /* compute remap info for glac */
          if (has_glac)
            full_search_nearest(nface_src, nx_src * ny_src, x_src, y_src, glac_count_src, nxc_dst, x_dst, y_dst,
                                glac_count_cold, idx_map_glac, face_map_glac);
          /* compute remap info for lake */
          if (has_lake)
            full_search_nearest(nface_src, nx_src * ny_src, x_src, y_src, lake_count_src, nxc_dst, x_dst, y_dst,
                                lake_count_cold, idx_map_lake, face_map_lake);
        }

        if (write_remap_file) { /* write out restart file */
          int *gdata = NULL;
          int fid;
          int dim_npts;
          int id_face_soil, id_index_soil;
          int id_face_glac, id_index_glac;
          int id_face_lake, id_index_lake;
          int dim_lon, dim_lat;

          if (mpp_pe() == mpp_root_pe()) printf("write remap information to remap_file\n");
          fid = mpp_open(file, MPP_WRITE);
          dim_lon = mpp_def_dim(fid, LON_NAME, nx_dst);
          dim_lat = mpp_def_dim(fid, LAT_NAME, ny_dst);
          dim_npts = mpp_def_dim(fid, "num_points", nx_dst * ny_dst);
          gdata = (int *)malloc(nx_dst * ny_dst * sizeof(int));

          id_face_soil =
            mpp_def_var(fid, "remap_face_soil", NC_INT, 1, &dim_npts, 1, "standard_name", "soil remap face");
          id_index_soil =
            mpp_def_var(fid, "remap_index_soil", NC_INT, 1, &dim_npts, 1, "standard_name", "soil remap index");
          if (has_glac) {
            id_face_glac =
              mpp_def_var(fid, "remap_face_glac", NC_INT, 1, &dim_npts, 1, "standard_name", "glac remap face");
            id_index_glac =
              mpp_def_var(fid, "remap_index_glac", NC_INT, 1, &dim_npts, 1, "standard_name", "glac remap index");
          }
          if (has_lake) {
            id_face_lake =
              mpp_def_var(fid, "remap_face_lake", NC_INT, 1, &dim_npts, 1, "standard_name", "lake remap face");
            id_index_lake =
              mpp_def_var(fid, "remap_index_lake", NC_INT, 1, &dim_npts, 1, "standard_name", "lake remap index");
          }
          mpp_def_global_att(fid, "history", history);

          mpp_end_def(fid);

          mpp_gather_field_int_root(nxc_dst, face_map_soil, gdata);
          mpp_put_var_value(fid, id_face_soil, gdata);
          mpp_gather_field_int_root(nxc_dst, idx_map_soil, gdata);
          mpp_put_var_value(fid, id_index_soil, gdata);

          if (has_glac) {
            mpp_gather_field_int_root(nxc_dst, face_map_glac, gdata);
            mpp_put_var_value(fid, id_face_glac, gdata);
            mpp_gather_field_int_root(nxc_dst, idx_map_glac, gdata);
            mpp_put_var_value(fid, id_index_glac, gdata);
          }

          if (has_lake) {
            mpp_gather_field_int_root(nxc_dst, face_map_lake, gdata);
            mpp_put_var_value(fid, id_face_lake, gdata);
            mpp_gather_field_int_root(nxc_dst, idx_map_lake, gdata);
            mpp_put_var_value(fid, id_index_lake, gdata);
          }
          mpp_close(fid);
          free(gdata);
        }
      }

      /* get soil, glac, lake and vegn data on the destination grid when
       * file_type is land */
      {
        nidx_dst = 0;
        for (i = 0; i < nxc_dst; i++) {
          int n, idx, p, count;

          if (soil_count_cold[i] > 0) {
            double totfrac;

            n = face_map_soil[i];
            if (n < 0) {
              printf("soil_count_cold=%d, face_map_soil=%d, pe=%d, i=%d\n", soil_count_cold[i], n, mpp_pe(), i);
              mpp_error("remap_land: soil_count_cold >0 but face_map_soil<0");
            }
            idx = idx_map_soil[i];
            p = n * nx_src * ny_src + idx;
            count = soil_count_src[p];
            if (filetype != GLACTYPE && filetype != LAKETYPE) {
              pos = land_count_dst[i];
              totfrac = 0;
              for (l = 0; l < count; l++) {
                soil_tag_dst[ntile_dst * i + pos] = soil_tag_src[ntile_src * p + l];
                vegn_tag_dst[ntile_dst * i + pos] = vegn_tag_src[ntile_src * p + l];
                land_frac_dst[ntile_dst * i + pos] = soil_frac_src[ntile_src * p + l] * soil_frac_cold[i];
                land_idx_map[ntile_dst * i + pos] = idx_soil_src[ntile_src * p + l];
                land_face_map[ntile_dst * i + pos] = n;
                totfrac += soil_frac_src[ntile_src * p + l];
                pos++;
              }
              pos = land_count_dst[i];
              for (l = pos; l < pos + count; l++) land_frac_dst[ntile_dst * i + l] /= totfrac;
              nidx_dst += count;
            }
            land_count_dst[i] += count;
          }
          if (has_glac) {
            if (glac_count_cold[i] > 0) {
              n = face_map_glac[i];
              idx = idx_map_glac[i];
              p = n * nx_src * ny_src + idx;
              pos = land_count_dst[i];
              count = glac_count_src[p];
              if (count != 1) mpp_error("remap_land: glac_count_src should be 1");
              if (filetype != SOILTYPE && filetype != VEGNTYPE && filetype != LAKETYPE) {
                glac_tag_dst[ntile_dst * i + pos] = glac_tag_src[p];
                land_frac_dst[ntile_dst * i + pos] = glac_frac_cold[i]; /* preserve fraction in cold restart file
                                                                         */
                land_idx_map[ntile_dst * i + pos] = idx_glac_src[p];
                land_face_map[ntile_dst * i + pos] = n;
                nidx_dst++;
              }
              land_count_dst[i]++;
            }
          }
          if (has_lake) {
            if (lake_count_cold[i] > 0) {
              n = face_map_lake[i];
              idx = idx_map_lake[i];
              p = n * nx_src * ny_src + idx;
              count = lake_count_src[p];
              pos = land_count_dst[i];
              if (count != 1) mpp_error("remap_land: lake_count_src should be 1");
              if (filetype != SOILTYPE && filetype != VEGNTYPE && filetype != GLACTYPE) {
                lake_tag_dst[ntile_dst * i + pos] = lake_tag_src[p];
                land_frac_dst[ntile_dst * i + pos] = lake_frac_cold[i]; /* preserve fraction in cold restart file*/
                land_idx_map[ntile_dst * i + pos] = idx_lake_src[p];
                land_face_map[ntile_dst * i + pos] = n;
                nidx_dst++;
              }
              land_count_dst[i]++;
            }
          }
        }
      }
      nidx_dst_global = nidx_dst;
      mpp_sum_int(1, &nidx_dst_global);

      /* compute the tile_index */
      for (i = 0; i < ntile_dst * nxc_dst; i++){
        idx_dst[i] = MPP_FILL_INT;
      }
      for (n = 0; n < ntile_dst; n++) {
        for (i = 0; i < nxc_dst; i++) {
          if (land_idx_map[ntile_dst * i + n] > -1) {
            idx_dst[ntile_dst * i + n] = n * nx_dst * ny_dst + i + isc_dst;
          }
        }
      }


      if (filetype == VEGNTYPE) {
        printf("Before compute cohort index1\n");
        ncoho_idx_dst = compute_dst_coho_idx
          (nxc_dst, ntile_src, ncohort, nx_src * ny_src, ntile_src, ncohort, nface_dst,
           idx_map_soil, face_map_soil, soil_count_cold, ncoho_idx_src, coho_idx_src,
           coho_idx_dst, coho_idx_d_to_s, coho_idx_face_d_to_s);
        ncoho_idx_dst_global = ncoho_idx_dst;
        printf("After compute cohort index. pos=%d\n", pos);
      }


      /* define the metadata for dst_restart_file */
      {
      	int dim_time, dim_cohort_index, dim_lat, dim_lon;
        int dim_tile_index, dim_cohort, dim_tile, dim_z;
        int dim_sc_cohort, dim_lc_cohort;
        int dim_nspecies;
        int dim_textlen;
        int dim_mdf_day;
        // int dim_sc_cohort_index,  dim_lc_cohort_index;

        dim_lon = mpp_def_dim(fid_dst, LON_NAME, nx_dst);
        dim_lat = mpp_def_dim(fid_dst, LAT_NAME, ny_dst);
        dim_tile = mpp_def_dim(fid_dst, TILE_NAME, ntile_dst);
        dim_tile_index = mpp_def_dim(fid_dst, TILE_INDEX_NAME, max(nidx_dst_global, 1));
        if (zaxis_exist) dim_z = mpp_def_dim(fid_dst, LEVEL_NAME, nz);
        if (filetype == VEGNTYPE) {
          dim_cohort = mpp_def_dim(fid_dst, COHORT_NAME, ncohort);
          dim_cohort_index = mpp_def_dim(fid_dst, COHORT_INDEX_NAME, max(ncoho_idx_dst_global, 1));
        }

        // TODO New 2020
        if ((filetype == VEGNTYPE) || (filetype == SOILTYPE)) {
          if (mpp_var_exist(fid_src[0], SC_COHORT_NAME)) {
            n_sc_cohort = mpp_get_dimlen(fid_src[0], SC_COHORT_NAME);
            dim_sc_cohort = mpp_def_dim(fid_dst, SC_COHORT_NAME, n_sc_cohort);
          }
          if (mpp_var_exist(fid_src[0], LC_COHORT_NAME)) {
            n_lc_cohort = mpp_get_dimlen(fid_src[0], LC_COHORT_NAME);
            dim_lc_cohort = mpp_def_dim(fid_dst, LC_COHORT_NAME, n_lc_cohort);
          }

          if (mpp_var_exist(fid_src[0], NSPECIES_NAME)) {
            n_nspecies = mpp_get_dimlen(fid_src[0], NSPECIES_NAME);
            dim_nspecies = mpp_def_dim(fid_dst, NSPECIES_NAME, n_nspecies);
          }

          if (mpp_var_exist(fid_src[0], TEXTLEN_NAME)) {
            n_textlen = mpp_get_dimlen(fid_src[0], TEXTLEN_NAME);
            dim_textlen = mpp_def_dim(fid_dst, TEXTLEN_NAME, n_textlen);
          }
        }
        if (filetype == FIRETYPE) {
          if (mpp_var_exist(fid_src[0], MDF_DAY_NAME)) {
            n_mdf_day = mpp_get_dimlen(fid_src[0], MDF_DAY_NAME);
            dim_mdf_day = mpp_def_dim(fid_dst, MDF_DAY_NAME, n_mdf_day);
          }
        }

        if (time_exist){
          dim_time = mpp_def_dim(fid_dst, timename, NC_UNLIMITED);
        }

        for (l = 0; l < nvar_src; l++) {
          char varname[64], dimname[64];
          int vid1, vid2, ndim, m, dims[4];

          mpp_get_varname(fid_src[0], l, varname);
          vid1 = mpp_get_varid(fid_src[0], varname);
          ndim = mpp_get_var_ndim(fid_src[0], vid1);

          printf("LR varname= %s vid1= %d ndim= %d\n", varname, vid1, ndim);

          for (m = 0; m < ndim; m++) {
            mpp_get_var_dimname(fid_src[0], vid1, m, dimname);

            if (!strcmp(dimname, timename))
              dims[m] = dim_time;
            else if (!strcmp(dimname, COHORT_INDEX_NAME))
              dims[m] = dim_cohort_index;
            else if (!strcmp(dimname, LAT_NAME))
              dims[m] = dim_lat;
            else if (!strcmp(dimname, LON_NAME))
              dims[m] = dim_lon;
            else if (!strcmp(dimname, TILE_INDEX_NAME))
              dims[m] = dim_tile_index;
            else if (!strcmp(dimname, LEVEL_NAME))
              dims[m] = dim_z;
            else if (!strcmp(dimname, COHORT_NAME))
              dims[m] = dim_cohort;
            else if (!strcmp(dimname, TILE_NAME))  // This and below are new Nov 2020
              dims[m] = dim_tile;
            else if (!strcmp(dimname, NSPECIES_NAME))
              dims[m] = dim_nspecies;
            else if (!strcmp(dimname, TEXTLEN_NAME))
              dims[m] = dim_textlen;
            else if (!strcmp(dimname, SC_COHORT_NAME))
              dims[m] = dim_sc_cohort;
            else if (!strcmp(dimname, LC_COHORT_NAME))
              dims[m] = dim_lc_cohort;
            else if (!strcmp(dimname,  MDF_DAY_NAME))
              dims[m] = dim_mdf_day;
            else {
              printf("REMAP_LAND: invalid dimension name %s ", dimname);
              mpp_error("REMAP_LAND: invalid dimension name ");
            }
          }

          vid2 = mpp_def_var(fid_dst, varname, var_type[l], ndim, dims, 0);
          mpp_copy_var_att(fid_src[0], vid1, fid_dst, vid2);
        }
        if (src_has_tile == 0 && cold_has_tile == 1) {
          int vid2, vid1;
          vid2 = mpp_def_var(fid_dst, TILE_NAME, MPP_INT, 1, &dim_tile, 0);
          vid1 = mpp_get_varid(fid_cold, TILE_NAME);
          mpp_copy_var_att(fid_cold, vid1, fid_dst, vid2);
        }
        if (src_has_cohort == 0 && cold_has_cohort == 1) {
          int vid2, vid1;
          vid2 = mpp_def_var(fid_dst, COHORT_NAME, MPP_INT, 1, &dim_cohort, 0);
          vid1 = mpp_get_varid(fid_cold, COHORT_NAME);
          mpp_copy_var_att(fid_cold, vid1, fid_dst, vid2);
        }
      }

      mpp_def_global_att(fid_dst, "history", history);

      mpp_end_def(fid_dst);

      /*-------------------------------------------------------------------------------
        Remap the data and write out to dst_restart_file
        It is assumed all the fields need to be remapped will have the first
        dimension name "tile_index" or "cohort_index" ( excludes "tile_index" or
        "cohort_index")
        -----------------------------------------------------------------------------*/

      rdata_global = (double *)malloc(nidx_dst_global * sizeof(double));
      idata_global = (int *)malloc(nidx_dst_global * sizeof(int));
      coho_idata_global = (int *)malloc(ncoho_idx_dst_global * sizeof(int));

      /* loop through each time level */
      for (t = 0; t < ntime; t++) {
        for (l = 0; l < nvar_src; l++) {
          char varname[64], dimname[64];
          int m;

          if (!has_taxis[l] && t > 0) continue;
          mpp_get_varname(fid_src[0], l, varname);

          if (nidx_dst_global == 0) {
            if (strcmp(varname, LON_NAME) && strcmp(varname, LAT_NAME) && strcmp(varname, LEVEL_NAME) &&
                strcmp(varname, TILE_NAME) && strcmp(varname, COHORT_NAME) && strcmp(varname, timename))
              continue;
          }

          vid_dst = mpp_get_varid(fid_dst, varname);
          vid_src = mpp_get_varid(fid_src[0], varname);

          if (time_exist) {
            if (strcmp(varname, timename) == 0) {
              /* copy the time data from src_restart_file to dst_restart_file */
              for (m = 0; m < 4; m++) {
                start[m] = 0;
                nwrite[m] = 1;
              }
              start[0] = t;
              mpp_put_var_value_block(fid_dst, vid_dst, start, nwrite, time_data + t);
              continue;
            }
          }

          pos = 0;
          for (m = 0; m < 4; m++) {
            start[m] = 0;
            nread[m] = 1;
          }

          if (!strcmp(varname, N_ACCUM_NAME) || !strcmp(varname, NMN_ACM_NAME)) { /* for n_accum, nmn_acm */
            int fid, scalar_data;
            fid = fid_src[0];
            if (nface_src == nface_dst) fid = fid_src[face_dst];
            vid = mpp_get_varid(fid, varname);
            mpp_get_var_value(fid, vid, &scalar_data);
            mpp_put_var_value(fid_dst, vid_dst, &scalar_data);
          } else if (!strcmp(varname, LON_NAME))
            mpp_put_var_value(fid_dst, vid_dst, lon_axis_dst);
          else if (!strcmp(varname, LAT_NAME))
            mpp_put_var_value(fid_dst, vid_dst, lat_axis_dst);
          else if (!strcmp(varname, TILE_NAME))
            mpp_put_var_value(fid_dst, vid_dst, tile_axis_data);
          else if (!strcmp(varname, LEVEL_NAME))
            mpp_put_var_value(fid_dst, vid_dst, z_axis_data);
          else if (!strcmp(varname, NSPECIES_NAME))
            mpp_put_var_value(fid_dst, vid_dst, nspecies_data);
          else if (!strcmp(varname, TEXTLEN_NAME))
            mpp_put_var_value(fid_dst, vid_dst, textlen_data);
          else if (!strcmp(varname, SC_COHORT_NAME)) {
            mpp_put_var_value(fid_dst, vid_dst, sc_cohort_data);
          } else if (!strcmp(varname, LC_COHORT_NAME)) {
            mpp_put_var_value(fid_dst, vid_dst, lc_cohort_data);
          } else if (!strcmp(varname, MDF_DAY_NAME)) {
            mpp_put_var_value(fid_dst, vid_dst, mdf_day_data);
          } else if (!strcmp(varname, COHORT_NAME)) {
            mpp_put_var_value(fid_dst, vid_dst, cohort_data);
          } else if (!strcmp(varname, TILE_INDEX_NAME)) {
            compress_int_data(ntile_dst, nxc_dst, nidx_dst, nidx_dst_global, land_count_dst, idx_dst, idata_global,
                              use_all_tile);
            mpp_put_var_value(fid_dst, vid_dst, idata_global);
          } else if (!strcmp(varname, COHORT_INDEX_NAME)) {
            printf("Before compressing cohort index\n");
            // compress_int_data(ncohort * ntile_dst, nxc_dst, ncoho_idx_dst, ncoho_idx_dst_global, land_count_dst,
            //coho_idx_dst,  coho_idata_global, use_all_tile);
            printf("Before saving cohort index\n");
            mpp_put_var_value(fid_dst, vid_dst, coho_idx_dst);
            // mpp_put_var_value(fid_dst, vid_dst, coho_idata_global);
            printf("After saving cohort index\n");
          } else if (!strcmp(varname, FRAC_NAME)) {
            compress_double_data(ntile_dst, nxc_dst, nidx_dst, nidx_dst_global, land_count_dst, land_frac_dst,
                                 rdata_global, use_all_tile);
            mpp_put_var_value(fid_dst, vid_dst, rdata_global);
          } else if (!strcmp(varname, GLAC_NAME)) {
            if (has_glac) {
              compress_int_data(ntile_dst, nxc_dst, nidx_dst, nidx_dst_global, land_count_dst, glac_tag_dst,
                                idata_global, use_all_tile);
              mpp_put_var_value(fid_dst, vid_dst, idata_global);
            }
          } else if (!strcmp(varname, LAKE_NAME)) {
            if (has_lake) {
              compress_int_data(ntile_dst, nxc_dst, nidx_dst, nidx_dst_global, land_count_dst, lake_tag_dst,
                                idata_global, use_all_tile);
              mpp_put_var_value(fid_dst, vid_dst, idata_global);
            }
          } else if (!strcmp(varname, SOIL_NAME)) {
            compress_int_data(ntile_dst, nxc_dst, nidx_dst, nidx_dst_global, land_count_dst, soil_tag_dst, idata_global,
                              use_all_tile);
            mpp_put_var_value(fid_dst, vid_dst, idata_global);
          } else if (!strcmp(varname, VEGN_NAME)) {
            compress_int_data(ntile_dst, nxc_dst, nidx_dst, nidx_dst_global, land_count_dst, vegn_tag_dst, idata_global,
                              use_all_tile);
            mpp_put_var_value(fid_dst, vid_dst, idata_global);

          } else if (!strcmp(varname, SPECIES_NAMES)) {
            // TODO: Note assumption that SPECIES_NAMES is the only character data.
            //  this will likely need to be removed in future versions

            // Data to be copied will be in  species_names[nspecies][textlen]
            int nd_nspecies = mpp_get_dimlen(fid_src[0], NSPECIES_NAME);
            int nd_textlen = mpp_get_dimlen(fid_src[0], TEXTLEN_NAME);
            char *vdata_src = (char *)malloc(nd_textlen * nd_nspecies * sizeof(char));
            size_t cstart[2], cnread[2];
            cstart[0] = 0;
            cnread[0] = nd_nspecies;
            cstart[1] = 0;
            cnread[1] = nd_textlen;
            mpp_get_var_value_block(fid_src[0], vid_src, cstart, cnread, vdata_src);
            mpp_put_var_value_block(fid_dst, vid_dst, cstart, cnread, vdata_src);
            free(vdata_src);
            vdata_src = NULL;
          } else if (nz_src[l] >= 0 && ((ndim_src[l] == 2) || ((ndim_src[l] == 3) && has_taxis[l]))) {
            /**
             * Read source data and do remapping for other fiels of original lm4.0 type.
             * This is the original lm4.0 code. The dimensions of the fields are either
             * [zfull, tile_index] of [time,zfull,tile_index]
             **/
            printf("\n*RL Other field lm4.0 varnameO=%s nface_src=%d nz_src[l]=%d\n", varname, nface_src, nz_src[l]);

            start_pos[0] = 0;
            for (m = 1; m < nface_src; m++) {
              start_pos[m] = start_pos[m - 1] + nidx_src[m - 1];
            }

            for (k = 0; k < nz_src[l]; k++) {
              pos = 0;
              int kid;
              /* read the source data */
              for (m = 0; m < nface_src; m++) {
                kid = 0;
                if (has_taxis[l]) {
                  start[0] = t;
                  kid = 1;
                }

                start[kid] = k;
                start[ndim_src[l] - 1] = 0;
                nread[ndim_src[l] - 1] = nidx_src[m];
                vid_src = mpp_get_varid(fid_src[m], varname);
                // printf("*RL C VBB varname=%s vid%d vartype=%d k=%d nz_src[l]=%d start=%lu, nread=%lu has_taxis=%d\n",
                //      varname, vid_src, var_type[l], k, nz_src[l], start[0], nread[ndim_src[l] - 1], has_taxis[l]);
                if (var_type[l] == MPP_INT) {
                  mpp_get_var_value_block(fid_src[m], vid_src, start, nread, idata_src + pos);
                } else if (var_type[l] == MPP_DOUBLE) {
                  mpp_get_var_value_block(fid_src[m], vid_src, start, nread, data_src + pos);
                } else {
                  mpp_error("remap_land : reading block for vartype other than INT or DOUBLE");
                }

                pos += nidx_src[m];
              }  // m loop
              for (m = 0; m < 4; m++) {
                start[m] = 0;
                nwrite[m] = 1;
              }
              kid = 0;
              if (has_taxis[l]) {
                start[0] = t;
                kid = 1;
              }

              start[kid] = k;
              start[ndim_src[l] - 1] = 0;
              nwrite[ndim_src[l] - 1] = nidx_dst_global;

              if (var_type[l] == MPP_INT) {
                for (m = 0; m < ntile_dst * nxc_dst; m++) {
                  int face, lll;
                  if (land_idx_map[m] < 0) {
                    idata_dst[m] = MPP_FILL_INT;
                    continue;
                  }
                  face = land_face_map[m];
                  lll = start_pos[face] + land_idx_map[m];
                  idata_dst[m] = idata_src[lll];
                }
                compress_int_data(ntile_dst, nxc_dst, nidx_dst, nidx_dst_global, land_count_dst, idata_dst,
                                  idata_global, use_all_tile);
                mpp_put_var_value_block(fid_dst, vid_dst, start, nwrite, idata_global);
              } else if (var_type[l] == MPP_DOUBLE) {
                for (m = 0; m < ntile_dst * nxc_dst; m++) {
                  int face, lll;
                  if (land_idx_map[m] < 0) {
                    data_dst[m] = MPP_FILL_DOUBLE;
                    continue;
                  }
                  face = land_face_map[m];
                  lll = start_pos[face] + land_idx_map[m];
                  data_dst[m] = data_src[lll];
                }
                compress_double_data(ntile_dst, nxc_dst, nidx_dst, nidx_dst_global, land_count_dst, data_dst,
                                     rdata_global, use_all_tile);
                mpp_put_var_value_block(fid_dst, vid_dst, start, nwrite, rdata_global);
              } else {
                mpp_error("remap_land : writing block for vartype other than INT or DOUBLE");
              }
            }  // k < nz_src[l]
          }    // end else other lm4.0 fields
          else if (ndim_src[l] == 1 && has_coho_idx[l] == 1) {  //COHO_IDX
            vid_src = mpp_get_varid(fid_src[0], varname);
            mpp_get_var_dimname(fid_src[0], vid_src, 0, dimname);
            if( strcmp(dimname, COHORT_INDEX_NAME)) {
              mpp_error("remap_land : 1D field but bot expected cohort_index :");
            }

            start_pos[0] = 0;
            for (m = 1; m < nface_src; m++) {
              start_pos[m] = start_pos[m - 1] + ncoho_idx_src[m - 1];
            }

            for (m = 0; m < 4; m++) {
              start[m] = 0; nread[m] = 1;  nwrite[m] = 1;
            }

            // Read all the data on all faces of current field.
            pos = 0;
            for (m = 0; m < nface_src; m++) {//for each face ?
              nread[0] = ncoho_idx_src[m]; //size of coho_index dim

              vid_src = mpp_get_varid(fid_src[m], varname);
              if (var_type[l] == MPP_INT) {
                mpp_get_var_value_block(fid_src[m], vid_src, start, nread, coho_idata_src + pos);
              } else if (var_type[l] == MPP_DOUBLE) {
                mpp_get_var_value_block(fid_src[m], vid_src, start, nread, coho_rdata_src + pos);
              } else {
                mpp_error("remap_land : reading block for vartype other than INT or DOUBLE");
              }
              pos += ncoho_idx_src[m];
            }  // m loop

            // Prepare for writing collected data
            nwrite[0] = ncoho_idx_dst_global;

            //Init the entire dst array
            if (var_type[l] == MPP_INT) {
              for (int k = 0; k < ncohort * ntile_dst * nxc_dst ; k++) {
                coho_idata_dst[k] = MPP_FILL_INT;
              }
            }else{
              for (int k = 0; k < ncohort * ntile_dst * nxc_dst  ; k++) {
                coho_rdata_dst[k] = MPP_FILL_DOUBLE;
              }
            }

            //Copy-map the data for the destination face
            for (int k = 0; k < ncoho_idx_dst_global ; k++) {
              //int i_dst = coho_idx_dst[k]; // the value of kth valid coho index
                int i_src = coho_idx_d_to_s [k]; // the corresponding source index
                int face = coho_idx_face_d_to_s [k]; // corresponding face in  src;
                int ipf_src = start_pos[face] + i_src;
                if (var_type[l] == MPP_INT) {
                  coho_idata_dst[k] = coho_idata_src[ipf_src];
               }else{
                  coho_rdata_dst[k] = coho_rdata_src[ipf_src];
                }
            }

            //Write the data
            if (var_type[l] == MPP_INT) {
              //TODO: verifiy its "compressed"
              mpp_put_var_value_block(fid_dst, vid_dst, start, nwrite, coho_idata_dst);
            }else{
              mpp_put_var_value_block(fid_dst, vid_dst, start, nwrite, coho_rdata_dst);
            }
          }
          else if (ndim_src[l] == 1) {
            vid_src = mpp_get_varid(fid_src[0], varname);
            mpp_get_var_dimname(fid_src[0], vid_src, 0, dimname);
            if (strcmp(dimname, TILE_INDEX_NAME)) {
              mpp_error("remap_land : 1D field but not tile_index");
            }

            start_pos[0] = 0;
            for (m = 1; m < nface_src; m++) {
              start_pos[m] = start_pos[m - 1] + nidx_src[m - 1];
            }

            for (m = 0; m < 4; m++) {
              start[m] = 0;
              nread[m] = 1;
              nwrite[m] = 1;
            }

            // Collect all the data
            pos = 0;
            for (m = 0; m < nface_src; m++) {
              nread[0] = nidx_src[m];

              vid_src = mpp_get_varid(fid_src[m], varname);
              if (var_type[l] == MPP_INT) {
                mpp_get_var_value_block(fid_src[m], vid_src, start, nread, idata_src + pos);
              } else if (var_type[l] == MPP_DOUBLE) {
                mpp_get_var_value_block(fid_src[m], vid_src, start, nread, data_src + pos);
              } else {
                mpp_error("remap_land : reading block for vartype other than INT or DOUBLE");
              }
              pos += nidx_src[m];
            }  // m loop

            // Prepare for writing collected data
            nwrite[0] = nidx_dst_global;

            if (var_type[l] == MPP_INT) {
              for (m = 0; m < ntile_dst * nxc_dst; m++) {
                int face, lll;
                if (land_idx_map[m] < 0) {
                  idata_dst[m] = MPP_FILL_INT;
                  continue;
                }
                face = land_face_map[m];
                lll = start_pos[face] + land_idx_map[m];
                idata_dst[m] = idata_src[lll];
              }
              compress_int_data(ntile_dst, nxc_dst, nidx_dst, nidx_dst_global, land_count_dst, idata_dst, idata_global,
                                use_all_tile);
              mpp_put_var_value_block(fid_dst, vid_dst, start, nwrite, idata_global);
            } else if (var_type[l] == MPP_DOUBLE) {
              for (m = 0; m < ntile_dst * nxc_dst; m++) {
                int face, lll;
                if (land_idx_map[m] < 0) {
                  data_dst[m] = MPP_FILL_DOUBLE;
                  continue;
                }
                face = land_face_map[m];
                lll = start_pos[face] + land_idx_map[m];
                data_dst[m] = data_src[lll];
              }
              compress_double_data(ntile_dst, nxc_dst, nidx_dst, nidx_dst_global, land_count_dst, data_dst,
                                   rdata_global, use_all_tile);
              mpp_put_var_value_block(fid_dst, vid_dst, start, nwrite, rdata_global);
            }
          }else if (ndim_src[l] == 2) {
            /**
             * The fields dimension should be [lc_cohort, tile_index] or [mdf_day,tile_inde]
             */
            char dimnameA[64], dimnameB[64];
            int lvid = mpp_get_varid(fid_src[0], varname);
            mpp_get_var_dimname(fid_src[0], lvid, 0, dimnameA);
            mpp_get_var_dimname(fid_src[0], lvid, 1, dimnameB);

            if (strcmp(dimnameB, TILE_INDEX_NAME) || (strcmp(dimnameA, LC_COHORT_NAME) && strcmp(dimnameA, MDF_DAY_NAME)) ) {
              mpp_error("remap_land : Expected 2D field with dims [lc_cohort,tile_index] or [mdf_name,tile_index] :");
            }

            start_pos[0] = 0;
            for (m = 1; m < nface_src; m++) {
              start_pos[m] = start_pos[m - 1] + nidx_src[m - 1];
            }

            int ldimlen = -1;
             if(!strcmp(dimnameA, LC_COHORT_NAME)){
               ldimlen = mpp_get_dimlen(fid_src[0], LC_COHORT_NAME);
             }else {
               ldimlen = mpp_get_dimlen(fid_src[0], MDF_DAY_NAME);
            }

            for (k = 0; k < ldimlen; k++) {
              pos = 0;
              // int kid;
              /* read the source data */
              for (m = 0; m < nface_src; m++) {
                // kid = 0;
                start[0] = k;
                nread[0] = 1;
                start[1] = 0;
                nread[1] = nidx_src[m];
                vid_src = mpp_get_varid(fid_src[m], varname);
                if (var_type[l] == MPP_INT) {
                  mpp_get_var_value_block(fid_src[m], vid_src, start, nread, idata_src + pos);
                } else if (var_type[l] == MPP_DOUBLE) {
                  mpp_get_var_value_block(fid_src[m], vid_src, start, nread, data_src + pos);
                } else {
                  mpp_error("remap_land : reading block for vartype other than INT or DOUBLE");
                }
                pos += nidx_src[m];
              }  // m loop

              // And write the data
              for (m = 0; m < 4; m++) {
                start[m] = 0;
                nwrite[m] = 1;
              }

              start[0] = k;
              nwrite[0] = 1;
              start[1] = 0;
              nwrite[1] = nidx_dst_global;

              if (var_type[l] == MPP_INT) {
                for (m = 0; m < ntile_dst * nxc_dst; m++) {
                  int face, lll;
                  if (land_idx_map[m] < 0) {
                    idata_dst[m] = MPP_FILL_INT;
                    continue;
                  }
                  face = land_face_map[m];
                  lll = start_pos[face] + land_idx_map[m];
                  idata_dst[m] = idata_src[lll];
                }
                compress_int_data(ntile_dst, nxc_dst, nidx_dst, nidx_dst_global, land_count_dst, idata_dst,
                                  idata_global, use_all_tile);
                mpp_put_var_value_block(fid_dst, vid_dst, start, nwrite, idata_global);
              } else if (var_type[l] == MPP_DOUBLE) {
                for (m = 0; m < ntile_dst * nxc_dst; m++) {
                  int face, lll;
                  if (land_idx_map[m] < 0) {
                    data_dst[m] = MPP_FILL_DOUBLE;
                    continue;
                  }
                  face = land_face_map[m];
                  lll = start_pos[face] + land_idx_map[m];
                  data_dst[m] = data_src[lll];
                }
                compress_double_data(ntile_dst, nxc_dst, nidx_dst, nidx_dst_global, land_count_dst, data_dst,
                                     rdata_global, use_all_tile);
                mpp_put_var_value_block(fid_dst, vid_dst, start, nwrite, rdata_global);
              } else {
                mpp_error("remap_land : writing block for vartype other than INT or DOUBLE");
              }
            }  // k < nz_src[l]

          } else if (ndim_src[l] == 3) {
            /**
             *  The fileds dimensions should be [soilCCohort,zfull,tile_index]
             **/
            // TODO:Simplify Combine sections that use tile_index as special cases of 3D ?
            char dimnameA[64], dimnameB[64], dimnameC[64];
            int lvid = mpp_get_varid(fid_src[0], varname);
            mpp_get_var_dimname(fid_src[0], lvid, 0, dimnameA);
            mpp_get_var_dimname(fid_src[0], lvid, 1, dimnameB);
            mpp_get_var_dimname(fid_src[0], lvid, 2, dimnameC);
            printf("LR 3D varname dimA dimB dimC %s %s %s %s\n", varname, dimnameA, dimnameB, dimnameC);

            if (strcmp(dimnameA, SC_COHORT_NAME) || strcmp(dimnameB, LEVEL_NAME) || strcmp(dimnameC, TILE_INDEX_NAME)) {
              mpp_error("remap_land : Expected 3D  field with dims [soilCCohort,zfull,tile_index]");
            }

            int sc_dimlen = mpp_get_dimlen(fid_src[0], SC_COHORT_NAME);
            start_pos[0] = 0;

            for (m = 1; m < nface_src; m++) {
              start_pos[m] = start_pos[m - 1] + nidx_src[m - 1];
            }

            for (m = 0; m < 4; m++) {
              start[m] = 0;
              nread[m] = 1;
              nwrite[m] = 1;
            }

            for (int is = 0; is < sc_dimlen; is++) {
              for (k = 0; k < nz_src[l]; k++) {
                pos = 0;
                for (m = 0; m < nface_src; m++) {//TODO: check start with 1 or 0 ?
                  start[0] = is;
                  start[1] = k;
                  start[2] = 0;
                  nread[2] = nidx_src[m];
                  vid_src = mpp_get_varid(fid_src[m], varname);
                  if (var_type[l] == MPP_INT) {
                    mpp_get_var_value_block(fid_src[m], vid_src, start, nread, idata_src + pos);
                  } else if (var_type[l] == MPP_DOUBLE) {
                    mpp_get_var_value_block(fid_src[m], vid_src, start, nread, data_src + pos);
                  } else {
                    mpp_error("remap_land : reading block for vartype other than INT or DOUBLE");
                  }
                  pos += nidx_src[m];
                }  // m loop

                for (m = 0; m < 4; m++) {
                  start[m] = 0;
                  nwrite[m] = 1;
                }

                start[0] = is;
                start[1] = k;
                start[2] = 0;
                nwrite[2] = nidx_dst_global;

                if (var_type[l] == MPP_INT) {
                  for (m = 0; m < ntile_dst * nxc_dst; m++) {
                    int face, lll;
                    if (land_idx_map[m] < 0) {
                      idata_dst[m] = MPP_FILL_INT;
                      continue;
                    }
                    face = land_face_map[m];
                    lll = start_pos[face] + land_idx_map[m];
                    idata_dst[m] = idata_src[lll];
                  }
                  compress_int_data(ntile_dst, nxc_dst, nidx_dst, nidx_dst_global, land_count_dst, idata_dst,
                                    idata_global, use_all_tile);
                  mpp_put_var_value_block(fid_dst, vid_dst, start, nwrite, idata_global);
                } else if (var_type[l] == MPP_DOUBLE) {
                  for (m = 0; m < ntile_dst * nxc_dst; m++) {
                    int face, lll;
                    if (land_idx_map[m] < 0) {
                      data_dst[m] = MPP_FILL_DOUBLE;
                      continue;
                    }
                    face = land_face_map[m];
                    lll = start_pos[face] + land_idx_map[m];
                    data_dst[m] = data_src[lll];
                  }
                  compress_double_data(ntile_dst, nxc_dst, nidx_dst, nidx_dst_global, land_count_dst, data_dst,
                                       rdata_global, use_all_tile);
                  mpp_put_var_value_block(fid_dst, vid_dst, start, nwrite, rdata_global);
                } else {
                  mpp_error("remap_land : writing block for vartype other than INT or DOUBLE");
                }
              } // k < nz_src[l]
            } // is < ls_dimlen
          } // end  (ndim_src[l] == 3)
        }  // end loop (l = 0; l < nvar_src; l++)

        if (t == 0 && src_has_tile == 0 && cold_has_tile == 1) {
          int *tile_data_cold = NULL;
          tile_data_cold  = (int *)malloc(ntile_dst * sizeof(int));
          for (int ic = 0; ic < ntile_dst; ic++) tile_data_cold[ic] = ic + 1;
          vid_dst = mpp_get_varid(fid_dst, TILE_NAME);
          mpp_put_var_value(fid_dst, vid_dst, tile_data_cold);
          free(tile_data_cold);
        }
        if (t == 0 && src_has_cohort == 0 && cold_has_cohort == 1) {
          int *cohort_data_cold = NULL;
          cohort_data_cold = (int *)malloc(ncohort * sizeof(int));
          for (int ic = 0; ic < ncohort; ic++) cohort_data[ic] = ic + 1;
          vid_dst = mpp_get_varid(fid_dst, COHORT_NAME);
          mpp_put_var_value(fid_dst, vid_dst, cohort_data_cold);
          free(cohort_data_cold);
        }
      }  // for ntime

      mpp_close(fid_dst);
      mpp_close(fid_cold);

      free(frac_cold);
      free(rdata_global);
      free(idata_global);
      free(coho_idata_global);

      if (mpp_pe() == mpp_root_pe()) printf("NOTE from remap_land: %s is created\n", file_dst);
      if (print_memory) {
        char mesg[128];
        sprintf(mesg, "End of loop face_dst = %d\n", face_dst);
        print_mem_usage(mesg);
      }
    }  // for (face_dst = 0; face_dst < nface_dst; face_dst++) {

    free(lon_axis_dst);
    free(lat_axis_dst);
    free(x_tmp);
    free(y_tmp);
    free(start_pos);
    free(idata_src);
    free(data_src);
    free(data_dst);
    free(idata_dst);
    free(x_dst);
    free(y_dst);
    free(land_idx_map);
    free(land_face_map);
    free(land_count_dst);
    free(idx_dst);
    free(soil_tag_dst);
    free(vegn_tag_dst);
    free(land_frac_dst);
    free(soil_count_cold);
    free(soil_frac_cold);
    free(tmp_frac_cold);
    free(soil_tag_cold);
    free(vegn_tag_cold);
    free(idx_map_soil);
    free(face_map_soil);
    free(tile_axis_data);
    free(nidx_land_src);
    free(soil_count_src);
    free(soil_frac_src);
    free(soil_tag_src);
    free(vegn_tag_src);
    free(idx_soil_src);

    free(coho_idx_dst);
    free(coho_idx_d_to_s);
    free(coho_idx_face_d_to_s);
    free(coho_idata_src);
    free(coho_rdata_src);
    free(coho_idata_dst);
    free(coho_rdata_dst);

    if (has_glac) {
      free(glac_tag_dst);
      free(glac_count_cold);
      free(glac_frac_cold);
      free(glac_tag_cold);
      free(idx_map_glac);
      free(face_map_glac);
    }
    if (has_lake) {
      free(lake_tag_dst);
      free(lake_count_cold);
      free(lake_frac_cold);
      free(lake_tag_cold);
      free(idx_map_lake);
      free(face_map_lake);
    }

    mpp_delete_domain2d(&Dom_dst);
  }  // block of pre-loop though each face

  {
    int n;
    for (n = 0; n < nface_src; n++) mpp_close(fid_src[n]);
  }

  /* release memory */
  free(x_src);
  free(y_src);

  free(fid_src);
  free(nidx_src);

  free(has_taxis);
  free(ndim_src);
  free(nz_src);
  free(var_type);
  if (time_exist) free(time_data);

  if (mpp_pe() == mpp_root_pe()) printf("\n******** Successfully run remap_land***********\n");

  mpp_end();

  return 0;
}  // main

void get_actual_file_name(int nface, int face, const char *file_orig, char *file) {
  if (nface == 1)
    strcpy(file, file_orig);
  else
    sprintf(file, "%s.tile%d.nc", file_orig, face + 1);
}

/********************************************************************************
void get_land_tile_info( )
MZ: determine count[], frac[], tag1[], tax2[] and idx[] for a given
********************************************************************************/
void get_land_tile_info(int fid, const char *name1, const char *name2, int nidx, const int *idx_in,
                        const double *frac_in, int nx, int ny, int ntile, int isc, int iec, int *count, double *frac,
                        int *tag1, int *tag2, int *idx, int all_tile) {
  int vid, l, i, j, k, p, nxc, pos;
  int *tmp1 = NULL;
  int *tmp2 = NULL;

  nxc = iec - isc + 1;

  tmp1 = (int *)malloc(nidx * sizeof(int));

  vid = mpp_get_varid(fid, name1);
  mpp_get_var_value(fid, vid, tmp1);
  if (tag2) {
    if (!name2) mpp_error("remap_land: name2 can not be NULL when tag2 is not NULL");
    tmp2 = (int *)malloc(nidx * sizeof(int));
    vid = mpp_get_varid(fid, name2);
    mpp_get_var_value(fid, vid, tmp2);
  }

  /* set count to 0 */
  for (i = 0; i < nxc; i++) count[i] = 0;
  for (i = 0; i < ntile * nxc; i++) {
    frac[i] = MPP_FILL_DOUBLE;
    tag1[i] = MPP_FILL_INT;
    if (tag2) tag2[i] = MPP_FILL_INT;
    if (idx) idx[i] = MPP_FILL_INT;
  }
  pos = 0;
  for (l = 0; l < nidx; l++) {
    if (tmp1[l] != MPP_FILL_INT) {
      i = idx_in[l] % nx;
      k = idx_in[l] / nx;
      j = k % ny;
      //if (i == 48 && j == 94) printf ...
      p = j * nx + i;
      if (p < isc || p > iec) continue;
      p = p - isc;
      if (count[p] > ntile)
        mpp_error("remap_land: number of tiles is greater than allowed ntiles on one grid cell");
      frac[ntile * p + count[p]] = frac_in[l];
      tag1[ntile * p + count[p]] = tmp1[l];
      if (tag2) tag2[ntile * p + count[p]] = tmp2[l];
      if (idx) {
        if (all_tile)
          idx[ntile * p + count[p]] = l;
        else
          idx[ntile * p + count[p]] = pos;
      }
      pos++;
      count[p]++;
    }
  }

  free(tmp1);
  if (tmp2) free(tmp2);
}


/********************************************************************
 void full_search_nearest
 search the nearest point from the first of source grid to the last.
MZ: For each destination grid point, find the nearest non-masked point of the
source grid , on any face of the source grid. The resultant idx_map[] and
face_map[] have the space index and the face index of the nearest.
********************************************************************/
void full_search_nearest(int nface_src, int npts_src, const double *lon_src, const double *lat_src, const int *mask_src,
                         int npts_dst, const double *lon_dst, const double *lat_dst, const int *mask_dst, int *idx_map,
                         int *face_map) {
  int i_dst, i_src, l, face_cur;
  int m, ind_cur; //pos
  double d_cur, d;
  double p1[2], p2[2];

  for (i_dst = 0; i_dst < npts_dst; i_dst++) {
    if (mask_dst[i_dst] == 0) {
      face_map[i_dst] = -1;
      idx_map[i_dst] = -1;
      continue;
    }
    d_cur = -1;
    //pos = 0;
    p1[0] = lon_dst[i_dst];
    p1[1] = lat_dst[i_dst];

    for (m = 0; m < nface_src; m++) {
      for (i_src = 0; i_src < npts_src; i_src++) {
        l = m * npts_src + i_src;
        /*
          if(face_cur==2 && i_dst == 16*192+168 && m==2 && (i_src==185 ||i_src==186 ||i_src==187||
          i_src==233 ||i_src==234 ||i_src==235)) {
          printf("at i_src=%d, mask_src=%d\n", i_src, mask_src[l]);
          }
        */
        if (mask_src[l] == 0) continue;
        p2[0] = lon_src[l];
        p2[1] = lat_src[l];
        d = great_circle_distance(p1, p2);
        /*
          if(face_cur==2 && i_dst == 16*192+168 && m==2 && (i_src==185 ||i_src==186 ||i_src==187||
          i_src==233 ||i_src==234 ||i_src==235)) {
          printf("grid=%g,%g, d=%g,d_cur=%g,ind_cur=%d\n", p2[0]*R2D, p2[1]*R2D, d, d_cur, ind_cur);
          }
        */
        if (d_cur < 0 || d < d_cur) {
          ind_cur = i_src;
          face_cur = m;
          d_cur = d;
        }
      }
    }
    if (d_cur < 0)
      mpp_error("remap_land(full_search_nearest): no nearest point is found");
    else {
      /*
        if(face_cur==2 && i_dst == 16*192+168) printf("at (169,17), i_dst=%d\n", ind_cur);
        if(face_cur==2 && i_dst == 16*192+169) printf("at (170,17), i_dst=%d\n", ind_cur);
        if(face_cur==2 && i_dst == 16*192+168) printf("at (169,17), grid=%g,%g\n", p1[0]*R2D,p1[1]*R2D);
      */
      face_map[i_dst] = face_cur;
      idx_map[i_dst] = ind_cur;
    }
  }

} /* full_search_nearest */

/*-------------------------------------------------------------------------
  void compress_double_data ()
  get global compressed data
  ------------------------------------------------------------------------*/
void compress_double_data(int ntile, int npts, int nidx, int nidx_global, const int *land_count, const double *data,
                          double *data_global, int all_tile) {
  int pos1, pos2, pos, n, i, count;
  double *data_local = NULL;

  if (nidx > 0) data_local = (double *)malloc(nidx * sizeof(double));
  for (i = 0; i < nidx_global; i++) data_global[i] = MPP_FILL_DOUBLE;

  pos1 = 0;
  pos2 = 0;
  for (n = 0; n < ntile; n++) {
    pos = pos1;
    for (i = 0; i < npts; i++) {
      if (n < land_count[i]) {
        if (data[ntile * i + n] != MPP_FILL_DOUBLE) {
          data_local[pos] = data[ntile * i + n];
          if (!all_tile) pos++;
        }
        if (all_tile) pos++;
      }
    }
    count = pos - pos1;
    mpp_gather_field_double_root(count, data_local + pos1, data_global + pos2);
    pos1 = pos;
    mpp_sum_int(1, &count);
    pos2 += count;
  }

  if (nidx > 0) free(data_local);
}

/*-------------------------------------------------------------------------
  void compress_int_data ()
  get global compressed data
//MZ copy into data_global[] all values in data[] that are not MPP_FILL_
//TODO: Note initialization of data_local and data_global varies from version of double data.
 ------------------------------------------------------------------------*/
void compress_int_data(int ntile, int npts, int nidx, int nidx_global, const int *land_count, const int *data,
                       int *data_global, int all_tile) {
  int pos1, pos2, pos, n, i, count;
  int *data_local = NULL;

  if (nidx > 0) data_local = (int *)malloc(nidx * sizeof(int));
  for (i = 0; i < nidx; i++) data_local[i] = MPP_FILL_INT;

  pos1 = 0;
  pos2 = 0;
  for (n = 0; n < ntile; n++) {
    pos = pos1;
    for (i = 0; i < npts; i++) {
      if (n < land_count[i]) {
        if (data[ntile * i + n] != MPP_FILL_INT) {
          data_local[pos] = data[ntile * i + n];
          if (!all_tile) pos++;
        }
        if (all_tile) pos++;
      }
    }
    count = pos - pos1;
    mpp_gather_field_int_root(count, data_local + pos1, data_global + pos2);

    pos1 = pos;
    mpp_sum_int(1, &count);
    pos2 += count;
  }
  if (nidx > 0) free(data_local);
}

/*
  Search for the key in the array, returning the location where it
  was found, and -1 if it was not found.
  Search starts at position key.
  Search stops when negative values are encounterd.
 */
int key_in_array_hints(int key, int* vec, int size, int hint){
  int result = -1;
  for (int i = hint; i  <size; i++){
    if(vec[i] < 0){
      break;
    }
    if(vec[i] == key){
      result = i;
      break;
    }
  }
  if(result < 0){
    for (int i = 0; i  < hint; i++){
      if(vec[i] < 0){
        break;
      }

      if(vec[i] == key){
        result = i;
        break;
      }
    }
  }

  return result;
}
/*
 Search for the key in the array, returning the location where it
was found, and -1 if it was not found.
 */
int key_in_array(int key, int* vec, int size){
  int result = -1;
  for (int i = 0; i  <size; i++){
    if(vec[i] == key){
      result = i;
      break;
    }
  }
  return result;
}


//Binary search
  int binary_search_(int k, int v[], int l, int r) {
    if (r >= l) {
        int midp = l + (r - l) / 2;

        if (v[midp] == k)
            return midp;

        if (v[midp] > k) {
          return binary_search_(k, v, l, midp - 1);
        }else{
          return binary_search_(k, v, midp + 1, r);
        }
    }
    return -1;
}
int binary_search(int k, int v[], int s) {
  return binary_search_(k,v,0,s-1);
}

void swap_vals(int *xp, int *yp)  {
    int temp = *xp;
    *xp = *yp;
    *yp = temp;
}

void bubble_sort(int* v, int n)  {

    for (int i = 0; i < n-1; i++){
      for (int j = 0; j < n-i-1; j++)  {
        if (v[j] > v[j+1]) { swap_vals(&v[j], &v[j+1]); }
}
}
}

int partition(int v[], int low, int high){
  int pivot = v[high];
  int i = (low - 1);
  for (int j = low; j <= high - 1; j++) {
    if (v[j] <= pivot) {
      i++;
      swap_vals(&(v[i]), &(v[j]));
    }
  }
  swap_vals(&(v[i + 1]), &(v[high]));
  return (i + 1);
}

//random pivot partition
int partition_rp(int v[], int low, int high)
{
  srand(time(NULL));
  int random = low + rand() % (high - low);
  swap_vals(&(v[random]), &(v[high]));
  return partition(v, low, high);
}


//quick_sort with random pivots
void quick_sort(int v[], int low, int high) {
  if (low < high)  {
    int pivot = partition_rp(v, low, high);
    quick_sort(v, low, pivot - 1);
    quick_sort(v, pivot + 1, high);
  }
}

void quick_sorti(int v[], int n) {
  quick_sort(v, 0, n-1);
}

bool is_sorted(int* v, int n){
  bool result = true;
  for (int i = 0; i < n-1; i++) {
    if(v[i+1] < v[i]){
    result = false;
    break;
    }
  }
  return result;
}

/**
   Return true if there are any duplicates in array v.
   Array v is assumed to be sorted.
 **/
bool has_dupes_sorted(int* v, int n){
  bool result = false;
  for (int i = 0; i < n-1; i++) {
    if(v[i+1] == v[i]){
      result = true;
      break;
    }
  }
  return result;
}



int make_sorted_set(int *v, int n) {
  printf("Etering make_sorted_set\n");
  //bubble_sort(v, n);
  quick_sorti(v , n);
  int* cp = (int *)malloc(n * sizeof(int));
  int ndupes = 0;
  int j = 0;
  int i = 1;
  cp[0] = v[0];
  while (i < n) {
    //If the next in V[] is not equal to the last one copied, do copy it:
    if (cp[j] != v[i]) {
      j++;
      cp[j] = v[i];
      i++;
    } else {
      ndupes++;
      i++;
    }
  }

  // copy back into v the resorted data, and without duplicates into original vec.
  for (i = 0; i < (n - ndupes); i++) {
    v[i] = cp[i];
  }
  for (i = n - ndupes; i < n; i++) {
    v[i] = MPP_FILL_INT;
  }
  return (n - ndupes);
}

/***
 * Compute the destination cohort index (coho_idx_dst) by mapping each src cohort index to a destination
 * cohort index (coho_idx_dst).
 * Cohort index dimensions are [cohort][tile][lat][lon], and the [tile][lat][lon] dimensions are already
 * mapped in idx_soil_src and face_map_soil. We call [lat][lon] together the space index (is), so we can also represent
 * its dimension as [cohort][tile][space]
 * Legend::
 *  src - source
 *  dst - destination
 *  nt -  number f tiles, or tile dimension length
 *  nc -  number of cohorts, or cohort dimension length
 *  ns -  spatial dimension length.  ns = nx ** ny
*   ic - cohort index
    it - tile index
    icts - cohort_tile_space index
    Note that array coho_idx_src has the indices for all the src faces; while coho_idx_dst
    only
   ***/
int compute_dst_coho_idx_V2(int ns_dst, int nt_dst, int nc_dst, int ns_src, int nt_src, int nc_src, int face_dst,
                         int *idx_map_soil, int *face_map_soil,  int* ncoho_idx_src, int *coho_idx_src,
                         int *coho_idx_dst, int * coho_idx_d_to_s,  int *coho_idx_face_d_to_s) {


  int n_face_coho_idx_dst = ns_dst * nt_dst * nc_dst;
  //DODO: use an expapndable data structure to grow these or ...
  for (int i = 0; i < n_face_coho_idx_dst; i++) {
    coho_idx_dst[i] = MPP_FILL_INT;
    coho_idx_d_to_s[i] = -1;
    coho_idx_face_d_to_s[i] = -1;
  }

  int ncoho_idx_dst = 0;

  //location of the src face indecies we need to map, assuming face_dst = face_src
  int pos = 0;
  int face_src = face_dst;
  for (int i = 0; i<face_src -1; i++){
    pos += ncoho_idx_src[i];
  }
  int *fcoho_idx_src = coho_idx_src + pos;

  int count = 0;
  for (int i = 0; i < ncoho_idx_src[face_src]; i++) {
    if (fcoho_idx_src[i] == MPP_FILL_INT) break;  // no more data. Rest should be MPP_FILL_INT. TODO bette check?
    int i_src = fcoho_idx_src[i];                 // value of ith index
    //With space index varying fastest - a NetCDF defined 4D index with row_major order  (C-like)
    int ic_src = i_src / (ns_src * nt_src);
    int it_src = (i_src / ns_src) % nt_src;
    int is_src = i_src % ns_src;
    //printf("P: [%d %d] %d %d %d \n",i, i_src, ic_src, it_src, is_src);
    // Now find the corresponing indecies in the destination cohort index:
    int ic_dst = ic_src;  // Verify with SM
    int it_dst = it_src;  // Are we sure tile id maps to same tile id ?
    // idx_map_soil [id_dst] has the source index that correspond to a destiation index.
    // We search the idx_map_soil to find the reverse.
    //is_dst is the location in the array idx_map_soil;
    int is_dst = key_in_array(is_src, idx_map_soil, ns_dst);
    if (is_dst < 0) {
     // printf(" is_dst < 0. Land/sea boundary point being missed ?");
    }else{
      coho_idx_face_d_to_s[count] = face_map_soil[is_src];
      coho_idx_d_to_s[count] = is_dst;// the index position in the src file given
      //The index value in the dst corresponding to src ith index (with value i_src):
      int i_dst = ic_dst * ns_dst * nt_dst + it_dst * ns_dst + is_dst;//row-major order
      coho_idx_dst[count] = i_dst;
      //face_src =;
      count++;
    }
    ///  missing stuff

  }
  ncoho_idx_dst = pos;
  return ncoho_idx_dst;

}

  int compute_dst_coho_idx(int ns_d, int nt_d, int nc_d, int ns_s, int nt_s, int nc_s, int nface_s,
                           int *idx_map_soil, int *face_map_soil,  int* soil_count_cold, int* ncoho_idx_s, int *coho_idx_s,
                           int *coho_idx_d, int * coho_idx_d_to_s,  int *coho_idx_face_d_to_s) {
    int ncoho_idx_d = 0;

    int n_face_coho_idx_d = nc_d * nt_d * nc_d;

    for (int i = 0; i < n_face_coho_idx_d; i++) {
      coho_idx_d[i] = MPP_FILL_INT;
    }

    //The strting point of the source indecies, for each face
    int** f_coho_idx_s = (int**)malloc(nface_s * sizeof(int*));
    f_coho_idx_s[0] = coho_idx_s;
    for (int m = 1; m < nface_s; m++) {
      f_coho_idx_s[m] = f_coho_idx_s[m - 1] + ncoho_idx_s[m - 1];
    }

    //With space index varying fastest - a NetCDF defined 4D index with row_major order  (C-like)
    int count = 0;
    for (int ic_d = 0; ic_d < nc_d; ic_d++) {
      for (int it_d = 0; it_d < nt_d; it_d++) {
        for (int is_d = 0; is_d < ns_d; is_d++) {  // loop over spatial index in dst
          if(soil_count_cold[is_d] > 0){//valid spatial point in cold restart. Verify
            int is_s = idx_map_soil[is_d];   // index in the src
            int if_s = face_map_soil[is_d];  //index of nearest face in the src. TODO:
            int it_s = it_d;  //tile field value source same as dest;
            int ic_s = ic_d;  //cohort field value source same as dst
            //int i_d = nc_d * nt_d * is + nc_d * it + ic; //this may be backwards
            //int i_s = nc_s * nt_s * is_s + nc_s * it + ic;
            //With space index varying fastest - a NetCDF defined 4D index with row_major order  (C-like)
            int i_d = ns_d * nt_d * ic_d + ns_d * it_d + is_d; //this may be backwards
            int i_s = ns_s * nt_s * ic_s + ns_s * it_s + is_s;
            // If i_s is present as a cohort_index value, save the corresponding i_d as next cohort index.

            //todo: Do binary search or hash instead.
            int key_loc = binary_search(i_s, f_coho_idx_s [if_s], ncoho_idx_s[if_s]);
            if (key_loc >= 0) {
              coho_idx_d[count] = i_d; //what is saved when index is saved
              coho_idx_face_d_to_s[count] = if_s;
              coho_idx_d_to_s[count] = key_loc;// the index *position* in the src file; need if_s also for access
              //printf(" i_d i_s face key_loc: %d %d %d %d\n", i_d, i_s, if_s, key_loc);
              count++;
            }else{
              //Most is_s shold not be found
            }
          }
        }
      }
    }
    ncoho_idx_d = count;
    return ncoho_idx_d;

    // Verify: If its possible to have duplicates,  then
    // coho_idx_dst needs to be made a sorted set (sorted + no duplicates)

    if(is_sorted(coho_idx_d, count) && !has_dupes_sorted(coho_idx_d, count)){
      ncoho_idx_d = count;
    }  else{
      printf("*** coho_idx_d not sorted yet **** \n");
      ncoho_idx_d = make_sorted_set(coho_idx_d, count);
      if (ncoho_idx_d != count) {
        printf("Indecies dropped? ncoho_idx_dst=%d count=%d\n", ncoho_idx_d, count);
      }
    }
    return ncoho_idx_d;
  }

/***
 *    *** This function is obsolete. ***
 * Compute the destination cohort index (coho_idx_dst) by looping through all possible values and only
* storing/saving the valid ones. Loop though the indices in the order that they would appear in the
 *  array of the cohort_index. (Otherwise it would be necessary to do a search with dictionary data
 * struct or brute force).
 * Cohort index dimensions are [cohort][tile][lat][lon], and the [tile][lat][lon] dimensions are already
 * mapped in idx_soil_src and face_map_soil. We call [lat][lon] together the space index (is), so we can also represent
 * its dimesnion as [cohort][tile][space]
 * Legend::
 *  src - source
 *  dst - destination
 *  nt -  number f tiles, or tile dimension length
 *  nc -  number of cohorts, or cohort dimension length
 *  ns -  spatial dimension length. @D ns = nx ** ny
*   ic - cohort index
    it - tile index
    icts - cohort_tile_space index
    Note that array coho_idx_src has the indices for all the src faces; while coho_idx_dst
    only
   ***/
int compute_dst_coho_idx_V1(int ns_dst, int nt_dst, int nc_dst, int ns_src, int nt_src, int nc_src, int *idx_map_soil,
                            int *face_map_soil, int face_dst, int *coho_idx_src, int *coho_idx_dst) {
  int ncoho_idx_dst = -1;

  int n_face_coho_idx_dst = nc_dst * nt_dst * nc_dst;
  int n_face_coho_idx_src = ns_src * nt_src * nc_src;

  for (int i = 0; i < n_face_coho_idx_dst; i++) {
    coho_idx_dst[i] = MPP_FILL_INT;
  }

  int *fcoho_idx_src = coho_idx_src + face_dst * n_face_coho_idx_src;
  int pos = 0;

  for (int is = 0; is < ns_dst; is++) {
    printf("is idx_map_soil[is] coho_idx_src[is] : %d  %d  %d\n", is, idx_map_soil[is], coho_idx_src[is]);
  }
  int kstart_hint = 0;
  for (int is = 0; is < ns_dst; is++) {  // loop over spatial index in dst
    int idx = idx_map_soil[is];          // index in the src
    //  int n = face_map_soil[is];  //index of nearest face in the src
    // int is_src = n * nx_src * ny_src + idx;  //spatial index in src; face considered
    int is_src = idx;                             // spatial index in src; face considered
    for (int it = 0; it < nt_dst; it++) {  // tile coord index in src
      //  if( land_idx_map[ntile_dst * is + it] > -1){ // ????
      for (int ic = 0; ic < nc_src; ic++) {  // cohort coord index in src
        int icts = nc_dst * nt_dst * is + nc_dst * it + ic;
        int icts_src = nc_src * nt_src * is_src + nc_src * it + ic;
        // If icst_src is present as a cohort_index value, save the corresponding
        // icst for the destination.
        // TODO: Should check to see if dst_index  cold start cohort_index? Guess not.
        // int next_index = coho_idx_src [pos];
        int key_loc = key_in_array_hints(icts_src, fcoho_idx_src, n_face_coho_idx_src, kstart_hint);
        if (key_loc >= 0) {
          kstart_hint = key_loc;
          coho_idx_dst[pos] = icts;
          printf(" icst icts_src coho_idx .. is is_src : %d  %d  %d %d %d\n", icts, icts_src, coho_idx_src[pos], is,
                 is_src);
          pos++;
        }
      }
    }
  }
  ncoho_idx_dst = pos;
  return ncoho_idx_dst;
}
