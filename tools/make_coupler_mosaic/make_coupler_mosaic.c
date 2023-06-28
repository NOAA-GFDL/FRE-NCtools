/***********************************************************************
 *                   GNU Lesser General Public License
 *
 * This file is part of the GFDL FRE NetCDF tools package (FRE-NCTools).
 *
 * FRE-NCtools is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at
 * your option) any later version.
 *
 * FRE-NCtools is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with FRE-NCTools.  If not, see
 * <http://www.gnu.org/licenses/>.
 **********************************************************************/
/*
 * Modify land grid to match ocean grid at coast and calculate atmos/land,
 * atmos/ocean, and land/ocean overlaps using the Sutherland-Hodgeman polygon
 * clipping algorithm (Sutherland, I. E. and G. W. Hodgeman, 1974:  Reentrant
 * polygon clipping, CACM, 17(1), 32-42).
 *  Warning, when the atmos grid is cubic grid, the number of model points should be
 *  even to avoid tiling error. I will come back to solve this issue in the future.
 *
 */
/* Notes on the Algorithm
	Read in ATM grid tiles xatm,yatm
	Calculate ATM grid cell area area_atm (either GCA or LEGACY)
	Read in LND grid if different from ATM grid
	Calculate LND grid cell area area_lnd, =area_atm if the same grid
	Read in OCN (super) grid xocn,yocn (I added read grid cell area)
	Calculate OCN grid cell area area_ocn
	Extend OCN grid south if it does not cover the southern cap
	*** Extending area_ocn[j=0,:]=area_ocn[j=1,:] is inaccurate if it's used
	Read OCN depth and set ocean mask omask=1 where wet

*/
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "constant.h"
#include "mpp.h"
#include "mpp_io.h"
#include "tool_util.h"
#include "mpp_domain.h"
#include "mosaic_util.h"
#include "create_xgrid.h"
#include "read_mosaic.h"
#define print_grid 0

char *usage[] = {
  "",
  "  make_coupler_mosaic --atmos_mosaic atmos_mosaic.nc --ocean_mosaic ocean_mosaic.nc      ",
  "              --ocean_topog ocean_topog.nc [--land_mosaic land_mosaic.nc] [--wave_mosaic]",
  "              [--sea_level #]  [--interp_order #] [--mosaic_name mosaic_name]           ",
  "              [--area_ratio_thresh #] [--check ] [--verbose] [--print_memory]            ",
  "              [--reproduce_siena]                                                        ",
  " ",
  "make_coupler_mosaic generates three exchange grids for the FMS coupler. The output ",
  "file includes exchange grid files for fluxes between atmosphere and surface (sea ice ",
  "and land), exchange grid files for runoff between land and sea ice. There might be more ",
  "than one exchange grid files between two model solo mosaic because there might be ",
  "multiple tiles in a solo mosaic. The exchange grid information is between model   ",
  "grid, not between supergrid. We assume the refinement ratio between model grid and ",
  "supergrid is 2. Currently we only output the exchange grid on T-cell. ",
  "Besides generating the exchange grid files, make_coupler_mosaic also generates the ",
  "coupler mosaic file (the file name will be mosaic_name.nc) which contains the atmos, ",
  "land and ocean mosaic path, ocean mosaic topog path and exchange grid file path. ",
  "make_coupler_mosaic expects NetCDF format input.",
  " ",
  "make_coupler_mosaic takes the following flags:",
  "",
  "REQUIRED:",
  "",
  "--atmos_mosaic atmos_mosaic.nc specify the atmosphere mosaic information. This file",
  "                               contains list of tile files which specify the grid ",
  "                               information for each tile. Each grid is required to be ",
  "                               logically rectangular grid. The file name can not be 'mosaic.nc' ",
  "",
  "--ocean_mosaic ocean_mosaic.nc specify the ocean mosaic information. This file",
  "                               contains list of tile files which specify the grid ",
  "                               information for each tile. The file name can not be 'mosaic.nc' ",
  " ",
  "--ocean_topog ocean_topog.nc   specify the topography information for ocean mosaic.",
  "                               The field name of the topography is depth_tile# or depth when ",
  "                               ntiles = 1, The topography data is positive down.",
  " ",
  "OPTIONAL FLAGS",
  "",
  "--land_mosaic land_mosaic.nc   specify the land mosaic information. This file",
  "                               contains list of tile files which specify the grid ",
  "                               information for each tile. Each grid is required to be ",
  "                               logically rectangular grid. When land_mosaic is not specified,",
  "                               atmosphere mosaic will be used to specify land mosaic.",
  "                               The file name can not be 'mosaic.nc'.",
  " ",
  "--wave_mosaic wave_mosaic.nc   specify the wave mosaic information. This file",
  "                               contains list of tile files which specify the grid ",
  "                               information for each tile. Each grid is required to be ",
  "                               logically rectangular grid. When wave_mosaic is specified, ",
  "                               exchange grid information between wave mosaic and ocean ",
  "                               mosaic will be generated. Otherwise none of both will be created.",
  " ",
  "--interp_order #               specify the order of conservative interplation. Its value ",
  "                               can be 1 ( linear order ) or 2 ( second order ) with default ",
  "                               value 2.                                                     ",
  "                                                                                            ",
  "--sea_level #                  specify the sea level ( in meters ) and its value will be used",
  "                               to determine land/sea mask. When topography of  ",
  "                               a grid cell is less than sea level, this grid cell will be land,",
  "                               otherwise it will be ocean. Default value is 0",
  " ",
  "--mosaic_name mosaic_name      coupler mosaic name. The output coupler mosaic file will be ",
  "                               mosaic_name.nc. default value is 'mosaic'. ",
  "                                                                                        ",
  "--area_ratio_thresh #          Criteria to decide if an overlap between any two model is an ",
  "                               exchange grid or not. When overlap area/model grid area is greater ",
  "                               than area_ratio_thresh, that overlap is saved as exchange grid. ",
  "                               The default value is 1.e-6                                      ",
  "                                                                                              ",
  "--check                        check the tiling error",
  " ",
  "--print_memory                 debug memory usage when it is set                       ",
  " ",
  "--reproduce_siena              Set to reproduce siena shared codes results              ",
  "    ",
  "--rotate_poly                  Set to calculate polar polygon areas by caculating the area of a copy ",
  "                               of the polygon, with the copy being rotated far away from the pole. ",
  " ",
  "--verbose                      Set --verbose to print out messages during running.        ",
  "",

  "A sample call to make_coupler_mosaic that makes exchange grids for atmosphere, land and ocean ",
  "mosaic (atmosphere and land are coincident) is: ",
  "",
  "  make_coupler_mosaic --atmos_mosaic atmos_mosaic.nc --ocean_mosaic ocean_mosaic.nc ",
  "                  --ocean_topog ocean_topog.nc",
  "",
  NULL };

#define MAXXGRIDFILE 100
#define MX 2000
#define TINY_VALUE (1.e-7)
#define TOLORENCE (1.e-4)
#define MIN_AREA_FRAC (1.e-4)
#define GREAT_CIRCLE_CLIP 1
#define LEGACY_CLIP 2


char grid_version[] = "0.2";
char tagname[] = "$Name: fre-nctools-bronx-10 $";

/* This file will get the directory that stores the file and the file name (without the dir path) */
void get_file_dir_and_name(char *file, char *filedir, char *filename)
{
  char *fptr=NULL;
  int siz;

  fptr = strrchr(file, '/');

    if(!fptr) {
      strcpy(filename, file);
      strcpy(filedir, "./");
    }
    else {
      ++fptr;
      siz = fptr - file;
      strcpy(filename, fptr);
      strncpy(filedir, file, siz);
    }
};


void get_global_grid(const char *grid_file, int nx, int ny, int x_refine, int y_refine, double *x, double *y)
{
  double *x_local, *y_local, *tmp;
  size_t start[4], nread[4];
  int nxc, nyc, isc, iec, jsc, jec;
  int isc2, iec2, jsc2, jec2;
  int nxc2, nyc2, layout[2];
  int fid, vid, i, j;
  domain2D Dom;

  /* define a temporary domain to read data on local and then
     use mpp_global_field to get global field.
     This may solve the IO issue for high resolution grid */

  mpp_define_layout(nx, ny, mpp_npes(), layout);
  mpp_define_domain2d(nx, ny, layout, 0, 0, &Dom);
  mpp_get_compute_domain2d(Dom, &isc, &iec, &jsc, &jec);
  nxc = iec - isc + 1;
  nyc = jec - jsc + 1;
  x_local = (double *)malloc((nxc+1)*(nyc+1)*sizeof(double));
  y_local = (double *)malloc((nxc+1)*(nyc+1)*sizeof(double));
  isc2 = x_refine*isc; iec2 = x_refine*iec + 1;
  jsc2 = y_refine*jsc; jec2 = y_refine*jec + 1;
  nxc2 = iec2 - isc2 + 1;
  nyc2 = jec2 - jsc2 + 1;
  tmp        = (double *)malloc((nxc2+1)*(nyc2+1)*sizeof(double));

  start[0] = jsc2;   start[1] = isc2;   start[2] = 0; start[3] = 0;
  nread[0] = nyc2+1; nread[1] = nxc2+1; nread[2] = 1; nread[3] = 1;
  fid = mpp_open(grid_file, MPP_READ);
  vid = mpp_get_varid(fid, "x");
  mpp_get_var_value_block(fid, vid, start, nread, tmp);
  for(j=0; j<=nyc; j++) for(i=0; i<=nxc; i++)
    x_local[j*(nxc+1)+i] = tmp[(j*y_refine)*(nxc2+1)+i*x_refine];
  mpp_global_field_all_double(Dom, nxc+1, nyc+1, x_local, x);

  vid = mpp_get_varid(fid, "y");
  mpp_get_var_value_block(fid, vid, start, nread, tmp);
  for(j=0; j<=nyc; j++) for(i=0; i<=nxc; i++)
    y_local[j*(nxc+1)+i] = tmp[(j*y_refine)*(nxc2+1)+i*x_refine];
  mpp_global_field_all_double(Dom, nxc+1, nyc+1, y_local, y);
  mpp_delete_domain2d(&Dom);
  mpp_close(fid);

  free(tmp);
  free(x_local);
  free(y_local);
}

void get_global_data(const char *data_file, const char *fieldname, int nx, int ny, double *data,
		     int ioff, int joff)
{
  double *data_local;
  size_t start[4], nread[4];
  int nxc, nyc, isc, iec, jsc, jec;
  int fid, vid, layout[2];
  domain2D Dom;

  mpp_define_layout(nx, ny, mpp_npes(), layout);
  mpp_define_domain2d(nx, ny, layout, 0, 0, &Dom);
  mpp_get_compute_domain2d(Dom, &isc, &iec, &jsc, &jec);
  nxc = iec - isc + 1;
  nyc = jec - jsc + 1;
  data_local    = (double *)malloc((nxc+ioff)*(nyc+joff)*sizeof(double));
  start[0] = jsc;      start[1] = isc;      start[2] = 0; start[3] = 0;
  nread[0] = nyc+joff; nread[1] = nxc+ioff; nread[2] = 1; nread[3] = 1;

  fid = mpp_open(data_file, MPP_READ);
  vid = mpp_get_varid(fid, fieldname);
  mpp_get_var_value_block(fid, vid, start, nread, data_local);
  mpp_close(fid);

  mpp_global_field_all_double(Dom, nxc+ioff, nyc+joff, data_local, data);
  free(data_local);

}

void get_grid_global_area(int nx, int ny, const double *x, const double *y, double *area)
{
  double *x_local, *y_local, *area_local;
  int nxc, nyc, isc, iec, jsc, jec;
  int isc2, iec2, jsc2, jec2;
  int nxc2, nyc2, layout[2];
  int i, j;
  domain2D Dom;

  mpp_define_layout(nx, ny, mpp_npes(), layout);
  mpp_define_domain2d(nx, ny, layout, 0, 0, &Dom);
  mpp_get_compute_domain2d(Dom, &isc, &iec, &jsc, &jec);
  nxc = iec - isc + 1;
  nyc = jec - jsc + 1;
  x_local    = (double *)malloc((nxc+1)*(nyc+1)*sizeof(double));
  y_local    = (double *)malloc((nxc+1)*(nyc+1)*sizeof(double));

  for(j=0; j<=nyc; j++) for(i=0; i<=nxc; i++) {
    x_local[j*(nxc+1)+i] = x[(j+jsc)*(nx+1)+i+isc];
    y_local[j*(nxc+1)+i] = y[(j+jsc)*(nx+1)+i+isc];
  }

  area_local = (double *)malloc(nxc*nyc*sizeof(double));
  get_grid_area(&nxc, &nyc, x_local, y_local, area_local);
  mpp_global_field_all_double(Dom, nxc, nyc, area_local, area);
  free(x_local);
  free(y_local);
  free(area_local);
  mpp_delete_domain2d(&Dom);


}

int get_nest_contact(const int *nx, const int *ny, int ncontacts, const int *tile1, const int *tile2,
		       const int *istart1, const int *iend1, const int *jstart1, const int *jend1,
		       const int *istart2, const int *iend2, const int *jstart2, const int *jend2,
		       int *tile_nest, int *tile_parent, int *is_nest, int *ie_nest, int *js_nest,
		       int *je_nest, int *is_parent, int *ie_parent, int *js_parent, int *je_parent)
{

  int n, nx1_contact, ny1_contact, nx2_contact, ny2_contact;
  int t1, t2;
  int nnest=0;

  for(n=0; n<ncontacts; n++) {
    /* same tile could not be nested */
    if(tile1[n] == tile2[n]) continue;
    nx1_contact = iend1[n]-istart1[n]+1;
    ny1_contact = jend1[n]-jstart1[n]+1;
    nx2_contact = iend2[n]-istart2[n]+1;
    ny2_contact = jend2[n]-jstart2[n]+1;
    t1 = tile1[n] - 1;
    t2 = tile2[n] - 1;
    /* For nesting, the contact index of one tile must match its global domain */
    if( (nx[t1] != nx1_contact || ny[t1] != ny1_contact ) &&
	(nx[t2] != nx2_contact || ny[t2] != ny2_contact ) ) continue;
    if(nx1_contact == nx2_contact && ny1_contact == ny2_contact)
      mpp_error("make_coupler_mosaic(get_nest_contact): There is no refinement for the overlapping region");
    nnest ++;
    if(nnest>1)mpp_error("make_coupler_mosaic(get_nest_contact): only support one nest region, contact developer");
    if(nx2_contact*ny2_contact > nx1_contact*ny1_contact) {
      if(nx2_contact%nx1_contact || ny2_contact%ny1_contact )
	if(mpp_pe()==mpp_root_pe()) mpp_error("make_coupler_mosaic(get_nest_contact):it is not a integer refinement");
      is_nest    [0] = istart2[n];
      ie_nest    [0] = iend2  [n];
      js_nest    [0] = jstart2[n];
      je_nest    [0] = jend2  [n];
      tile_nest  [0] = tile2  [n]-1;
      is_parent  [0] = istart1[n];
      ie_parent  [0] = iend1  [n];
      js_parent  [0] = jstart1[n];
      je_parent  [0] = jend1  [n];
      tile_parent[0] = tile1  [n]-1;
    }
    else {
      if(nx1_contact%nx2_contact || ny1_contact%ny2_contact )
	mpp_error("make_coupler_mosaic(get_nest_contact):it is not a integer refinement");
      is_nest    [0] = istart1[n];
      ie_nest    [0] = iend1  [n];
      js_nest    [0] = jstart1[n];
      je_nest    [0] = jend1  [n];
      tile_nest  [0] = tile1  [n]-1;
      is_parent  [0] = istart2[n];
      ie_parent  [0] = iend2  [n];
      js_parent  [0] = jstart2[n];
      je_parent  [0] = jend2  [n];
      tile_parent[0] = tile2  [n]-1;
    }
  }

  return nnest;

}

int main (int argc, char *argv[])
{
  int c, i, same_mosaic;
  extern char *optarg;
  char *omosaic  = NULL;
  char *amosaic  = NULL;
  char *lmosaic  = NULL;
  char *wmosaic  = NULL;
  char *otopog   = NULL;
  char mosaic_name[STRING] = "mosaic", mosaic_file[STRING];
  char omosaic_name[STRING], amosaic_name[STRING], lmosaic_name[STRING], wmosaic_name[STRING];
  char **otile_name=NULL, **atile_name=NULL, **ltile_name=NULL, **wtile_name=NULL;
  int x_refine = 2, y_refine = 2;
  int  interp_order = 2;
  double area_ratio_thresh = 1.0e-6;
  unsigned int check = 0;
  unsigned int verbose = 0;
  int errflg = (argc == 1);
  char history[512];
  int nfile_lxo=0, nfile_axo=0, nfile_axl=0;
  int nfile_wxo=0;
  char lxo_file[MAXXGRIDFILE][STRING];
  char axo_file[MAXXGRIDFILE][STRING];
  char axl_file[MAXXGRIDFILE][STRING];
  char wxo_file[MAXXGRIDFILE][STRING];
  char amosaic_dir[STRING], amosaic_file[STRING];
  char lmosaic_dir[STRING], lmosaic_file[STRING];
  char omosaic_dir[STRING], omosaic_file[STRING];
  char wmosaic_dir[STRING], wmosaic_file[STRING];
  char otopog_dir[STRING], otopog_file[STRING];
  int    ntile_ocn, ntile_atm, ntile_lnd, ntile_wav;
  int    *nxo = NULL, *nyo = NULL, *nxa = NULL, *nya = NULL, *nxl = NULL, *nyl = NULL;
  int    *nxw = NULL, *nyw = NULL;
  double **xocn = NULL, **yocn = NULL, **xatm = NULL, **yatm = NULL, **xlnd = NULL, **ylnd = NULL;
  double **xwav = NULL, **ywav = NULL;
  double **cart_xatm=NULL, **cart_yatm=NULL, **cart_zatm=NULL;
  double **cart_xocn=NULL, **cart_yocn=NULL, **cart_zocn=NULL;
  double **cart_xlnd=NULL, **cart_ylnd=NULL, **cart_zlnd=NULL;
  double **cart_xwav=NULL, **cart_ywav=NULL, **cart_zwav=NULL;
  double **area_ocn = NULL, **area_lnd = NULL, **area_atm = NULL, **area_wav = NULL;
  double **atm_xarea=NULL;
  double **omask = NULL;
  double sea_level = 0.;
  int    clip_method = LEGACY_CLIP;
  int    atm_great_circle_algorithm=0;
  int    lnd_great_circle_algorithm=0;
  int    ocn_great_circle_algorithm=0;

  int    lnd_same_as_atm = 0;
  int    ocn_same_as_atm = 0;
  int    wav_same_as_ocn = 0;
  int    option_index = 0;
  double axo_area_sum = 0, axl_area_sum = 0;
  double axo_area_sum_nest = 0, axl_area_sum_nest = 0;
  double wxo_area_sum = 0;
  int    ocn_south_ext = 0;
  int    tile_nest, is_nest, ie_nest, js_nest, je_nest;
  int    tile_parent, is_parent, ie_parent, js_parent, je_parent;
  int    print_memory=0;
  int    reproduce_siena=0;
  int    rotate_poly=0;

  static struct option long_options[] = {
    {"atmos_mosaic",         required_argument, NULL, 'a'},
    {"land_mosaic",          required_argument, NULL, 'l'},
    {"ocean_mosaic",         required_argument, NULL, 'o'},
    {"wave_mosaic",          required_argument, NULL, 'w'},
    {"ocean_topog",          required_argument, NULL, 't'},
    {"sea_level",            required_argument, NULL, 's'},
    {"interp_order",         required_argument, NULL, 'i'},
    {"mosaic_name",          required_argument, NULL, 'm'},
    {"area_ratio_thresh",    required_argument, NULL, 'r'},
    {"check",                no_argument,       NULL, 'n'},
    {"verbose",              no_argument,       NULL, 'v'},
    {"print_memory",         no_argument,       NULL, 'p'},
    {"reproduce_siena",      no_argument,       NULL, 'q'},
    {"rotate_poly",          no_argument,       NULL, 'u'},
    {NULL, 0, NULL, 0}
  };

  mpp_init(&argc, &argv);
  mpp_domain_init();

  /*
   * process command line
   */

  while ((c = getopt_long(argc, argv, "i:", long_options, &option_index) ) != -1)
    switch (c) {
    case 'a':
      amosaic = optarg;
      break;
    case 'l':
      lmosaic = optarg;
      break;
    case 'o':
      omosaic = optarg;
      break;
    case 'w':
      wmosaic = optarg;
      break;
    case 't':
      otopog = optarg;
      break;
    case 'i':
      interp_order = atoi(optarg);
      break;
    case 's':
      sea_level = atof(optarg);
      break;
    case 'm':
      strcpy(mosaic_name,optarg);
      break;
    case 'r':
      area_ratio_thresh = atof(optarg);
    case 'n':
      check = 1;
      break;
    case 'v':
      verbose = 1;
      break;
    case 'p':
      print_memory = 1;
      break;
    case 'q':
      reproduce_siena = 1;
      break;
    case 'u':
      rotate_poly = 1;
      break;
    case '?':
      errflg++;
    }
  if (errflg || !amosaic || !omosaic || !otopog) {
    char **u = usage;
    while (*u) { fprintf(stderr, "%s\n", *u); u++; }
    exit(2);
  }

  /* interp_order should be 1 or 2 */
  if(interp_order != 1 && interp_order !=2 )mpp_error("make_coupler_mosaic: interp_order should be 1 or 2");

  strcpy(history,argv[0]);

  for(i=1;i<argc;i++) {
    strcat(history, " ");
    strcat(history, argv[i]);
  }

  /*if lmosaic is not specifiied, assign amosaic value to it */
  if(!lmosaic) lmosaic = amosaic;

  if(reproduce_siena) set_reproduce_siena_true();

  if(rotate_poly) set_rotate_poly_true();

  /*mosaic_file can not have the same name as amosaic, lmosaic or omosaic, also the file name of
    amosaic, lmosaic, omosaic can not be "mosaic.nc"
  */
  sprintf(mosaic_file, "%s.nc", mosaic_name);
  get_file_dir_and_name(amosaic, amosaic_dir, amosaic_file);
  get_file_dir_and_name(lmosaic, lmosaic_dir, lmosaic_file);
  get_file_dir_and_name(omosaic, omosaic_dir, omosaic_file);
  if(wmosaic) get_file_dir_and_name(wmosaic, wmosaic_dir, wmosaic_file);
  get_file_dir_and_name(otopog, otopog_dir, otopog_file);
  if( !strcmp(mosaic_file, amosaic_file) || !strcmp(mosaic_file, lmosaic_file) || !strcmp(mosaic_file, omosaic_file) )
    mpp_error("make_coupler_mosaic: mosaic_file can not have the same name as amosaic, lmosaic or omosaic");
  if( !strcmp(amosaic_file, "mosaic.nc") || !strcmp(lmosaic_file, "mosaic.nc") || !strcmp(omosaic_file, "mosaic.nc") )
    mpp_error("make_coupler_mosaic: the file name of amosaic, lmosaic or omosaic can not be mosaic.nc");

  if(print_memory) print_mem_usage("before read atmosphere grid");

  /*
   * Read atmosphere grid
   */
  {
    int n, m_fid, g_fid, vid, gid, tid;
    int ncontacts, nnest;
    int *tile1=NULL, *istart1=NULL, *iend1=NULL, *jstart1=NULL, *jend1=NULL;
    int *tile2=NULL, *istart2=NULL, *iend2=NULL, *jstart2=NULL, *jend2=NULL;
    size_t start[4], nread[4];
    char dir[STRING], filename[STRING], file[2*STRING];

    for(n=0; n<4; n++) {
      start[n] = 0;
      nread[n] = 1;
    }

    m_fid = mpp_open(amosaic, MPP_READ);
    vid = mpp_get_varid(m_fid, "mosaic");
    mpp_get_var_value(m_fid, vid, amosaic_name);
    ntile_atm  = mpp_get_dimlen(m_fid, "ntiles");
    nxa        = (int *) malloc (ntile_atm*sizeof(int));
    nya        = (int *) malloc (ntile_atm*sizeof(int));
    xatm       = (double **) malloc( ntile_atm*sizeof(double *));
    yatm       = (double **) malloc( ntile_atm*sizeof(double *));
    area_atm   = (double **) malloc( ntile_atm*sizeof(double *));
    if(check) atm_xarea  = (double **) malloc( ntile_atm*sizeof(double *));
    atile_name = (char **)malloc(ntile_atm*sizeof(char *));
    /* grid should be located in the same directory of mosaic file */
    get_file_path(amosaic, dir);
    gid = mpp_get_varid(m_fid, "gridfiles");
    tid = mpp_get_varid(m_fid, "gridtiles");
    for(n=0; n<ntile_atm; n++) {
      int i, j;

      start[0] = n; start[1] = 0; nread[0] = 1; nread[1] = STRING;
      mpp_get_var_value_block(m_fid, gid, start, nread, filename);
      atile_name[n] = (char *)malloc(STRING*sizeof(char));
      mpp_get_var_value_block(m_fid, tid, start, nread,  atile_name[n]);
      sprintf(file, "%s/%s", dir, filename);
      g_fid = mpp_open(file, MPP_READ);
      nxa[n] = mpp_get_dimlen(g_fid, "nx");
      nya[n] = mpp_get_dimlen(g_fid, "ny");
      /* check if use great_circle_algorithm */
      {
	int great_circle_algorithm=0;
	great_circle_algorithm = get_great_circle_algorithm(g_fid);
	if(n>0) {
	  if( atm_great_circle_algorithm != great_circle_algorithm)
	    mpp_error("make_topog: atribute 'great_circle_algorithm' of field 'tile' have different value for different tile");
	}
	atm_great_circle_algorithm = great_circle_algorithm;
      }

      mpp_close(g_fid);
      if(nxa[n]%x_refine != 0 ) mpp_error("make_coupler_mosaic: atmos supergrid x-size can not be divided by x_refine");
      if(nya[n]%y_refine != 0 ) mpp_error("make_coupler_mosaic: atmos supergrid y-size can not be divided by y_refine");
      nxa[n] /= x_refine;
      nya[n] /= y_refine;
      xatm[n]     = (double *)malloc((nxa[n]+1)*(nya[n]+1)*sizeof(double));
      yatm[n]     = (double *)malloc((nxa[n]+1)*(nya[n]+1)*sizeof(double));
      area_atm[n] = (double *)malloc((nxa[n]  )*(nya[n]  )*sizeof(double));
      if(check) {
	atm_xarea[n]= (double *)malloc((nxa[n]  )*(nya[n]  )*sizeof(double));
        for(i=0; i<nxa[n]*nya[n]; i++) atm_xarea[n][i] = 0;
      }
      get_global_grid(file, nxa[n], nya[n], x_refine, y_refine, xatm[n], yatm[n]);
      /* convert to radians */
      for(i=0; i<(nxa[n]+1)*(nya[n]+1); i++) {
	xatm[n][i] *= D2R;
	yatm[n][i] *= D2R;
      }

    }

    if(atm_great_circle_algorithm)
      clip_method = GREAT_CIRCLE_CLIP;
    else
      clip_method = LEGACY_CLIP;

    /* Currenly only implement interp_order = 1 when clip_method = "conserve_great_circle" */
    if( interp_order != 1 && clip_method == GREAT_CIRCLE_CLIP)
      mpp_error("make_coupler_mosaic:  Currenly only implement interp_order = 1 when clip_method = 'conserve_great_circle', contact developer");


    /* compute atm_area */
    if(clip_method == GREAT_CIRCLE_CLIP) {
      cart_xatm       = (double **) malloc( ntile_atm*sizeof(double *));
      cart_yatm       = (double **) malloc( ntile_atm*sizeof(double *));
      cart_zatm =  (double **) malloc( ntile_atm*sizeof(double *));
      for(n=0; n<ntile_atm; n++) {
	cart_xatm[n]     = (double *)malloc((nxa[n]+1)*(nya[n]+1)*sizeof(double));
	cart_yatm[n]     = (double *)malloc((nxa[n]+1)*(nya[n]+1)*sizeof(double));
	cart_zatm[n]     = (double *)malloc((nxa[n]+1)*(nya[n]+1)*sizeof(double));
	latlon2xyz((nxa[n]+1)*(nya[n]+1), xatm[n], yatm[n], cart_xatm[n], cart_yatm[n], cart_zatm[n]);
	get_grid_great_circle_area(&(nxa[n]), &(nya[n]), xatm[n], yatm[n], area_atm[n]);
      }
    }
    else {
      for(n=0; n<ntile_atm; n++) {
	int i, j;
        //Calculate the ATM grid cell area (not read from grid files)
	get_grid_global_area(nxa[n], nya[n], xatm[n], yatm[n], area_atm[n]);
	for(j=0; j<nya[n]; j++) for(i=0; i<nxa[n]; i++) {
	  if(area_atm[n][j*nxa[n]+i] <= 0 ||area_atm[n][j*nxa[n]+i] >1.0e14  ) printf("Possible error in the calculated area_atm in tile %d, i=%d, j=%d, area=%f\n",n+1,i,j,area_atm[n][j*nxa[n]+i]);
	}
      }
    }

    mpp_close(m_fid);
    /* read the contact information in atmos_mosaic to see if there is a nested grid */
    ncontacts = read_mosaic_ncontacts( amosaic );
    is_nest = -1; ie_nest = -1; js_nest = -1; je_nest = -1;
    js_parent = -1; ie_parent = -1; js_parent = -1; je_parent = -1;
    tile_nest = -1; tile_parent = -1;
    if(ncontacts >0) {
      tile1   = (int *)malloc(ncontacts*sizeof(int));
      tile2   = (int *)malloc(ncontacts*sizeof(int));
      istart1 = (int *)malloc(ncontacts*sizeof(int));
      iend1   = (int *)malloc(ncontacts*sizeof(int));
      jstart1 = (int *)malloc(ncontacts*sizeof(int));
      jend1   = (int *)malloc(ncontacts*sizeof(int));
      istart2 = (int *)malloc(ncontacts*sizeof(int));
      iend2   = (int *)malloc(ncontacts*sizeof(int));
      jstart2 = (int *)malloc(ncontacts*sizeof(int));
      jend2   = (int *)malloc(ncontacts*sizeof(int));
      read_mosaic_contact(amosaic, tile1, tile2, istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2);
      nnest = get_nest_contact( nxa, nya, ncontacts, tile1, tile2, istart1, iend1, jstart1, jend1,
				istart2, iend2, jstart2, jend2, &tile_nest, &tile_parent, &is_nest, &ie_nest,
				&js_nest, &je_nest, &is_parent, &ie_parent, &js_parent, &je_parent);
      free(tile1);
      free(tile2);
      free(istart1);
      free(iend1);
      free(jstart1);
      free(jend1);
      free(istart2);
      free(iend2);
      free(jstart2);
      free(jend2);

    }
  }
  if(print_memory)print_mem_usage("after read atmosphere grid");
  /*
   * Read land grid
   */
  if (strcmp(lmosaic, amosaic) ) { /* land mosaic is different from atmosphere mosaic */
    int n, m_fid, g_fid, vid, gid, tid;
    size_t start[4], nread[4];
    char dir[STRING], filename[STRING], file[2*STRING];

    for(n=0; n<4; n++) {
      start[n] = 0;
      nread[n] = 1;
    }
    m_fid = mpp_open(lmosaic, MPP_READ);
    vid = mpp_get_varid(m_fid, "mosaic");
    mpp_get_var_value(m_fid, vid, lmosaic_name);
    ntile_lnd  = mpp_get_dimlen(m_fid, "ntiles");
    nxl        = (int *) malloc (ntile_lnd*sizeof(int));
    nyl        = (int *) malloc (ntile_lnd*sizeof(int));
    xlnd       = (double **) malloc( ntile_lnd*sizeof(double *));
    ylnd       = (double **) malloc( ntile_lnd*sizeof(double *));
    area_lnd   = (double **) malloc( ntile_lnd*sizeof(double *));
    ltile_name = (char **)malloc(ntile_lnd*sizeof(char *));
    /* grid should be located in the same directory of mosaic file */
    get_file_path(lmosaic, dir);
    gid = mpp_get_varid(m_fid, "gridfiles");
    tid = mpp_get_varid(m_fid, "gridtiles");
    for(n=0; n<ntile_lnd; n++) {
      int i, j;

      start[0] = n; start[1] = 0; nread[0] = 1; nread[1] = STRING;
      mpp_get_var_value_block(m_fid, gid, start, nread, filename);
      ltile_name[n] = (char *)malloc(STRING*sizeof(char));
      mpp_get_var_value_block(m_fid, tid, start, nread,  ltile_name[n]);
      sprintf(file, "%s/%s", dir, filename);
      g_fid = mpp_open(file, MPP_READ);
      nxl[n] = mpp_get_dimlen(g_fid, "nx");
      nyl[n] = mpp_get_dimlen(g_fid, "ny");
      /* check if use great_circle_algorithm */
      {
	int great_circle_algorithm=0;
	great_circle_algorithm = get_great_circle_algorithm(g_fid);
	if(n>0) {
	  if( lnd_great_circle_algorithm != great_circle_algorithm)
	    mpp_error("make_topog: atribute 'great_circle_algorithm' of field 'tile' have different value for different tile");
	}
	lnd_great_circle_algorithm = great_circle_algorithm;
      }

      mpp_close(g_fid);
      if(nxl[n]%x_refine != 0 ) mpp_error("make_coupler_mosaic: land supergrid x-size can not be divided by x_refine");
      if(nyl[n]%y_refine != 0 ) mpp_error("make_coupler_mosaic: land supergrid y-size can not be divided by y_refine");
      nxl[n]      /= x_refine;
      nyl[n]      /= y_refine;
      xlnd[n]     = (double *)malloc((nxl[n]+1)*(nyl[n]+1)*sizeof(double));
      ylnd[n]     = (double *)malloc((nxl[n]+1)*(nyl[n]+1)*sizeof(double));
      area_lnd[n] = (double *)malloc((nxl[n]  )*(nyl[n]  )*sizeof(double));
      get_global_grid(file, nxl[n], nyl[n], x_refine, y_refine, xlnd[n], ylnd[n]);
      /*scale grid from degree to radian, because create_xgrid assume the grid is in radians */
      for(i=0; i<(nxl[n]+1)*(nyl[n]+1); i++) {
	xlnd[n][i] *= D2R;
	ylnd[n][i] *= D2R;
      }

    }

    /* compute lnd_area */
    if(clip_method == GREAT_CIRCLE_CLIP) {
      cart_xlnd = (double **) malloc( ntile_lnd*sizeof(double *));
      cart_ylnd = (double **) malloc( ntile_lnd*sizeof(double *));
      cart_zlnd =  (double **) malloc( ntile_lnd*sizeof(double *));
      for(n=0; n<ntile_lnd; n++) {
	cart_xlnd[n]     = (double *)malloc((nxl[n]+1)*(nyl[n]+1)*sizeof(double));
	cart_ylnd[n]     = (double *)malloc((nxl[n]+1)*(nyl[n]+1)*sizeof(double));
	cart_zlnd[n]     = (double *)malloc((nxl[n]+1)*(nyl[n]+1)*sizeof(double));
	latlon2xyz((nxl[n]+1)*(nyl[n]+1), xlnd[n], ylnd[n], cart_xlnd[n], cart_ylnd[n], cart_zlnd[n]);
	get_grid_great_circle_area(&(nxl[n]), &(nyl[n]), xlnd[n], ylnd[n], area_lnd[n]);
      }
    }
    else {
      for(n=0; n<ntile_lnd; n++) {
	get_grid_global_area(nxl[n], nyl[n], xlnd[n], ylnd[n], area_lnd[n]);
      }
    }
    mpp_close(m_fid);
  }
  else { /* land mosaic is same as atmosphere mosaic */
    if(tile_nest>=0)mpp_error("make_coupler_mosaic: land mosaic must be different from atmos mosaic "
			       "when there is nest region in atmosphere mosaic");
    ntile_lnd = ntile_atm;
    nxl = nxa;
    nyl = nya;
    if(clip_method == GREAT_CIRCLE_CLIP) {
      cart_xlnd = cart_xatm;
      cart_ylnd = cart_yatm;
      cart_zlnd = cart_zatm;
    }
    else {
      xlnd = xatm;
      ylnd = yatm;
    }
    area_lnd = area_atm;
    lnd_same_as_atm = 1;
    strcpy(lmosaic_name, amosaic_name);
    ltile_name = atile_name;
    lnd_great_circle_algorithm = atm_great_circle_algorithm;
  }
  int n;
  for(n=0; n<ntile_lnd; n++) {
    if(mpp_pe()==mpp_root_pe() && verbose) printf("Number of ATM and LND cells for tile %d are %d, %d.\n",n+1,nxa[n]*nya[n],nxl[n]*nyl[n]);
    if(nxa[n]*nya[n] != nxl[n]*nyl[n]) printf("Warning: Number of ATM and LND cells for tile %d are not equal %d, %d.\n",n+1,nxa[n]*nya[n],nxl[n]*nyl[n]);
  }
  if(print_memory)print_mem_usage("after read land grid");

  if (strcmp(omosaic, amosaic) == 0 ) ocn_same_as_atm = 1;
  /*
   * Read ocean grid boundaries and mask (where water is) for each tile within the mosaic.
   */
  {
    int n, ntiles, m_fid, g_fid, t_fid, vid, gid, tid;
    size_t start[4], nread[4];
    char dir[STRING], filename[STRING], file[2*STRING];

    for(n=0; n<4; n++) {
      start[n] = 0;
      nread[n] = 1;
    }
    m_fid = mpp_open(omosaic, MPP_READ);
    vid = mpp_get_varid(m_fid, "mosaic");
    mpp_get_var_value(m_fid, vid, omosaic_name);
    ntile_ocn  = mpp_get_dimlen(m_fid, "ntiles");
    nxo      = (int     *) malloc(ntile_ocn*sizeof(int));
    nyo      = (int     *) malloc(ntile_ocn*sizeof(int));
    xocn     = (double **) malloc(ntile_ocn*sizeof(double *));
    yocn     = (double **) malloc(ntile_ocn*sizeof(double *));
    area_ocn = (double **) malloc(ntile_ocn*sizeof(double *));
    otile_name = (char **) malloc(ntile_ocn*sizeof(char *));
    /* grid should be located in the same directory of mosaic file */
    get_file_path(omosaic, dir);
    gid = mpp_get_varid(m_fid, "gridfiles");
    tid = mpp_get_varid(m_fid, "gridtiles");

    /* For the purpose of reproducing between processor count, the layout
       is set to (1, npes). */

    for(n=0; n<ntile_ocn; n++) {
      double *tmpx, *tmpy;
      int i, j;
      double min_atm_lat, min_lat;
      int nyo_old;

      start[0] = n; start[1] = 0; nread[0] = 1; nread[1] = STRING;
      mpp_get_var_value_block(m_fid, gid, start, nread, filename);
      otile_name[n] = (char *)malloc(STRING*sizeof(char));
      mpp_get_var_value_block(m_fid, tid, start, nread,  otile_name[n]);
      sprintf(file, "%s/%s", dir, filename);
      g_fid = mpp_open(file, MPP_READ);
      nxo[n] = mpp_get_dimlen(g_fid, "nx");
      nyo[n] = mpp_get_dimlen(g_fid, "ny");
      /* check if use great_circle_algorithm */
      {
	int great_circle_algorithm=0;
	great_circle_algorithm = get_great_circle_algorithm(g_fid);
	if(n>0) {
	  if( ocn_great_circle_algorithm != great_circle_algorithm)
	    mpp_error("make_topog: atribute 'great_circle_algorithm' of field 'tile' have different value for different tile");
	}
	ocn_great_circle_algorithm = great_circle_algorithm;
      }

      mpp_close(g_fid);

      if(nxo[n]%x_refine != 0 ) mpp_error("make_coupler_mosaic: ocean supergrid x-size can not be divided by x_refine");
      if(nyo[n]%y_refine != 0 ) mpp_error("make_coupler_mosaic: ocean supergrid y-size can not be divided by y_refine");
      nxo[n] /= x_refine;
      nyo[n] /= y_refine;

      tmpx    = (double *)malloc((nxo[n]+1)*(nyo[n]+1)*sizeof(double));
      tmpy    = (double *)malloc((nxo[n]+1)*(nyo[n]+1)*sizeof(double));
      get_global_grid(file, nxo[n], nyo[n], x_refine, y_refine, tmpx, tmpy);

      /* sometimes the ocean grid only covers part of the globe, especially may not cover
	 the south pole region. In order to get all the exchange grid between atmosXland,
	 we need to extend one point to cover the whole atmosphere. This needs the
	 assumption of one-tile ocean. Also we assume the latitude is the along j=0
      */
      if(ntile_ocn == 1) {
	int na;
	int is_uniform;

	/* check if the latitude is uniform or not at j=1 */
        is_uniform = 1;
	for(i=1; i<=nxo[n]; i++) {
          if(tmpy[i] != tmpy[i-1]) {
	    is_uniform = 0;
	  }
	}
	if(!is_uniform && mpp_pe()==mpp_root_pe()) printf("\nNOTE from make_coupler_mosaic: ocean grid latitude is not uniform along j = 1\n");

	/* calculate the minimum of latitude of atmosphere grid */
	min_atm_lat = -90.*D2R;
	if(tmpy[0]*D2R > min_atm_lat + TINY_VALUE) { /* extend one point in south direction*/
	  ocn_south_ext = 1;
	  if(mpp_pe()==mpp_root_pe())printf("make_coupler_mosaic: one row is added to the south end to cover the globe\n");
	}
      }
      nyo_old = nyo[n];
      nyo[n] += ocn_south_ext;
      xocn[n]     = (double *)malloc((nxo[n]+1)*(nyo[n]+1)*sizeof(double));
      yocn[n]     = (double *)malloc((nxo[n]+1)*(nyo[n]+1)*sizeof(double));
      area_ocn[n] = (double *)malloc((nxo[n]  )*(nyo[n]  )*sizeof(double));
      /* assign the global latitude */
      for(j = 0; j < nyo_old+1; j++) for(i = 0; i < nxo[n]+1; i++) {
	xocn[n][(j+ocn_south_ext)*(nxo[n]+1)+i] = tmpx[j*(nxo[n]+1)+i] * D2R;
	yocn[n][(j+ocn_south_ext)*(nxo[n]+1)+i] = tmpy[j*(nxo[n]+1)+i] * D2R;
      }
      if(ocn_south_ext==1) {
	for(i=0; i<nxo[n]+1; i++) {
	  xocn[n][i] = xocn[n][nxo[n]+1+i];
	  yocn[n][i] = min_atm_lat;
	}
      }
      free(tmpx);
      free(tmpy);
    }

    if(clip_method == GREAT_CIRCLE_CLIP) {
      cart_xocn = (double **) malloc( ntile_lnd*sizeof(double *));
      cart_yocn = (double **) malloc( ntile_lnd*sizeof(double *));
      cart_zocn =  (double **) malloc( ntile_ocn*sizeof(double *));
      for(n=0; n<ntile_ocn; n++) {
	cart_xocn[n]     = (double *)malloc((nxo[n]+1)*(nyo[n]+1)*sizeof(double));
	cart_yocn[n]     = (double *)malloc((nxo[n]+1)*(nyo[n]+1)*sizeof(double));
	cart_zocn[n]     = (double *)malloc((nxo[n]+1)*(nyo[n]+1)*sizeof(double));
	latlon2xyz((nxo[n]+1)*(nyo[n]+1), xocn[n], yocn[n], cart_xocn[n], cart_yocn[n], cart_zocn[n]);
	get_grid_great_circle_area(&(nxo[n]), &(nyo[n]), xocn[n], yocn[n], area_ocn[n]);
      }
    }
    else {
      if(mpp_pe()==mpp_root_pe())printf("make_coupler_mosaic: Calculating area_ocn based on lat-lon segments between adjacent grid points.\n");
      for(n=0; n<ntile_ocn; n++) {
        get_grid_global_area(nxo[n], nyo[n], xocn[n], yocn[n], area_ocn[n]);
      }
    }
    mpp_close(m_fid);

    /* read ocean topography */
    omask = (double **)malloc(ntile_ocn*sizeof(double *));
    for(n=0; n<ntile_ocn; n++) {
      char name[128];
      char depth_name[128], mask_name[128];
      int nx, ny, i, j;
      double *depth;
      int mask_name_exist;

      t_fid = mpp_open(otopog, MPP_READ);
      if(n==0) {
         if(mpp_dim_exist(t_fid, "ntiles"))
            ntiles = mpp_get_dimlen(t_fid, "ntiles");
         else
	    ntiles = 1;

         if(ntile_ocn != ntiles) mpp_error("make_coupler_mosaic: dimlen ntiles in mosaic file is not the same as dimlen in topog file");
      }

      if(ntiles == 1)
	strcpy(name, "nx");
      else
	sprintf(name, "nx_tile%d", n+1);
      nx = mpp_get_dimlen(t_fid, name);
      if(ntiles == 1)
	strcpy(name, "ny");
      else
	sprintf(name, "ny_tile%d", n+1);
      ny = mpp_get_dimlen(t_fid, name);
      if( nx != nxo[n] || ny+ocn_south_ext != nyo[n]) mpp_error("make_coupler_mosaic: grid size mismatch between mosaic file and topog file");
      if(ntiles == 1) {
	strcpy(depth_name, "depth");
	strcpy(mask_name, "area_frac");
      }
      else {
	sprintf(depth_name, "depth_tile%d", n+1);
        sprintf(mask_name, "area_frac_tile%d", n+1);
      }
      omask[n] = (double *)malloc(nxo[n]*nyo[n]*sizeof(double));
      for(i=0; i<nxo[n]*nyo[n]; i++) omask[n][i] = 0;
      mask_name_exist = mpp_var_exist(t_fid, mask_name);
      mpp_close(t_fid);
      if(mask_name_exist) {
	if(mpp_pe() == mpp_root_pe() && verbose ) printf("\nNOTE from make_coupler_mosaic: the ocean land/sea mask will be "
					     "determined by field area_frac from file %s\n", otopog);
	get_global_data(otopog, mask_name, nx, ny, omask[n]+ocn_south_ext*nx, 0, 0);
      }
      else {
	if(mpp_pe() == mpp_root_pe() && verbose) printf("\nNOTE from make_coupler_mosaic: the ocean land/sea mask will be "
					     "determined by field depth from file %s\n", otopog);
	depth    = (double *)malloc(nx*ny*sizeof(double));
	get_global_data(otopog, depth_name, nx, ny, depth, 0, 0);
	for(j=0; j<ny; j++) for(i=0; i<nx; i++) {
	  if(depth[j*nx+i] >sea_level) omask[n][(j+ocn_south_ext)*nx+i] = 1;
	}
	free(depth);
      }
    }
  }
  if(print_memory)print_mem_usage("after read ocean grid");

  /* when atm_great_circle_algorithm is 0, lnd_great_circle_algorithm/ocn_great_circle_algorithm must be 0 */
  if( atm_great_circle_algorithm == 0 && (lnd_great_circle_algorithm || ocn_great_circle_algorithm))
    mpp_error("make_coupler_mosaic: when atm does not use great_circle_algorithm, lnd/ocn can not use great_circle_algorithm");

  /*
   * Read wave grid
   */
  if ( wmosaic ) {
    int n, m_fid, g_fid, vid, gid, tid;
    size_t start[4], nread[4];
    char dir[STRING], filename[STRING], file[2*STRING];

    if (strcmp(wmosaic, omosaic) == 0 ) wav_same_as_ocn = 1;

    for(n=0; n<4; n++) {
      start[n] = 0;
      nread[n] = 1;
    }
    m_fid = mpp_open(wmosaic, MPP_READ);
    vid = mpp_get_varid(m_fid, "mosaic");
    mpp_get_var_value(m_fid, vid, wmosaic_name);
    ntile_wav  = mpp_get_dimlen(m_fid, "ntiles");
    nxw        = (int *) malloc (ntile_wav*sizeof(int));
    nyw        = (int *) malloc (ntile_wav*sizeof(int));
    xwav       = (double **) malloc( ntile_wav*sizeof(double *));
    ywav       = (double **) malloc( ntile_wav*sizeof(double *));
    area_wav   = (double **) malloc( ntile_wav*sizeof(double *));
    wtile_name = (char **)malloc(ntile_wav*sizeof(char *));
    /* grid should be located in the same directory of mosaic file */
    get_file_path(wmosaic, dir);
    gid = mpp_get_varid(m_fid, "gridfiles");
    tid = mpp_get_varid(m_fid, "gridtiles");
    for(n=0; n<ntile_wav; n++) {
      int i, j;

      start[0] = n; start[1] = 0; nread[0] = 1; nread[1] = STRING;
      mpp_get_var_value_block(m_fid, gid, start, nread, filename);
      wtile_name[n] = (char *)malloc(STRING*sizeof(char));
      mpp_get_var_value_block(m_fid, tid, start, nread,  wtile_name[n]);
      sprintf(file, "%s/%s", dir, filename);
      g_fid = mpp_open(file, MPP_READ);
      nxw[n] = mpp_get_dimlen(g_fid, "nx");
      nyw[n] = mpp_get_dimlen(g_fid, "ny");
      mpp_close(g_fid);
      if(nxw[n]%x_refine != 0 ) mpp_error("make_coupler_mosaic: wave supergrid x-size can not be divided by x_refine");
      if(nyw[n]%y_refine != 0 ) mpp_error("make_coupler_mosaic: wave supergrid y-size can not be divided by y_refine");
      nxw[n]      /= x_refine;
      nyw[n]      /= y_refine;
      xwav[n]     = (double *)malloc((nxw[n]+1)*(nyw[n]+1)*sizeof(double));
      ywav[n]     = (double *)malloc((nxw[n]+1)*(nyw[n]+1)*sizeof(double));
      get_global_grid(file, nxw[n], nyw[n], x_refine, y_refine, xwav[n], ywav[n]);
      area_wav[n] = (double *)malloc((nxw[n]  )*(nyw[n]  )*sizeof(double));

      /*scale grid from degree to radian, because create_xgrid assume the grid is in radians */
      for(i=0; i<(nxw[n]+1)*(nyw[n]+1); i++) {
	xwav[n][i] *= D2R;
	ywav[n][i] *= D2R;
      }
    }

    if(clip_method == GREAT_CIRCLE_CLIP) {
      cart_xwav = (double **) malloc( ntile_wav*sizeof(double *));
      cart_ywav = (double **) malloc( ntile_wav*sizeof(double *));
      cart_zwav = (double **) malloc( ntile_wav*sizeof(double *));
      for(n=0; n<ntile_wav; n++) {
	cart_xwav[n]     = (double *)malloc((nxw[n]+1)*(nyw[n]+1)*sizeof(double));
	cart_ywav[n]     = (double *)malloc((nxw[n]+1)*(nyw[n]+1)*sizeof(double));
	cart_zwav[n]     = (double *)malloc((nxw[n]+1)*(nyw[n]+1)*sizeof(double));
	latlon2xyz((nxw[n]+1)*(nyw[n]+1), xwav[n], ywav[n], cart_xwav[n], cart_ywav[n], cart_zwav[n]);
	get_grid_great_circle_area(&(nxw[n]), &(nyw[n]), xwav[n], ywav[n], area_wav[n]);
      }
    }
    else {
      for(n=0; n<ntile_wav; n++) {
        get_grid_global_area(nxw[n], nyw[n], xwav[n], ywav[n], area_wav[n]);
      }
    }
    mpp_close(m_fid);
  }

  /* remove longitude and latitude data when clip_method is 'great_circle' */
  if(!print_grid){
    int n;
    if(clip_method == GREAT_CIRCLE_CLIP){
      for(n=0; n<ntile_atm; n++) {
	free(xatm[n]);
	free(yatm[n]);
      }
      free(xatm);
      free(yatm);
      for(n=0; n<ntile_ocn; n++) {
	free(xocn[n]);
	free(yocn[n]);
      }
      free(xocn);
      free(yocn);
      if( xlnd ) {
	for(n=0; n<ntile_lnd; n++) {
	  free(xlnd[n]);
	  free(ylnd[n]);
	}
	free(xlnd);
	free(ylnd);
      }
      if( xwav ) {
	for(n=0; n<ntile_wav; n++) {
	  free(xwav[n]);
	  free(ywav[n]);
	}
	free(xwav);
	free(ywav);
      }
    }
  }
  /* Either omosaic is different from both lmosaic and amosaic,
     or all the three mosaic are the same
  */
  if(strcmp(omosaic_name, amosaic_name)) { /* omosaic is different from amosaic */
    if(!strcmp(omosaic_name, lmosaic_name)) mpp_error("make_coupler_mosaic: omosaic is the same as lmosaic, "
						      "but different from amosaic.");
    same_mosaic = 0;
  }
  else { /* omosaic is same as amosaic */
    if(strcmp(omosaic_name, lmosaic_name)) mpp_error("make_coupler_mosaic: omosaic is the same as amosaic, "
						      "but different from lmosaic.");
    same_mosaic = 1;
  }

  /***************************************************************************************
     First generate the exchange grid between atmos mosaic and land/ocean mosaic
  ***************************************************************************************/
  nfile_axo = 0;
  nfile_axl = 0;
  nfile_lxo = 0;
  int nbad=0;
  {
    int no, nl, na, n;
    size_t  **naxl, **naxo;
    int     ***atmxlnd_ia,   ***atmxlnd_ja,   ***atmxlnd_il,   ***atmxlnd_jl;
    int     ***atmxocn_ia,   ***atmxocn_ja,   ***atmxocn_io,   ***atmxocn_jo;
    double  ***atmxlnd_area, ***atmxlnd_dia,  ***atmxlnd_dja,  ***atmxlnd_dil,  ***atmxlnd_djl;
    double  ***atmxocn_area, ***atmxocn_dia,  ***atmxocn_dja,  ***atmxocn_dio,  ***atmxocn_djo;
    double  ***atmxocn_clon, ***atmxocn_clat, ***atmxlnd_clon, ***atmxlnd_clat;
    double   min_area;
    time_t time_start, time_end;

    naxl         = (size_t ** )malloc(ntile_atm*sizeof(size_t *));
    naxo         = (size_t ** )malloc(ntile_atm*sizeof(size_t *));
    atmxlnd_area = (double ***)malloc(ntile_atm*sizeof(double **));
    atmxlnd_ia   = (int    ***)malloc(ntile_atm*sizeof(int    **));
    atmxlnd_ja   = (int    ***)malloc(ntile_atm*sizeof(int    **));
    atmxlnd_il   = (int    ***)malloc(ntile_atm*sizeof(int    **));
    atmxlnd_jl   = (int    ***)malloc(ntile_atm*sizeof(int    **));
    atmxocn_area = (double ***)malloc(ntile_atm*sizeof(double **));
    atmxocn_ia   = (int    ***)malloc(ntile_atm*sizeof(int    **));
    atmxocn_ja   = (int    ***)malloc(ntile_atm*sizeof(int    **));
    atmxocn_io   = (int    ***)malloc(ntile_atm*sizeof(int    **));
    atmxocn_jo   = (int    ***)malloc(ntile_atm*sizeof(int    **));

    if(interp_order == 2 ) {
      atmxlnd_dia  = (double ***)malloc(ntile_atm*sizeof(double **));
      atmxlnd_dja  = (double ***)malloc(ntile_atm*sizeof(double **));
      atmxlnd_dil  = (double ***)malloc(ntile_atm*sizeof(double **));
      atmxlnd_djl  = (double ***)malloc(ntile_atm*sizeof(double **));
      atmxocn_dia  = (double ***)malloc(ntile_atm*sizeof(double **));
      atmxocn_dja  = (double ***)malloc(ntile_atm*sizeof(double **));
      atmxocn_dio  = (double ***)malloc(ntile_atm*sizeof(double **));
      atmxocn_djo  = (double ***)malloc(ntile_atm*sizeof(double **));
      atmxlnd_clon = (double ***)malloc(ntile_atm*sizeof(double **));
      atmxlnd_clat = (double ***)malloc(ntile_atm*sizeof(double **));
      atmxocn_clon = (double ***)malloc(ntile_atm*sizeof(double **));
      atmxocn_clat = (double ***)malloc(ntile_atm*sizeof(double **));
    }

    for(na=0; na<ntile_atm; na++) {
      naxl[na]         = (size_t * )malloc(ntile_lnd*sizeof(size_t));
      naxo[na]         = (size_t * )malloc(ntile_ocn*sizeof(size_t));
      atmxlnd_area[na] = (double **)malloc(ntile_lnd*sizeof(double *));
      atmxlnd_ia[na]   = (int    **)malloc(ntile_lnd*sizeof(int    *));
      atmxlnd_ja[na]   = (int    **)malloc(ntile_lnd*sizeof(int    *));
      atmxlnd_il[na]   = (int    **)malloc(ntile_lnd*sizeof(int    *));
      atmxlnd_jl[na]   = (int    **)malloc(ntile_lnd*sizeof(int    *));
      atmxocn_area[na] = (double **)malloc(ntile_ocn*sizeof(double *));
      atmxocn_ia[na]   = (int    **)malloc(ntile_ocn*sizeof(int    *));
      atmxocn_ja[na]   = (int    **)malloc(ntile_ocn*sizeof(int    *));
      atmxocn_io[na]   = (int    **)malloc(ntile_ocn*sizeof(int    *));
      atmxocn_jo[na]   = (int    **)malloc(ntile_ocn*sizeof(int    *));

      if(interp_order == 2 ) {
	atmxlnd_dia [na] = (double **)malloc(ntile_lnd*sizeof(double *));
	atmxlnd_dja [na] = (double **)malloc(ntile_lnd*sizeof(double *));
	atmxlnd_dil [na] = (double **)malloc(ntile_lnd*sizeof(double *));
	atmxlnd_djl [na] = (double **)malloc(ntile_lnd*sizeof(double *));
	atmxocn_dia [na] = (double **)malloc(ntile_ocn*sizeof(double *));
	atmxocn_dja [na] = (double **)malloc(ntile_ocn*sizeof(double *));
	atmxocn_dio [na] = (double **)malloc(ntile_ocn*sizeof(double *));
	atmxocn_djo [na] = (double **)malloc(ntile_ocn*sizeof(double *));
	atmxlnd_clon[na] = (double **)malloc(ntile_lnd*sizeof(double *));
	atmxlnd_clat[na] = (double **)malloc(ntile_lnd*sizeof(double *));
	atmxocn_clon[na] = (double **)malloc(ntile_ocn*sizeof(double *));
	atmxocn_clat[na] = (double **)malloc(ntile_ocn*sizeof(double *));
      }

      for(nl=0; nl<ntile_lnd; nl++) {
	atmxlnd_area[na][nl] = (double *)malloc(MAXXGRID*sizeof(double));
	atmxlnd_ia  [na][nl] = (int    *)malloc(MAXXGRID*sizeof(int   ));
	atmxlnd_ja  [na][nl] = (int    *)malloc(MAXXGRID*sizeof(int   ));
	atmxlnd_il  [na][nl] = (int    *)malloc(MAXXGRID*sizeof(int   ));
	atmxlnd_jl  [na][nl] = (int    *)malloc(MAXXGRID*sizeof(int   ));
	if(interp_order == 2 ) {
	  atmxlnd_clon[na][nl] = (double *)malloc(MAXXGRID*sizeof(double));
	  atmxlnd_clat[na][nl] = (double *)malloc(MAXXGRID*sizeof(double));
	  atmxlnd_dia [na][nl] = (double *)malloc(MAXXGRID*sizeof(double));
	  atmxlnd_dja [na][nl] = (double *)malloc(MAXXGRID*sizeof(double));
	  atmxlnd_dil [na][nl] = (double *)malloc(MAXXGRID*sizeof(double));
	  atmxlnd_djl [na][nl] = (double *)malloc(MAXXGRID*sizeof(double));
	}
      }

      for(no=0; no<ntile_ocn; no++) {
	atmxocn_area[na][no] = (double *)malloc(MAXXGRID*sizeof(double));
	atmxocn_ia  [na][no] = (int    *)malloc(MAXXGRID*sizeof(int   ));
	atmxocn_ja  [na][no] = (int    *)malloc(MAXXGRID*sizeof(int   ));
	atmxocn_io  [na][no] = (int    *)malloc(MAXXGRID*sizeof(int   ));
	atmxocn_jo  [na][no] = (int    *)malloc(MAXXGRID*sizeof(int   ));
	if(interp_order == 2 ) {
	  atmxocn_clon[na][no] = (double *)malloc(MAXXGRID*sizeof(double));
	  atmxocn_clat[na][no] = (double *)malloc(MAXXGRID*sizeof(double));
	  atmxocn_dia [na][no] = (double *)malloc(MAXXGRID*sizeof(double));
	  atmxocn_dja [na][no] = (double *)malloc(MAXXGRID*sizeof(double));
	  atmxocn_dio [na][no] = (double *)malloc(MAXXGRID*sizeof(double));
	  atmxocn_djo [na][no] = (double *)malloc(MAXXGRID*sizeof(double));
	}
      }
    }

  if(print_memory)print_mem_usage("before calcuting exchange grid");


    time_start = time(NULL);
    for(na=0; na<ntile_atm; na++) {

      int      k,l, is, ie, js, je, la, ia, ja, il, jl, io, jo, layout[2];
      int      n0, n1, n2, n3, na_in, nl_in, no_in, n_out, n_out2;
      double   xa_min, ya_min, xo_min, yo_min, xl_min, yl_min, xa_avg;
      double   xa_max, ya_max, xo_max, yo_max, xl_max, yl_max;
      double   xarea;
      double   xa[MV], ya[MV], za[MV];
      double   xl[MV], yl[MV], zl[MV];
      double   xo[MV], yo[MV], zo[MV];
      double   x_out[MV], y_out[MV], z_out[MV];
      double   y_out_max, y_out_min;
      double   atmxlnd_x[MX][MV], atmxlnd_y[MX][MV], atmxlnd_z[MX][MV];
      int      atmxlnd_c[MX][MV];
      int      num_v[MX];
      int      axl_i[MX], axl_j[MX], axl_t[MX];
      double   axl_xmin[MX], axl_xmax[MX], axl_ymin[MX], axl_ymax[MX];
      double   axl_area[MX], axl_clon[MX], axl_clat[MX];
      size_t   count;
      int one, intc,atmxlnd_count;
      domain2D Dom;
      double   yy;
      int      *js_lnd, *je_lnd;
      int      *js_ocn, *je_ocn;
      int      *is_lnd, *ie_lnd;
      int      *is_ocn, *ie_ocn;
      char     mesg[256];

      if(print_memory) {
        sprintf(mesg, "start of loop na=%d", na);
        print_mem_usage(mesg);
      }
      for(nl=0; nl<ntile_lnd; nl++) naxl[na][nl] = 0;
      for(no=0; no<ntile_ocn; no++) naxo[na][no] = 0;
      layout[0] = mpp_npes();
      layout[1] = 1;

      mpp_define_domain2d(nxa[na]*nya[na], 1, layout, 0, 0, &Dom);
      mpp_get_compute_domain2d(Dom, &is, &ie, &js, &je );
      /* find the js_ocn, je_ocn, js_lnd, je_lnd
         and is_ocn, ie_ocn, is_lnd, ie_lnd
         In x-direction, cyclic condition will be considered */
      is_lnd = (int *)malloc(ntile_lnd*sizeof(int));
      ie_lnd = (int *)malloc(ntile_lnd*sizeof(int));
      is_ocn = (int *)malloc(ntile_ocn*sizeof(int));
      ie_ocn = (int *)malloc(ntile_ocn*sizeof(int));
      js_lnd = (int *)malloc(ntile_lnd*sizeof(int));
      je_lnd = (int *)malloc(ntile_lnd*sizeof(int));
      js_ocn = (int *)malloc(ntile_ocn*sizeof(int));
      je_ocn = (int *)malloc(ntile_ocn*sizeof(int));

      if( clip_method == GREAT_CIRCLE_CLIP ) { /* we just use is = 0 and ie = nx-1,
						we may change it to improve performace */
	for(nl=0; nl<ntile_lnd; nl++) {
	  is_lnd[nl] = 0;
	  ie_lnd[nl] = nxl[nl]-1;
	  js_lnd[nl] = 0;
	  je_lnd[nl] = nyl[nl]-1;
	}
	for(no=0; no<ntile_ocn; no++) {
	  is_ocn[no] = 0;
	  ie_ocn[no] = nxo[no]-1;
	  js_ocn[no] = 0;
	  je_ocn[no] = nyo[no]-1;
	}
      }
      else {
	ya_min = 9999;
	ya_max = -9999;
	for(la=is;la<=ie;la++) {

	  ia = la%nxa[na];
	  ja = la/nxa[na];
	  n0 = ja    *(nxa[na]+1) + ia;
	  n1 = ja    *(nxa[na]+1) + ia+1;
	  n2 = (ja+1)*(nxa[na]+1) + ia+1;
	  n3 = (ja+1)*(nxa[na]+1) + ia;

	  ya[0] = yatm[na][n0];
	  ya[1] = yatm[na][n1];
	  ya[2] = yatm[na][n2];
	  ya[3] = yatm[na][n3];
	  if(ya[0] > ya_max) ya_max = ya[0];
	  if(ya[1] > ya_max) ya_max = ya[1];
	  if(ya[2] > ya_max) ya_max = ya[2];
	  if(ya[3] > ya_max) ya_max = ya[3];
	  if(ya[0] < ya_min) ya_min = ya[0];
	  if(ya[1] < ya_min) ya_min = ya[1];
	  if(ya[2] < ya_min) ya_min = ya[2];
	  if(ya[3] < ya_min) ya_min = ya[3];
	}

	for(nl=0; nl<ntile_lnd; nl++) {
	  js_lnd[nl] = nyl[nl]; je_lnd[nl] = -1;
	  for(jl = 0; jl <= nyl[nl]; jl ++) for(il = 0; il <= nxl[nl]; il++) {
	    yy = ylnd[nl][jl    *(nxl[nl]+1) + il];
            if( yy > ya_min ) {
               if(jl < js_lnd[nl] ) js_lnd[nl] = jl;
            }
            if( yy < ya_max ) {
               if(jl > je_lnd[nl] ) je_lnd[nl] = jl;
            }
	  }
	  js_lnd[nl] = max(0, js_lnd[nl]-1);
	  je_lnd[nl] = min(nyl[nl]-1, je_lnd[nl]+1);
	  if(nl==na ||  !lnd_same_as_atm ) {
	    is_lnd[nl] = 0;
	    ie_lnd[nl] = nxl[nl]-1;
	  }
	  else {
	    is_lnd[nl] = nxl[nl]-1;
	    ie_lnd[nl] = 0;
	  }

	}

	for(no=0; no<ntile_ocn; no++) {
	  js_ocn[no] = nyo[no]; je_ocn[no] = -1;
	  for(jo = 0; jo <= nyo[no]; jo ++) for(io = 0; io <= nxo[no]; io++) {
	    yy = yocn[no][jo    *(nxo[no]+1) + io];
            if( yy > ya_min ) {
               if(jo < js_ocn[no] ) js_ocn[no] = jo;
            }
            if( yy < ya_max ) {
               if(jo > je_ocn[no] ) je_ocn[no] = jo;
            }
	  }
	  js_ocn[no] = max(0, js_ocn[no]-1);
	  je_ocn[no] = min(nyo[no]-1, je_ocn[no]+1);

	  if(no==na ||  !ocn_same_as_atm ) {
	    is_ocn[no] = 0;
	    ie_ocn[no] = nxo[no]-1;
	  }
	  else {
	    is_ocn[no] = nxo[no]-1;
	    ie_ocn[no] = 0;
	  }

	}
      }

     if(mpp_pe()==mpp_root_pe() && verbose)printf("na = %d, la = %d, is=%d, ie = %d\n", na, la, is, ie);
      atmxlnd_count=0;
      for(la=is;la<=ie;la++) {
	ia = la%nxa[na];
	ja = la/nxa[na];

        if(print_grid) {
          n0 = ja    *(nxa[na]+1) + ia;
          n1 = ja    *(nxa[na]+1) + ia+1;
          n2 = (ja+1)*(nxa[na]+1) + ia+1;
          n3 = (ja+1)*(nxa[na]+1) + ia;
          xa[0] = xatm[na][n0]; ya[0] = yatm[na][n0];
          xa[1] = xatm[na][n1]; ya[1] = yatm[na][n1];
          xa[2] = xatm[na][n2]; ya[2] = yatm[na][n2];
          xa[3] = xatm[na][n3]; ya[3] = yatm[na][n3];
          printf("atm grid is \n");
          printf("%15.11f, %15.11f \n", xa[0]*R2D, ya[0]*R2D);
          printf("%15.11f, %15.11f \n", xa[1]*R2D, ya[1]*R2D);
          printf("%15.11f, %15.11f \n", xa[2]*R2D, ya[2]*R2D);
          printf("%15.11f, %15.11f \n", xa[3]*R2D, ya[3]*R2D);
          printf("%15.11f, %15.11f \n", xa[0]*R2D, ya[0]*R2D);
        }

	if(clip_method == GREAT_CIRCLE_CLIP) { /*clockwise*/
	  n0 = ja    *(nxa[na]+1) + ia;
	  n1 = (ja+1)*(nxa[na]+1) + ia;
	  n2 = (ja+1)*(nxa[na]+1) + ia+1;
	  n3 = ja    *(nxa[na]+1) + ia+1;
	  xa[0] = cart_xatm[na][n0]; ya[0] = cart_yatm[na][n0]; za[0] = cart_zatm[na][n0];
	  xa[1] = cart_xatm[na][n1]; ya[1] = cart_yatm[na][n1]; za[1] = cart_zatm[na][n1];
	  xa[2] = cart_xatm[na][n2]; ya[2] = cart_yatm[na][n2]; za[2] = cart_zatm[na][n2];
	  xa[3] = cart_xatm[na][n3]; ya[3] = cart_yatm[na][n3]; za[3] = cart_zatm[na][n3];
	}
	else {
	  n0 = ja    *(nxa[na]+1) + ia;
	  n1 = ja    *(nxa[na]+1) + ia+1;
	  n2 = (ja+1)*(nxa[na]+1) + ia+1;
	  n3 = (ja+1)*(nxa[na]+1) + ia;
	  xa[0] = xatm[na][n0]; ya[0] = yatm[na][n0];
	  xa[1] = xatm[na][n1]; ya[1] = yatm[na][n1];
	  xa[2] = xatm[na][n2]; ya[2] = yatm[na][n2];
	  xa[3] = xatm[na][n3]; ya[3] = yatm[na][n3];
	  ya_min  = minval_double(4, ya);
	  ya_max  = maxval_double(4, ya);

	  na_in   = fix_lon(xa, ya, 4, M_PI);

	  xa_min  = minval_double(na_in, xa);
	  xa_max  = maxval_double(na_in, xa);
	  xa_avg  = avgval_double(na_in, xa);
	}
	count = 0;
      	for(nl=0; nl<ntile_lnd; nl++) {
       	  for(jl = js_lnd[nl]; jl <= je_lnd[nl]; jl ++) for(il = is_lnd[nl]; il <= ie_lnd[nl]; il++) {

            if(print_grid) {
              n0 = jl    *(nxl[nl]+1) + il;
              n1 = jl    *(nxl[nl]+1) + il+1;
              n2 = (jl+1)*(nxl[nl]+1) + il+1;
              n3 = (jl+1)*(nxl[nl]+1) + il;
              xl[0] = xlnd[nl][n0]; yl[0] = ylnd[nl][n0];
              xl[1] = xlnd[nl][n1]; yl[1] = ylnd[nl][n1];
              xl[2] = xlnd[nl][n2]; yl[2] = ylnd[nl][n2];
              xl[3] = xlnd[nl][n3]; yl[3] = ylnd[nl][n3];

              printf("land grid is \n");
              printf("%15.11f, %15.11f \n", xl[0]*R2D, yl[0]*R2D);
              printf("%15.11f, %15.11f \n", xl[1]*R2D, yl[1]*R2D);
              printf("%15.11f, %15.11f \n", xl[2]*R2D, yl[2]*R2D);
              printf("%15.11f, %15.11f \n", xl[3]*R2D, yl[3]*R2D);
              printf("%15.11f, %15.11f \n", xl[0]*R2D, yl[0]*R2D);
            }

       	   if(clip_method == GREAT_CIRCLE_CLIP) { /* clockwise */
	      n0 = jl    *(nxl[nl]+1) + il;
	      n1 = (jl+1)*(nxl[nl]+1) + il;
	      n2 = (jl+1)*(nxl[nl]+1) + il+1;
	      n3 = jl    *(nxl[nl]+1) + il+1;
	      xl[0] = cart_xlnd[nl][n0]; yl[0] = cart_ylnd[nl][n0]; zl[0] = cart_zlnd[nl][n0];
	      xl[1] = cart_xlnd[nl][n1]; yl[1] = cart_ylnd[nl][n1]; zl[1] = cart_zlnd[nl][n1];
	      xl[2] = cart_xlnd[nl][n2]; yl[2] = cart_ylnd[nl][n2]; zl[2] = cart_zlnd[nl][n2];
	      xl[3] = cart_xlnd[nl][n3]; yl[3] = cart_ylnd[nl][n3]; zl[3] = cart_zlnd[nl][n3];
	      if(lnd_same_as_atm) {
		if(na==nl && ja==jl && ia==il) {
		  n_out = 4;
		  for(n=0; n<n_out; n++) {
		    x_out[n] = xl[n];
		    y_out[n] = yl[n];
		    z_out[n] = zl[n];
		  }
		}
		else
		  n_out = 0;
	      }
              else {
		n_out = clip_2dx2d_great_circle(xa, ya, za, 4, xl, yl, zl, 4,
						x_out, y_out, z_out);
	      }
	   }
	    else {
	      n0 = jl    *(nxl[nl]+1) + il;
	      n1 = jl    *(nxl[nl]+1) + il+1;
	      n2 = (jl+1)*(nxl[nl]+1) + il+1;
	      n3 = (jl+1)*(nxl[nl]+1) + il;
	      xl[0] = xlnd[nl][n0]; yl[0] = ylnd[nl][n0];
	      xl[1] = xlnd[nl][n1]; yl[1] = ylnd[nl][n1];
	      xl[2] = xlnd[nl][n2]; yl[2] = ylnd[nl][n2];
	      xl[3] = xlnd[nl][n3]; yl[3] = ylnd[nl][n3];
	      yl_min = minval_double(4, yl);
	      yl_max = maxval_double(4, yl);
	      if(yl_min >= ya_max || yl_max <= ya_min ) continue;
	      nl_in  = fix_lon(xl, yl, 4, xa_avg);
	      xl_min = minval_double(nl_in, xl);
	      xl_max = maxval_double(nl_in, xl);
	      /* xl should in the same range as xa after lon_fix, so no need to
		 consider cyclic condition
	      */

	      if(xa_min >= xl_max || xa_max <= xl_min ) continue;
	      if(lnd_same_as_atm) {
	        if(na==nl && ja==jl && ia==il) {
		  if(na_in != nl_in){
                     printf("Error: LND and AM grids are the same but na_in != nl_in, n_out,%d,%d,%d",na_in, nl_in, n_out);
                     mpp_error("make_coupler_mosaic: inconsistent number of grid box coreners");
                  }
		  n_out = nl_in;
		  for(n=0; n<n_out; n++) {
		    x_out[n] = xl[n];
		    y_out[n] = yl[n];
		  }
		 }
		 else
		   n_out = 0;
	      }
              else {
	        n_out = clip_2dx2d( xa, ya, na_in, xl, yl, nl_in, x_out, y_out );
	      }
	   }

	    if (  n_out > 0 ) {
	      if(clip_method == GREAT_CIRCLE_CLIP)
		xarea=great_circle_area ( n_out, x_out, y_out, z_out);
	      else
		xarea = poly_area(x_out, y_out, n_out);
	      min_area = min(area_lnd[nl][jl*nxl[nl]+il], area_atm[na][la]);
	      y_out_min = minval_double(n_out, y_out);
	      y_out_max = maxval_double(n_out, y_out);
	      if(fabs(y_out_min+0.5*M_PI) < 0.00003 && verbose) {
		printf("Near South Pole ATMxLND grid cell,  ATMxLND_area/ATM_area =%f, ATM_area=%f, LND_area=%f \n", xarea/area_atm[na][la], area_atm[na][la], area_lnd[nl][jl*nxl[nl]+il]);
                printf("longitudes and latitudes of the %d corners are:\n", n_out);
		for(n=0; n<n_out; n++) printf("%7.3f, ", x_out[n]*R2D);
                printf("\n");
		for(n=0; n<n_out; n++) printf("%7.3f, ", y_out[n]*R2D);
                printf("\n");
	      }
	      if(fabs(y_out_max-0.5*M_PI) < 0.00003 && verbose) {
		printf("Near North Pole ATMxLND grid cell,  ATMxLND_area/ATM_area =%f, ATM_area=%f, LND_area=%f\n", xarea/area_atm[na][la], area_atm[na][la], area_lnd[nl][jl*nxl[nl]+il]);
                printf("longitudes and latitudes of the %d corners are:\n", n_out);
		for(n=0; n<n_out; n++) printf("%7.3f, ", x_out[n]*R2D);
                printf("\n");
		for(n=0; n<n_out; n++) printf("%7.3f, ", y_out[n]*R2D);
                printf("\n");
	      }
	      if( xarea/min_area > area_ratio_thresh ) {

                if(print_grid) {
                  double xtmp[20],ytmp[20];
                  printf("n_axl is %d\n", n_out);
                  /* convert to lon-lat */
                  xyz2latlon(n_out, x_out, y_out, z_out, xtmp, ytmp);
                  for(n=0; n<n_out; n++) printf("%15.11f, %15.11f \n", xtmp[n]*R2D, ytmp[n]*R2D);
                }

		axl_i[count]    = il;
		axl_j[count]    = jl;
		axl_t[count]    = nl;
		num_v[count]    = n_out;
		axl_area[count] = 0;
		if(interp_order == 2) {
		  axl_clon[count] = 0;
		  axl_clat[count] = 0;
		}

		/*  remember the exchange grid vertices */
		if( clip_method == GREAT_CIRCLE_CLIP) {
		  for(n=0; n<n_out; n++) {
		    atmxlnd_x[count][n] = x_out[n];
		    atmxlnd_y[count][n] = y_out[n];
		    atmxlnd_z[count][n] = z_out[n];
		  }
		}
		else {
		  for(n=0; n<n_out; n++) {
		    atmxlnd_x[count][n] = x_out[n];
		    atmxlnd_y[count][n] = y_out[n];
		  }
		  axl_xmin[count] = minval_double(n_out, x_out);
		  axl_xmax[count] = maxval_double(n_out, x_out);
		  axl_ymin[count] = minval_double(n_out, y_out);
		  axl_ymax[count] = maxval_double(n_out, y_out);
		}

		++count;
		if(count>MX) mpp_error("make_coupler_mosaic: count is greater than MX, increase MX");
      	      }/*if( xarea/min_area > area_ratio_thresh )*/
	    }/*if(nout>0)*/
	  }/*ni,nj loop*/
	}/*nl loop*/
        atmxlnd_count=atmxlnd_count+count;
	/* calculate atmos/ocean x-cells */
	for(no=0; no<ntile_ocn; no++) {
	  for(jo = js_ocn[no]; jo <= je_ocn[no]; jo++) for(io = is_ocn[no]; io <= ie_ocn[no]; io++) {
	    double ocn_frac, lnd_frac;

            if(print_grid) {
              n0 = jo    *(nxo[no]+1) + io;
              n1 = jo    *(nxo[no]+1) + io+1;
              n2 = (jo+1)*(nxo[no]+1) + io+1;
              n3 = (jo+1)*(nxo[no]+1) + io;
              xo[0] = xocn[no][n0]; yo[0] = yocn[no][n0];
              xo[1] = xocn[no][n1]; yo[1] = yocn[no][n1];
              xo[2] = xocn[no][n2]; yo[2] = yocn[no][n2];
              xo[3] = xocn[no][n3]; yo[3] = yocn[no][n3];
              printf("ocean grid is \n");
              printf("%15.11f, %15.11f \n", xo[0]*R2D+360, yo[0]*R2D);
              printf("%15.11f, %15.11f \n", xo[1]*R2D+360, yo[1]*R2D);
              printf("%15.11f, %15.11f \n", xo[2]*R2D+360, yo[2]*R2D);
              printf("%15.11f, %15.11f \n", xo[3]*R2D+360, yo[3]*R2D);
              printf("%15.11f, %15.11f \n", xo[0]*R2D+360, yo[0]*R2D);
            }


     	    if(clip_method == GREAT_CIRCLE_CLIP) { /* clockwise */
	      n0 = jo    *(nxo[no]+1) + io;
	      n1 = (jo+1)*(nxo[no]+1) + io;
	      n2 = (jo+1)*(nxo[no]+1) + io+1;
	      n3 = jo    *(nxo[no]+1) + io+1;
	      xo[0] = cart_xocn[no][n0]; yo[0] = cart_yocn[no][n0]; zo[0] = cart_zocn[no][n0];
	      xo[1] = cart_xocn[no][n1]; yo[1] = cart_yocn[no][n1]; zo[1] = cart_zocn[no][n1];
	      xo[2] = cart_xocn[no][n2]; yo[2] = cart_yocn[no][n2]; zo[2] = cart_zocn[no][n2];
	      xo[3] = cart_xocn[no][n3]; yo[3] = cart_yocn[no][n3]; zo[3] = cart_zocn[no][n3];
	    }
	    else {
	      n0 = jo    *(nxo[no]+1) + io;
	      n1 = jo    *(nxo[no]+1) + io+1;
	      n2 = (jo+1)*(nxo[no]+1) + io+1;
	      n3 = (jo+1)*(nxo[no]+1) + io;
	      xo[0] = xocn[no][n0]; yo[0] = yocn[no][n0];
	      xo[1] = xocn[no][n1]; yo[1] = yocn[no][n1];
	      xo[2] = xocn[no][n2]; yo[2] = yocn[no][n2];
	      xo[3] = xocn[no][n3]; yo[3] = yocn[no][n3];
	      yo_min = minval_double(4, yo);
	      yo_max = maxval_double(4, yo);
	      no_in  = fix_lon(xo, yo, 4, xa_avg);
	      xo_min = minval_double(no_in, xo);
	      xo_max = maxval_double(no_in, xo);
	    }
	    ocn_frac = omask[no][jo*nxo[no]+io];
	    lnd_frac = 1 - ocn_frac;

	    if(ocn_frac > MIN_AREA_FRAC) { /* over sea/ice */
	      /* xo should in the same range as xa after lon_fix, so no need to
		 consider cyclic condition
	      */

	      if( clip_method == GREAT_CIRCLE_CLIP ) {
		n_out = clip_2dx2d_great_circle(xa, ya, za, 4, xo, yo, zo, 4,
						x_out, y_out, z_out);
	      }
	      else {
		if(xa_min >= xo_max || xa_max <= xo_min || yo_min >= ya_max || yo_max <= ya_min ) continue;
		n_out = clip_2dx2d( xa, ya, na_in, xo, yo, no_in, x_out, y_out );
	      }
	      if(fabs(y_out_max-0.5*M_PI) < 0.00003 && fabs(yo_max-0.5*M_PI) < 0.00003 && verbose) {
		printf("Near North Pole ATMxLND grid cell\n");
		for(n=0; n<na_in; n++) printf("%7.3f, ", xa[n]*R2D);
                printf("\n");
		for(n=0; n<na_in; n++) printf("%7.3f, ", ya[n]*R2D);
                printf("\n");
		printf("Near North Pole OCN grid cell\n");
		for(n=0; n<no_in; n++) printf("%7.3f, ", xo[n]*R2D);
                printf("\n");
		for(n=0; n<no_in; n++) printf("%7.3f, ", yo[n]*R2D);
                printf("\n");
		printf("Near North Pole ATMxLNDxOCN grid cell\n");
		for(n=0; n<n_out; n++) printf("%7.3f, ", x_out[n]*R2D);
                printf("\n");
		for(n=0; n<n_out; n++) printf("%7.3f, ", y_out[n]*R2D);
                printf("\n");
	      }
	      if (  n_out > 0) {


		if( clip_method == GREAT_CIRCLE_CLIP )
		  xarea=great_circle_area ( n_out, x_out, y_out, z_out)*ocn_frac;
		else
		  xarea = poly_area(x_out, y_out, n_out )*ocn_frac;

		if(xarea<0) printf("error: xarea<0, %f",xarea);
		min_area = min(area_ocn[no][jo*nxo[no]+io], area_atm[na][la]);
	        if(fabs(y_out_max-0.5*M_PI) < 0.00003 && fabs(yo_max-0.5*M_PI) < 0.00003 && verbose) {
	           printf("Near North Pole: ATMxLNDxOCN_area/ATM_area =%f \n",xarea/min_area);}
		if(xarea/min_area > area_ratio_thresh) {

		  atmxocn_area[na][no][naxo[na][no]] = xarea;
		  atmxocn_io[na][no][naxo[na][no]]   = io;
		  atmxocn_jo[na][no][naxo[na][no]]   = jo;
		  atmxocn_ia[na][no][naxo[na][no]]   = ia;
		  atmxocn_ja[na][no][naxo[na][no]]   = ja;
		  if(interp_order == 2) {
		    atmxocn_clon[na][no][naxo[na][no]] = poly_ctrlon ( x_out, y_out, n_out, xa_avg)*ocn_frac;
		    atmxocn_clat[na][no][naxo[na][no]] = poly_ctrlat ( x_out, y_out, n_out )*ocn_frac;
		  }
		  ++(naxo[na][no]);
		  if(naxo[na][no] > MAXXGRID) mpp_error("naxo is greater than MAXXGRID, increase MAXXGRID");
		}
	      }
	    }
	    if(lnd_frac > MIN_AREA_FRAC) { /* over land */
	      /* find the overlap of atmxlnd and ocean cell */
	      for(l=0; l<count; l++) {
		if( clip_method == GREAT_CIRCLE_CLIP )
		  n_out = clip_2dx2d_great_circle(atmxlnd_x[l], atmxlnd_y[l], atmxlnd_z[l], num_v[l], xo, yo, zo, 4,
						  x_out, y_out, z_out);
		else {
		  if(axl_xmin[l] >= xo_max || axl_xmax[l] <= xo_min || axl_ymin[l] >= ya_max || axl_ymax[l] <= ya_min ) continue;
		  n_out = clip_2dx2d( atmxlnd_x[l], atmxlnd_y[l], num_v[l], xo, yo, no_in, x_out, y_out );
		}
		if( n_out > 0) {
		  if( clip_method == GREAT_CIRCLE_CLIP )
		    xarea=great_circle_area ( n_out, x_out, y_out, z_out)*lnd_frac;
		  else
		    xarea = poly_area(x_out, y_out, n_out )*lnd_frac;
		  min_area = min(area_lnd[axl_t[l]][axl_j[l]*nxl[axl_t[l]]+axl_i[l]], area_atm[na][la]);
		  if(xarea/min_area > area_ratio_thresh) {

                    if(print_grid) {
                      double xtmp[20],ytmp[20];
                      printf("num exchange grid between ocean and axl is %d\n", n_out);
                      /* convert to lon-lat */
                      xyz2latlon(n_out, x_out, y_out, z_out, xtmp, ytmp);

                      for(n=0; n<n_out; n++) printf("%15.11f, %15.11f \n", xtmp[n]*R2D, ytmp[n]*R2D);
                    }

		    axl_area[l] += xarea;
		    if(interp_order == 2) {
		      axl_clon[l] += poly_ctrlon ( x_out, y_out, n_out, xa_avg)*lnd_frac;
		      axl_clat[l] += poly_ctrlat ( x_out, y_out, n_out)*lnd_frac;
		    }
		  }
		}
	      }//for(l=0; l<count; l++)
	    }//if(lnd_frac > MIN_AREA_FRAC)
	  }//for(jo = js_ocn[no]; jo <= je_ocn[no]; jo++) for(io = is_ocn[no]; io <= ie_ocn[no]; io++)
	}//for(no=0; no<ntile_ocn; no++)
	/* get the exchange grid between land and atmos. */
	for(l=0; l<count; l++) {
	  nl = axl_t[l];
	  min_area = min(area_lnd[nl][axl_j[l]*nxl[nl]+axl_i[l]], area_atm[na][la]);
	  if(fabs(axl_ymin[l]+0.5*M_PI) < 0.00003 && verbose){
	    printf("Near South Pole: ATMxLNDxOCN_area/ATM_area =%f \n",axl_area[l]/min_area);}
	  if(axl_area[l]/min_area > area_ratio_thresh) {
	    atmxlnd_area[na][nl][naxl[na][nl]] = axl_area[l];
	    atmxlnd_ia  [na][nl][naxl[na][nl]] = ia;
	    atmxlnd_ja  [na][nl][naxl[na][nl]] = ja;
	    atmxlnd_il  [na][nl][naxl[na][nl]] = axl_i[l];
	    atmxlnd_jl  [na][nl][naxl[na][nl]] = axl_j[l];
	    if(interp_order == 2) {
	      atmxlnd_clon[na][nl][naxl[na][nl]] = axl_clon[l];
	      atmxlnd_clat[na][nl][naxl[na][nl]] = axl_clat[l];
	    }
	    ++(naxl[na][nl]);
	    if(naxl[na][nl] > MAXXGRID) mpp_error("naxl is greater than MAXXGRID, increase MAXXGRID");
	  }
	}

      }/* end of la loop */
      mpp_sum_int(1, &atmxlnd_count);
      if(mpp_pe()==mpp_root_pe() && verbose)printf("Number of atmxlnd cells for ATM tile %d is %i .\n",na+1,atmxlnd_count);

      mpp_delete_domain2d(&Dom);
      free(is_lnd);
      free(ie_lnd);
      free(js_lnd);
      free(je_lnd);
      free(is_ocn);
      free(ie_ocn);
      free(js_ocn);
      free(je_ocn);
      if(print_memory) {
        sprintf(mesg, "end of loop na=%d", na);
        print_mem_usage(mesg);
      }

    } /* end of na loop */
   if(print_memory)print_mem_usage("after calcuting exchange grid");
    time_end = time(NULL);
    if(verbose) printf("one pe %d, The loop used %f seconds.\n", mpp_pe(), difftime(time_end, time_start));
    /* calculate the centroid of model grid, as well as land_mask and ocean_mask */
    {
      double **l_area, **o_area;
      int    nl, no, ll, lo;
      l_area = (double **)malloc(ntile_lnd*sizeof(double *));
      o_area = (double **)malloc(ntile_ocn*sizeof(double *));
      for(nl =0; nl<ntile_lnd; nl++) {
	l_area[nl] = (double *)malloc(nxl[nl]*nyl[nl]*sizeof(double));
	for(ll=0; ll<nxl[nl]*nyl[nl]; ll++) {
	  l_area[nl][ll] = 0;
	}
      }
      for(no =0; no<ntile_ocn; no++) {
	o_area[no] = (double *)malloc(nxo[no]*nyo[no]*sizeof(double));
	for(lo=0; lo<nxo[no]*nyo[no]; lo++) {
	  o_area[no][lo] = 0;
	}
      }

      if(interp_order == 1) {
	for(na=0; na<ntile_atm; na++) {
	   /* could not add the exchange grid between land and nest atmosphere */
	  if( na != tile_nest ) {
	    for(nl=0; nl<ntile_lnd; nl++) {
	      int nxgrid;

	      nxgrid = naxl[na][nl];
	      mpp_sum_int(1, &nxgrid);
	      if(nxgrid > 0) {
		double *g_area;
		int    *g_il, *g_jl;
		int    ii;
		g_il = (int    *)malloc(nxgrid*sizeof(int   ));
		g_jl = (int    *)malloc(nxgrid*sizeof(int   ));
		g_area = (double *)malloc(nxgrid*sizeof(double));
		mpp_gather_field_int   (naxl[na][nl], atmxlnd_il[na][nl], g_il);
		mpp_gather_field_int   (naxl[na][nl], atmxlnd_jl[na][nl], g_jl);
		mpp_gather_field_double(naxl[na][nl], atmxlnd_area[na][nl], g_area);
		for(i=0; i<nxgrid; i++) {
		  ii = g_jl[i]*nxl[nl]+g_il[i];
		  l_area[nl][ii] += g_area[i];
		}
		free(g_il);
		free(g_jl);
		free(g_area);
	      }
	    }
	  }
	  /* could not add the exchange grid between ocean and nest atmosphere */
	  if( na != tile_nest ) {
	    for(no=0; no<ntile_ocn; no++) {
	      int nxgrid;
	      nxgrid = naxo[na][no];
	      mpp_sum_int(1, &nxgrid);
	      if(nxgrid > 0) {
		double *g_area;
		int    *g_io, *g_jo;
		int    ii;
		g_io = (int    *)malloc(nxgrid*sizeof(int   ));
		g_jo = (int    *)malloc(nxgrid*sizeof(int   ));
		g_area = (double *)malloc(nxgrid*sizeof(double));
		mpp_gather_field_int   (naxo[na][no], atmxocn_io[na][no], g_io);
		mpp_gather_field_int   (naxo[na][no], atmxocn_jo[na][no], g_jo);
		mpp_gather_field_double(naxo[na][no], atmxocn_area[na][no], g_area);
		for(i=0; i<nxgrid; i++) {
		  ii = g_jo[i]*nxo[no]+g_io[i];
		  o_area[no][ii] += g_area[i];
		}
		free(g_io);
		free(g_jo);
		free(g_area);
	      }
	    }
	  }
	}
      }
      else { /* interp_order == 2 */
	double **l_clon, **l_clat;
	double **o_clon, **o_clat;
	double  *a_area,  *a_clon,  *a_clat;
	int la;

	l_clon = (double **)malloc(ntile_lnd*sizeof(double *));
	l_clat = (double **)malloc(ntile_lnd*sizeof(double *));
	for(nl =0; nl<ntile_lnd; nl++) {
	  l_clon[nl] = (double *)malloc(nxl[nl]*nyl[nl]*sizeof(double));
	  l_clat[nl] = (double *)malloc(nxl[nl]*nyl[nl]*sizeof(double));
	  for(ll=0; ll<nxl[nl]*nyl[nl]; ll++) {
	    l_clon[nl][ll] = 0;
	    l_clat[nl][ll] = 0;
	  }
	}
	o_clon = (double **)malloc(ntile_ocn*sizeof(double *));
	o_clat = (double **)malloc(ntile_ocn*sizeof(double *));
	for(no =0; no<ntile_ocn; no++) {
	  o_clon[no] = (double *)malloc(nxo[no]*nyo[no]*sizeof(double));
	  o_clat[no] = (double *)malloc(nxo[no]*nyo[no]*sizeof(double));
	  for(lo=0; lo<nxo[no]*nyo[no]; lo++) {
	    o_clon[no][lo] = 0;
	    o_clat[no][lo] = 0;
	  }
	}
	for(na=0; na<ntile_atm; na++) {
	  //	double *area, *clon, *clat;

	  a_area = (double *)malloc(nxa[na]*nya[na]*sizeof(double));
	  a_clon = (double *)malloc(nxa[na]*nya[na]*sizeof(double));
	  a_clat = (double *)malloc(nxa[na]*nya[na]*sizeof(double));
	  for(la=0; la<nxa[na]*nya[na]; la++) {
	    a_area[la] = 0;
	    a_clon[la] = 0;
	    a_clat[la] = 0;
	  }
	   /* could not add the exchange grid between land and nest atmosphere */
	  if( na != tile_nest ) {
	    for(nl=0; nl<ntile_lnd; nl++) {
	      int nxgrid;

	      nxgrid = naxl[na][nl];
	      mpp_sum_int(1, &nxgrid);
	      if(nxgrid > 0) {
		double *g_area, *g_clon, *g_clat;
		int    *g_ia,   *g_ja,   *g_il, *g_jl;
		int    ii;
		g_ia = (int    *)malloc(nxgrid*sizeof(int   ));
		g_ja = (int    *)malloc(nxgrid*sizeof(int   ));
		g_il = (int    *)malloc(nxgrid*sizeof(int   ));
		g_jl = (int    *)malloc(nxgrid*sizeof(int   ));
		g_area = (double *)malloc(nxgrid*sizeof(double));
		g_clon = (double *)malloc(nxgrid*sizeof(double));
		g_clat = (double *)malloc(nxgrid*sizeof(double));
		mpp_gather_field_int   (naxl[na][nl], atmxlnd_ia[na][nl], g_ia);
		mpp_gather_field_int   (naxl[na][nl], atmxlnd_ja[na][nl], g_ja);
		mpp_gather_field_int   (naxl[na][nl], atmxlnd_il[na][nl], g_il);
		mpp_gather_field_int   (naxl[na][nl], atmxlnd_jl[na][nl], g_jl);
		mpp_gather_field_double(naxl[na][nl], atmxlnd_area[na][nl], g_area);
		mpp_gather_field_double(naxl[na][nl], atmxlnd_clon[na][nl], g_clon);
		mpp_gather_field_double(naxl[na][nl], atmxlnd_clat[na][nl], g_clat);
		for(i=0; i<nxgrid; i++) {
		  ii = g_ja[i]*nxa[na]+g_ia[i];
		  a_area[ii] += g_area[i];
		  a_clon[ii] += g_clon[i];
		  a_clat[ii] += g_clat[i];
		  ii = g_jl[i]*nxl[nl]+g_il[i];
		  l_area[nl][ii] += g_area[i];
		  l_clon[nl][ii] += g_clon[i];
		  l_clat[nl][ii] += g_clat[i];
		}
		free(g_ia);
		free(g_ja);
		free(g_il);
		free(g_jl);
		free(g_area);
		free(g_clon);
		free(g_clat);
	      }
	    }
	  }

	  /* could not add the exchange grid between ocean and nest atmosphere */
	  if( na != tile_nest ) {
	    for(no=0; no<ntile_ocn; no++) {
	      int nxgrid;
	      nxgrid = naxo[na][no];
	      mpp_sum_int(1, &nxgrid);
	      if(nxgrid > 0) {
		double *g_area, *g_clon, *g_clat;
		int    *g_ia,   *g_ja,   *g_io, *g_jo;
		int    ii;
		g_ia = (int    *)malloc(nxgrid*sizeof(int   ));
		g_ja = (int    *)malloc(nxgrid*sizeof(int   ));
		g_io = (int    *)malloc(nxgrid*sizeof(int   ));
		g_jo = (int    *)malloc(nxgrid*sizeof(int   ));
		g_area = (double *)malloc(nxgrid*sizeof(double));
		g_clon = (double *)malloc(nxgrid*sizeof(double));
		g_clat = (double *)malloc(nxgrid*sizeof(double));
		mpp_gather_field_int   (naxo[na][no], atmxocn_ia[na][no], g_ia);
		mpp_gather_field_int   (naxo[na][no], atmxocn_ja[na][no], g_ja);
		mpp_gather_field_int   (naxo[na][no], atmxocn_io[na][no], g_io);
		mpp_gather_field_int   (naxo[na][no], atmxocn_jo[na][no], g_jo);
		mpp_gather_field_double(naxo[na][no], atmxocn_area[na][no], g_area);
		mpp_gather_field_double(naxo[na][no], atmxocn_clon[na][no], g_clon);
		mpp_gather_field_double(naxo[na][no], atmxocn_clat[na][no], g_clat);
		for(i=0; i<nxgrid; i++) {
		  ii = g_ja[i]*nxa[na]+g_ia[i];
		  a_area[ii] += g_area[i];
		  a_clon[ii] += g_clon[i];
		  a_clat[ii] += g_clat[i];
		  ii = g_jo[i]*nxo[no]+g_io[i];
		  o_area[no][ii] += g_area[i];
		  o_clon[no][ii] += g_clon[i];
		  o_clat[no][ii] += g_clat[i];
		}
		free(g_ia);
		free(g_ja);
		free(g_io);
		free(g_jo);
		free(g_area);
		free(g_clon);
		free(g_clat);
	      }
	    }
	  }
	  for(la=0; la<nxa[na]*nya[na]; la++) {
	    if(a_area[la] > 0) {
	      a_clon[la] /= a_area[la];
	      a_clat[la] /= a_area[la];
	    }
	  }

	  /* substract atmos centroid to get the centroid distance between atmos grid and exchange grid. */
	  for(nl=0; nl<ntile_lnd; nl++) {
	    for(i=0; i<naxl[na][nl]; i++) {
	      la = atmxlnd_ja[na][nl][i]*nxa[na] + atmxlnd_ia[na][nl][i];
	      atmxlnd_dia[na][nl][i] = atmxlnd_clon[na][nl][i]/atmxlnd_area[na][nl][i] - a_clon[la];
	      atmxlnd_dja[na][nl][i] = atmxlnd_clat[na][nl][i]/atmxlnd_area[na][nl][i] - a_clat[la];
	    }
	  }
	  for(no=0; no<ntile_ocn; no++) {
	    for(i=0; i<naxo[na][no]; i++) {
	      la = atmxocn_ja[na][no][i]*nxa[na] + atmxocn_ia[na][no][i];
	      atmxocn_dia[na][no][i] = atmxocn_clon[na][no][i]/atmxocn_area[na][no][i] - a_clon[la];
	      atmxocn_dja[na][no][i] = atmxocn_clat[na][no][i]/atmxocn_area[na][no][i] - a_clat[la];
	    }
	  }

	  free(a_area);
	  free(a_clon);
	  free(a_clat);
	}


	/* centroid distance from exchange grid to land grid */
	for(nl=0; nl<ntile_lnd; nl++) {
	  for(ll=0; ll<nxl[nl]*nyl[nl]; ll++) {
	    if(l_area[nl][ll] > 0) {
	      l_clon[nl][ll] /= l_area[nl][ll];
	      l_clat[nl][ll] /= l_area[nl][ll];
	    }
	  }
	  for(na=0; na<ntile_atm; na++) {
	    for(i=0; i<naxl[na][nl]; i++) {
	      ll = atmxlnd_jl[na][nl][i]*nxl[nl] + atmxlnd_il[na][nl][i];
	      atmxlnd_dil[na][nl][i] = atmxlnd_clon[na][nl][i]/atmxlnd_area[na][nl][i] - l_clon[nl][ll];
	      atmxlnd_djl[na][nl][i] = atmxlnd_clat[na][nl][i]/atmxlnd_area[na][nl][i] - l_clat[nl][ll];
	    }
	  }
	  free(l_clon[nl]);
	  free(l_clat[nl]);
	}

	/* centroid distance from exchange grid to ocean grid */
	for(no=0; no<ntile_ocn; no++) {
	  for(lo=0; lo<nxo[no]*nyo[no]; lo++) {
	    if(o_area[no][lo] > 0) {
	      o_clon[no][lo] /= o_area[no][lo];
	      o_clat[no][lo] /= o_area[no][lo];
	    }
	  }
	  for(na=0; na<ntile_atm; na++) {
	    for(i=0; i<naxo[na][no]; i++) {
	      lo = atmxocn_jo[na][no][i]*nxo[no] + atmxocn_io[na][no][i];
	      atmxocn_dio[na][no][i] = atmxocn_clon[na][no][i]/atmxocn_area[na][no][i] - o_clon[no][lo];
	      atmxocn_djo[na][no][i] = atmxocn_clat[na][no][i]/atmxocn_area[na][no][i] - o_clat[no][lo];
	    }
	  }
	  free(o_clon[no]);
	  free(o_clat[no]);
	}
	free(o_clon);
	free(o_clat);
	free(l_clon);
	free(l_clat);
      }

      /* calculate ocean_frac and compare ocean_frac with omask */
      /* also write out ocn_frac */
      if(mpp_pe()==mpp_root_pe()){
	int    io, jo;
	double ocn_frac;
	int    id_mask, fid, dims[2];
	int    id_areao, id_areax;
	char ocn_mask_file[STRING];
	double *mask;
	double *areao,*areax;
	int ny;

	for(no=0; no<ntile_ocn; no++) {
	  ny = nyo[no]-ocn_south_ext;
	  mask = (double *)malloc(nxo[no]*ny*sizeof(double));
	  areao = (double *)malloc(nxo[no]*ny*sizeof(double));
	  areax = (double *)malloc(nxo[no]*ny*sizeof(double));
	  for(jo=0; jo<ny; jo++) for(io=0; io<nxo[no]; io++) {
	    i = (jo+ocn_south_ext)*nxo[no]+io;
	    ocn_frac = o_area[no][i]/area_ocn[no][i];
	    if( fabs(omask[no][i] - ocn_frac) > TOLORENCE ) {
              nbad++;
	      printf("at ocean point (%d,%d), omask = %f, ocn_frac = %f, diff = %f\n",
		     io, jo, omask[no][i], ocn_frac, omask[no][i] - ocn_frac);
	    }
	    mask[jo*nxo[no]+io] = ocn_frac;
	    areao[jo*nxo[no]+io] = area_ocn[no][i];
	    areax[jo*nxo[no]+io] = o_area[no][i];
	  }
	  if(ntile_ocn > 1)
	    sprintf(ocn_mask_file, "ocean_mask_tile%d.nc", no+1);
	  else
	    strcpy(ocn_mask_file, "ocean_mask.nc");

	  fid = mpp_open(ocn_mask_file, MPP_WRITE);
    print_provenance_gv_gca(fid, history, grid_version, clip_method == GREAT_CIRCLE_CLIP );

	  dims[1] = mpp_def_dim(fid, "nx", nxo[no]);
	  dims[0] = mpp_def_dim(fid, "ny", ny);
	  id_mask = mpp_def_var(fid, "mask", MPP_DOUBLE, 2, dims,  2, "standard_name",
				"ocean fraction at T-cell centers", "units", "none");
	  id_areao = mpp_def_var(fid, "areaO", MPP_DOUBLE, 2, dims,  2, "standard_name",
				"ocean grid area", "units", "none");
	  id_areax = mpp_def_var(fid, "areaX", MPP_DOUBLE, 2, dims,  2, "standard_name",
				"ocean exchange grid area", "units", "none");
	  mpp_end_def(fid);
	  mpp_put_var_value(fid, id_mask, mask);
	  mpp_put_var_value(fid, id_areao, areao);
	  mpp_put_var_value(fid, id_areax, areax);
	  mpp_close(fid);
	  free(mask);
	}
        if(nbad>0) {
          printf("make_coupler_mosaic: number of points with omask != ofrac is %d\n", nbad);
        }
      }

      /* calculate land_frac and  write out land_frac */
      {
	int    il, jl;
	int    id_mask,id_a_l,id_l_a,id_a_a, fid, dims[2];
	char lnd_mask_file[STRING];
	double *mask,*l_a,*a_l,*a_a;

	for(nl=0; nl<ntile_lnd; nl++) {
	  mask = (double *)malloc(nxl[nl]*nyl[nl]*sizeof(double));
	  l_a = (double *)malloc(nxl[nl]*nyl[nl]*sizeof(double));
	  a_l = (double *)malloc(nxl[nl]*nyl[nl]*sizeof(double));
	  a_a = (double *)malloc(nxl[nl]*nyl[nl]*sizeof(double));
	  for(jl=0; jl<nyl[nl]; jl++) for(il=0; il<nxl[nl]; il++) {
	    i = jl*nxl[nl]+il;
            mask[i] = l_area[nl][i]/area_lnd[nl][i];
            //Land mask being greater than 1 is meaningless and causes trouble
	    //mask[i] = min(l_area[nl][i]/area_lnd[nl][i],1.0);
	    l_a[i] = l_area[nl][i];
	    a_l[i] = area_lnd[nl][i];
	    a_a[i] = area_atm[nl][i];
	  }
	  if(ntile_lnd > 1)
	    sprintf(lnd_mask_file, "land_mask_tile%d.nc", nl+1);
	  else
	    strcpy(lnd_mask_file, "land_mask.nc");
	  fid = mpp_open(lnd_mask_file, MPP_WRITE);
    print_provenance_gv_gca(fid, history, grid_version, clip_method == GREAT_CIRCLE_CLIP );

	  dims[1] = mpp_def_dim(fid, "nx", nxl[nl]);
	  dims[0] = mpp_def_dim(fid, "ny", nyl[nl]);
	  id_mask = mpp_def_var(fid, "mask", MPP_DOUBLE, 2, dims,  2, "standard_name",
				"land fraction at T-cell centers", "units", "none");
	  id_a_a = mpp_def_var(fid, "area_atm", MPP_DOUBLE, 2, dims,  2, "standard_name",
				"area atm ", "units", "none");
	  id_a_l = mpp_def_var(fid, "area_lnd", MPP_DOUBLE, 2, dims,  2, "standard_name",
				"area land ", "units", "none");
	  id_l_a = mpp_def_var(fid, "l_area", MPP_DOUBLE, 2, dims,  2, "standard_name",
				"land x area", "units", "none");
	  mpp_end_def(fid);
	  mpp_put_var_value(fid, id_mask, mask);
	  mpp_put_var_value(fid, id_a_l, a_l);
	  mpp_put_var_value(fid, id_a_a, a_a);
	  mpp_put_var_value(fid, id_l_a, l_a);
	  free(mask);
	  mpp_close(fid);
	}
      }

      for(nl=0; nl<ntile_lnd; nl++) free(l_area[nl]);
      for(no=0; no<ntile_ocn; no++) free(o_area[no]);
      free(o_area);
      free(l_area);
    }


    for(na=0; na<ntile_atm; na++) {
    /* write out atmXlnd data*/
      for(nl = 0; nl < ntile_lnd; nl++) {
	int nxgrid;
	nxgrid = naxl[na][nl];
	mpp_sum_int(1, &nxgrid);
	if(nxgrid>0) {
	  size_t start[4], nwrite[4];
	  int *gdata_int;
	  double *gdata_dbl;

	  int fid, dim_string, dim_ncells, dim_two, dims[4];
	  int id_xgrid_area, id_contact, n;
	  int id_tile1_cell, id_tile2_cell, id_tile1_dist, id_tile2_dist;
	  char contact[STRING];

	  for(i=0; i<4; i++) {
	    start[i] = 0; nwrite[i] = 1;
	  }
	  if(same_mosaic)
	    sprintf(axl_file[nfile_axl], "atm_%s_%sXlnd_%s_%s.nc", amosaic_name, atile_name[na], lmosaic_name, ltile_name[nl]);
	  else
	    sprintf(axl_file[nfile_axl], "%s_%sX%s_%s.nc", amosaic_name, atile_name[na], lmosaic_name, ltile_name[nl]);
	  sprintf(contact, "%s:%s::%s:%s", amosaic_name, atile_name[na], lmosaic_name, ltile_name[nl]);

	  fid = mpp_open(axl_file[nfile_axl], MPP_WRITE);
    print_provenance_gv_gca(fid, history, grid_version,  clip_method == GREAT_CIRCLE_CLIP);

	  dim_string = mpp_def_dim(fid, "string", STRING);
	  dim_ncells = mpp_def_dim(fid, "ncells", nxgrid);
	  dim_two    = mpp_def_dim(fid, "two", 2);
	  if(interp_order == 2) {
	    id_contact = mpp_def_var(fid, "contact", MPP_CHAR, 1, &dim_string, 7, "standard_name", "grid_contact_spec",
				     "contact_type", "exchange", "parent1_cell",
				     "tile1_cell", "parent2_cell", "tile2_cell", "xgrid_area_field", "xgrid_area",
				     "distant_to_parent1_centroid", "tile1_distance", "distant_to_parent2_centroid", "tile2_distance");
	  }
	  else {
	    id_contact = mpp_def_var(fid, "contact", MPP_CHAR, 1, &dim_string, 5, "standard_name", "grid_contact_spec",
				     "contact_type", "exchange", "parent1_cell",
				     "tile1_cell", "parent2_cell", "tile2_cell", "xgrid_area_field", "xgrid_area");
	  }

	  dims[0] = dim_ncells; dims[1] = dim_two;
	  id_tile1_cell = mpp_def_var(fid, "tile1_cell", MPP_INT, 2, dims, 1, "standard_name", "parent_cell_indices_in_mosaic1");
	  id_tile2_cell = mpp_def_var(fid, "tile2_cell", MPP_INT, 2, dims, 1, "standard_name", "parent_cell_indices_in_mosaic2");
	  id_xgrid_area = mpp_def_var(fid, "xgrid_area", MPP_DOUBLE, 1, &dim_ncells, 2, "standard_name",
				      "exchange_grid_area", "units", "m2");
	  if(interp_order == 2) {
	    id_tile1_dist = mpp_def_var(fid, "tile1_distance", MPP_DOUBLE, 2, dims, 1, "standard_name", "distance_from_parent1_cell_centroid");
	    id_tile2_dist = mpp_def_var(fid, "tile2_distance", MPP_DOUBLE, 2, dims, 1, "standard_name", "distance_from_parent2_cell_centroid");
	  }
	  mpp_end_def(fid);

	  /* the index will start from 1, instead of 0 ( fortran index) */
	  for(i = 0;i < naxl[na][nl]; i++) {
	    ++(atmxlnd_ia[na][nl][i]);
	    ++(atmxlnd_ja[na][nl][i]);
	    ++(atmxlnd_il[na][nl][i]);
	    ++(atmxlnd_jl[na][nl][i]);
	  }
	  nwrite[0] = strlen(contact);
	  mpp_put_var_value_block(fid, id_contact, start, nwrite, contact);
     	  nwrite[0] = nxgrid;

	  gdata_int = (int *)malloc(nxgrid*sizeof(int));
	  gdata_dbl = (double *)malloc(nxgrid*sizeof(double));

	  mpp_gather_field_double(naxl[na][nl], atmxlnd_area[na][nl], gdata_dbl);
	  if(check) {
	    int *gdata_ia=NULL, *gdata_ja=NULL;
	    int ia, ja;

	    gdata_ia = (int *)malloc(nxgrid*sizeof(int));
	    gdata_ja = (int *)malloc(nxgrid*sizeof(int));
            mpp_gather_field_int(naxl[na][nl], atmxlnd_ia[na][nl], gdata_ia);
            mpp_gather_field_int(naxl[na][nl], atmxlnd_ja[na][nl], gdata_ja);
	    for(n=0; n<nxgrid; n++) {
	      ja = gdata_ja[n] - 1;
	      ia = gdata_ia[n] - 1;
	      atm_xarea[na][ja*nxa[na]+ia] += gdata_dbl[n];
	    }
	    free(gdata_ia);
	    free(gdata_ja);

	    if( na == tile_nest)
	      for(n=0; n<nxgrid; n++) axl_area_sum_nest += gdata_dbl[n];
	    else
	      for(n=0; n<nxgrid; n++) axl_area_sum += gdata_dbl[n];
	  }
	  mpp_put_var_value(fid, id_xgrid_area, gdata_dbl);
	  mpp_gather_field_int(naxl[na][nl], atmxlnd_ia[na][nl], gdata_int);
	  mpp_put_var_value_block(fid, id_tile1_cell, start, nwrite, gdata_int);
	  mpp_gather_field_int(naxl[na][nl], atmxlnd_il[na][nl], gdata_int);
	  mpp_put_var_value_block(fid, id_tile2_cell, start, nwrite, gdata_int);
	  start[1] = 1;
	  mpp_gather_field_int(naxl[na][nl], atmxlnd_ja[na][nl], gdata_int);
	  mpp_put_var_value_block(fid, id_tile1_cell, start, nwrite, gdata_int);
	  mpp_gather_field_int(naxl[na][nl], atmxlnd_jl[na][nl], gdata_int);
	  mpp_put_var_value_block(fid, id_tile2_cell, start, nwrite, gdata_int);
	  if(interp_order == 2) {
	    start[1] = 0;
  	    mpp_gather_field_double(naxl[na][nl], atmxlnd_dia[na][nl], gdata_dbl);
	    mpp_put_var_value_block(fid, id_tile1_dist, start, nwrite, gdata_dbl);
	    mpp_gather_field_double(naxl[na][nl], atmxlnd_dil[na][nl], gdata_dbl);
	    mpp_put_var_value_block(fid, id_tile2_dist, start, nwrite, gdata_dbl);
	    start[1] = 1;
  	    mpp_gather_field_double(naxl[na][nl], atmxlnd_dja[na][nl], gdata_dbl);
	    mpp_put_var_value_block(fid, id_tile1_dist, start, nwrite, gdata_dbl);
	    mpp_gather_field_double(naxl[na][nl], atmxlnd_djl[na][nl], gdata_dbl);
	    mpp_put_var_value_block(fid, id_tile2_dist, start, nwrite, gdata_dbl);
	  }
	  mpp_close(fid);
	  free(gdata_int);
	  free(gdata_dbl);
	  ++nfile_axl;
	}
      } /* end of nl loop */

      /* write out atmXocn data */
      for(no = 0; no < ntile_ocn; no++) {
	int nxgrid;

	nxgrid = naxo[na][no];
	mpp_sum_int(1, &nxgrid);
	if(nxgrid>0) {
	  size_t start[4], nwrite[4];
	  int *gdata_int;
	  double *gdata_dbl;
	  int fid, dim_string, dim_ncells, dim_two, dims[4];
	  int id_xgrid_area, id_contact, n;
	  int id_tile1_cell, id_tile2_cell, id_tile1_dist, id_tile2_dist;
	  char contact[STRING];

	  for(i=0; i<4; i++) {
	    start[i] = 0; nwrite[i] = 1;
	  }

	  if(same_mosaic)
	    sprintf(axo_file[nfile_axo], "atm_%s_%sXocn_%s_%s.nc", amosaic_name, atile_name[na], omosaic_name, otile_name[no]);
	  else
	    sprintf(axo_file[nfile_axo], "%s_%sX%s_%s.nc", amosaic_name, atile_name[na], omosaic_name, otile_name[no]);

	  sprintf(contact, "%s:%s::%s:%s", amosaic_name, atile_name[na], omosaic_name, otile_name[no]);
	  fid = mpp_open(axo_file[nfile_axo], MPP_WRITE);
    print_provenance_gv_gca(fid, history, grid_version, clip_method == GREAT_CIRCLE_CLIP );

	  dim_string = mpp_def_dim(fid, "string", STRING);
	  dim_ncells = mpp_def_dim(fid, "ncells", nxgrid);
	  dim_two    = mpp_def_dim(fid, "two", 2);
	  if(interp_order == 2) {
	    id_contact = mpp_def_var(fid, "contact", MPP_CHAR, 1, &dim_string, 7, "standard_name", "grid_contact_spec",
				   "contact_type", "exchange", "parent1_cell",
				   "tile1_cell", "parent2_cell", "tile2_cell", "xgrid_area_field", "xgrid_area",
				   "distant_to_parent1_centroid", "tile1_distance", "distant_to_parent2_centroid", "tile2_distance");
	  }
	  else {
	    id_contact = mpp_def_var(fid, "contact", MPP_CHAR, 1, &dim_string, 5, "standard_name", "grid_contact_spec",
				   "contact_type", "exchange", "parent1_cell",
				   "tile1_cell", "parent2_cell", "tile2_cell", "xgrid_area_field", "xgrid_area" );
	  }
	  dims[0] = dim_ncells; dims[1] = dim_two;
	  id_tile1_cell = mpp_def_var(fid, "tile1_cell", MPP_INT, 2, dims, 1, "standard_name", "parent_cell_indices_in_mosaic1");
	  id_tile2_cell = mpp_def_var(fid, "tile2_cell", MPP_INT, 2, dims, 1, "standard_name", "parent_cell_indices_in_mosaic2");
	  id_xgrid_area = mpp_def_var(fid, "xgrid_area", MPP_DOUBLE, 1, &dim_ncells, 2, "standard_name",
				      "exchange_grid_area", "units", "m2");
	  if(interp_order == 2) {
	    id_tile1_dist = mpp_def_var(fid, "tile1_distance", MPP_DOUBLE, 2, dims, 1, "standard_name",
					"distance_from_parent1_cell_centroid");
	    id_tile2_dist = mpp_def_var(fid, "tile2_distance", MPP_DOUBLE, 2, dims, 1, "standard_name",
					"distance_from_parent2_cell_centroid");
	  }
	  mpp_end_def(fid);

	  /* the index will start from 1, instead of 0 ( fortran index) */
	  for(i = 0;i < naxo[na][no]; i++) {
	    ++(atmxocn_ia[na][no][i]);
	    ++(atmxocn_ja[na][no][i]);
	    ++(atmxocn_io[na][no][i]);
	    atmxocn_jo[na][no][i] += 1-ocn_south_ext; /* possible one artificial j-level is added at south end */
	  }

          nwrite[0] = strlen(contact);
	  mpp_put_var_value_block(fid, id_contact, start, nwrite, contact);

	  nwrite[0] = nxgrid;

	  gdata_int = (int *)malloc(nxgrid*sizeof(int));
	  gdata_dbl = (double *)malloc(nxgrid*sizeof(double));

	  mpp_gather_field_double(naxo[na][no], atmxocn_area[na][no], gdata_dbl);

	  if(check) {
	    int *gdata_ia=NULL, *gdata_ja=NULL;
            int ia, ja;

	    gdata_ia = (int *)malloc(nxgrid*sizeof(int));
	    gdata_ja = (int *)malloc(nxgrid*sizeof(int));
            mpp_gather_field_int(naxo[na][no], atmxocn_ia[na][no], gdata_ia);
            mpp_gather_field_int(naxo[na][no], atmxocn_ja[na][no], gdata_ja);
	    for(n=0; n<nxgrid; n++) {
	      ja = gdata_ja[n] - 1;
	      ia = gdata_ia[n] - 1;
	      atm_xarea[na][ja*nxa[na]+ia] += gdata_dbl[n];
	    }
	    free(gdata_ia);
	    free(gdata_ja);

	    if( na == tile_nest)
	      for(n=0; n<nxgrid; n++) axo_area_sum_nest += gdata_dbl[n];
	    else
	      for(n=0; n<nxgrid; n++) axo_area_sum += gdata_dbl[n];
	  }
	  mpp_put_var_value_block(fid, id_xgrid_area, start, nwrite, gdata_dbl);
	  mpp_gather_field_int(naxo[na][no], atmxocn_ia[na][no], gdata_int);
	  mpp_put_var_value_block(fid, id_tile1_cell, start, nwrite, gdata_int);
	  mpp_gather_field_int(naxo[na][no], atmxocn_io[na][no], gdata_int);
	  mpp_put_var_value_block(fid, id_tile2_cell, start, nwrite, gdata_int);
	  start[1] = 1;
	  mpp_gather_field_int(naxo[na][no], atmxocn_ja[na][no], gdata_int);
	  mpp_put_var_value_block(fid, id_tile1_cell, start, nwrite, gdata_int);
	  mpp_gather_field_int(naxo[na][no], atmxocn_jo[na][no], gdata_int);
	  mpp_put_var_value_block(fid, id_tile2_cell, start, nwrite, gdata_int);
	  if(interp_order == 2) {
	    start[1] = 0;
  	    mpp_gather_field_double(naxo[na][no], atmxocn_dia[na][no], gdata_dbl);
	    mpp_put_var_value_block(fid, id_tile1_dist, start, nwrite, gdata_dbl);
	    mpp_gather_field_double(naxo[na][no], atmxocn_dio[na][no], gdata_dbl);
	    mpp_put_var_value_block(fid, id_tile2_dist, start, nwrite, gdata_dbl);
	    start[1] = 1;
  	    mpp_gather_field_double(naxo[na][no], atmxocn_dja[na][no], gdata_dbl);
	    mpp_put_var_value_block(fid, id_tile1_dist, start, nwrite, gdata_dbl);
	    mpp_gather_field_double(naxo[na][no], atmxocn_djo[na][no], gdata_dbl);
	    mpp_put_var_value_block(fid, id_tile2_dist, start, nwrite, gdata_dbl);
	  }
	  mpp_close(fid);
     	  free(gdata_int);
	  free(gdata_dbl);
	  ++nfile_axo;
	}
      } /* end of no loop */

    } /* end of na loop */

    /*release the memory */
    for(na=0; na<ntile_atm; na++) {
      for(nl=0; nl<ntile_lnd; nl++) {
	free(atmxlnd_area[na][nl]);
	free(atmxlnd_ia  [na][nl]);
	free(atmxlnd_ja  [na][nl]);
	free(atmxlnd_il  [na][nl]);
	free(atmxlnd_jl  [na][nl]);
	if(interp_order == 2) {
	  free(atmxlnd_clon[na][nl]);
	  free(atmxlnd_clat[na][nl]);
	  free(atmxlnd_dia [na][nl]);
	  free(atmxlnd_dja [na][nl]);
	  free(atmxlnd_dil [na][nl]);
	  free(atmxlnd_djl [na][nl]);
	}
      }
      free(atmxlnd_area[na]);
      free(atmxlnd_ia  [na]);
      free(atmxlnd_ja  [na]);
      free(atmxlnd_il  [na]);
      free(atmxlnd_jl  [na]);
      if(interp_order == 2) {
	free(atmxlnd_clon[na]);
	free(atmxlnd_clat[na]);
	free(atmxlnd_dia [na]);
	free(atmxlnd_dja [na]);
	free(atmxlnd_dil [na]);
	free(atmxlnd_djl [na]);
      }
      for(no=0; no<ntile_ocn; no++) {
	free(atmxocn_area[na][no]);
	free(atmxocn_ia  [na][no]);
	free(atmxocn_ja  [na][no]);
	free(atmxocn_io  [na][no]);
	free(atmxocn_jo  [na][no]);
	if(interp_order == 2) {
	  free(atmxocn_clon[na][no]);
	  free(atmxocn_clat[na][no]);
	  free(atmxocn_dia [na][no]);
	  free(atmxocn_dja [na][no]);
	  free(atmxocn_dio [na][no]);
	  free(atmxocn_djo [na][no]);
	}
      }
      free(atmxocn_area[na]);
      free(atmxocn_ia  [na]);
      free(atmxocn_ja  [na]);
      free(atmxocn_io  [na]);
      free(atmxocn_jo  [na]);
      if(interp_order == 2) {
	free(atmxocn_clon[na]);
	free(atmxocn_clat[na]);
	free(atmxocn_dia [na]);
	free(atmxocn_dja [na]);
	free(atmxocn_dio [na]);
	free(atmxocn_djo [na]);
      }
      free(naxl[na]);
      free(naxo[na]);
    }
    free(atmxlnd_area);
    free(atmxlnd_ia  );
    free(atmxlnd_ja  );
    free(atmxlnd_il  );
    free(atmxlnd_jl  );
    free(atmxocn_area);
    free(atmxocn_ia  );
    free(atmxocn_ja  );
    free(atmxocn_io  );
    free(atmxocn_jo  );
    if(interp_order == 2) {
      free(atmxlnd_clon);
      free(atmxlnd_clat);
      free(atmxlnd_dja );
      free(atmxlnd_dia );
      free(atmxlnd_dil );
      free(atmxlnd_djl );
      free(atmxocn_clon);
      free(atmxocn_clat);
      free(atmxocn_dia );
      free(atmxocn_dja );
      free(atmxocn_dio );
      free(atmxocn_djo );
    }
    free(naxl);
    free(naxo);
  }
  if(mpp_pe() == mpp_root_pe() && verbose) printf("\nNOTE from make_coupler_mosaic: Complete the process to create exchange grids "
				       "for fluxes between atmosphere and surface (sea ice and land)\n" );

  /***************************************************************************************
     Then generate the exchange grid between land mosaic and ocean mosaic
     if land mosaic is different from atmos mosaic
  ***************************************************************************************/
  nfile_lxo = 0;
  if( !lnd_same_as_atm ) {
    int     no, nl, ll, lo;
    size_t  **nlxo;
    int     ***lndxocn_il, ***lndxocn_jl, ***lndxocn_io, ***lndxocn_jo;
    double  ***lndxocn_area, ***lndxocn_dil, ***lndxocn_djl, ***lndxocn_dio, ***lndxocn_djo;
    double  ***lndxocn_clon, ***lndxocn_clat;
    double  min_area;

    nlxo         = (size_t ** )malloc(ntile_lnd*sizeof(size_t * ));
    lndxocn_area = (double ***)malloc(ntile_lnd*sizeof(double **));
    lndxocn_il   = (int    ***)malloc(ntile_lnd*sizeof(int    **));
    lndxocn_jl   = (int    ***)malloc(ntile_lnd*sizeof(int    **));
    lndxocn_io   = (int    ***)malloc(ntile_lnd*sizeof(int    **));
    lndxocn_jo   = (int    ***)malloc(ntile_lnd*sizeof(int    **));
    if(interp_order == 2) {
      lndxocn_dil  = (double ***)malloc(ntile_lnd*sizeof(double **));
      lndxocn_djl  = (double ***)malloc(ntile_lnd*sizeof(double **));
      lndxocn_dio  = (double ***)malloc(ntile_lnd*sizeof(double **));
      lndxocn_djo  = (double ***)malloc(ntile_lnd*sizeof(double **));
      lndxocn_clon = (double ***)malloc(ntile_lnd*sizeof(double **));
      lndxocn_clat = (double ***)malloc(ntile_lnd*sizeof(double **));
    }
    for(nl=0; nl<ntile_lnd; nl++) {
      nlxo        [nl] = (size_t * )malloc(ntile_ocn*sizeof(size_t  ));
      lndxocn_area[nl] = (double **)malloc(ntile_ocn*sizeof(double *));
      lndxocn_il  [nl] = (int    **)malloc(ntile_ocn*sizeof(int    *));
      lndxocn_jl  [nl] = (int    **)malloc(ntile_ocn*sizeof(int    *));
      lndxocn_io  [nl] = (int    **)malloc(ntile_ocn*sizeof(int    *));
      lndxocn_jo  [nl] = (int    **)malloc(ntile_ocn*sizeof(int    *));
      if(interp_order == 2) {
	lndxocn_dil [nl] = (double **)malloc(ntile_ocn*sizeof(double *));
	lndxocn_djl [nl] = (double **)malloc(ntile_ocn*sizeof(double *));
	lndxocn_dio [nl] = (double **)malloc(ntile_ocn*sizeof(double *));
	lndxocn_djo [nl] = (double **)malloc(ntile_ocn*sizeof(double *));
	lndxocn_clon[nl] = (double **)malloc(ntile_ocn*sizeof(double *));
	lndxocn_clat[nl] = (double **)malloc(ntile_ocn*sizeof(double *));
      }
      for(no=0; no<ntile_ocn; no++) {
	lndxocn_area[nl][no] = (double *)malloc(MAXXGRID*sizeof(double));
	lndxocn_il  [nl][no] = (int    *)malloc(MAXXGRID*sizeof(int   ));
	lndxocn_jl  [nl][no] = (int    *)malloc(MAXXGRID*sizeof(int   ));
	lndxocn_io  [nl][no] = (int    *)malloc(MAXXGRID*sizeof(int   ));
	lndxocn_jo  [nl][no] = (int    *)malloc(MAXXGRID*sizeof(int   ));
	if(interp_order == 2 ) {
	  lndxocn_dil[nl][no]  = (double *)malloc(MAXXGRID*sizeof(double));
	  lndxocn_djl[nl][no]  = (double *)malloc(MAXXGRID*sizeof(double));
	  lndxocn_dio[nl][no]  = (double *)malloc(MAXXGRID*sizeof(double));
	  lndxocn_djo[nl][no]  = (double *)malloc(MAXXGRID*sizeof(double));
	  lndxocn_clon[nl][no] = (double *)malloc(MAXXGRID*sizeof(double *));
	  lndxocn_clat[nl][no] = (double *)malloc(MAXXGRID*sizeof(double *));
	}
      }
    }

    for(nl=0; nl<ntile_lnd; nl++) {
      int      il, jl, io, jo, is, ie, js, je, layout[2];
      int      n0, n1, n2, n3, nl_in, no_in, n_out, nxgrid;
      double   xarea, xctrlon, xctrlat;
      double   xl_min, yl_min, xo_min, yo_min, xl_avg;
      double   xl_max, yl_max, xo_max, yo_max;
      double   xl[MV], yl[MV], zl[MV], xo[MV], yo[MV], zo[MV], x_out[MV], y_out[MV], z_out[MV];
      domain2D Dom;
      int      *is_ocn, *ie_ocn;
      int      *js_ocn, *je_ocn;
      double   yy;

      for(no=0; no<ntile_ocn; no++) nlxo[nl][no] = 0;

      layout[0] = mpp_npes();
      layout[1] = 1;

      mpp_define_domain2d(nxl[nl]*nyl[nl], 1, layout, 0, 0, &Dom);
      mpp_get_compute_domain2d(Dom, &is, &ie, &js, &je );

      js_ocn = (int *)malloc(ntile_ocn*sizeof(int));
      je_ocn = (int *)malloc(ntile_ocn*sizeof(int));
      is_ocn = (int *)malloc(ntile_ocn*sizeof(int));
      ie_ocn = (int *)malloc(ntile_ocn*sizeof(int));

      if(clip_method == GREAT_CIRCLE_CLIP) {
	for(no=0; no<ntile_ocn; no++) {
	  is_ocn[no] = 0;
	  ie_ocn[no] = nxo[no]-1;
	  js_ocn[no] = 0;
	  je_ocn[no] = nyo[no]-1;
	}
      }
      else {
	yl_min = 9999;
	yl_max = -9999;
	for(ll=is;ll<=ie;ll++) {

	  il = ll%nxl[nl];
	  jl = ll/nxl[nl];
	  n0 = jl    *(nxl[nl]+1) + il;
	  n1 = jl    *(nxl[nl]+1) + il+1;
	  n2 = (jl+1)*(nxl[nl]+1) + il+1;
	  n3 = (jl+1)*(nxl[nl]+1) + il;

	  yl[0] = ylnd[nl][n0];
	  yl[1] = ylnd[nl][n1];
	  yl[2] = ylnd[nl][n2];
	  yl[3] = ylnd[nl][n3];
	  if(yl[0] > yl_max) yl_max = yl[0];
	  if(yl[1] > yl_max) yl_max = yl[1];
	  if(yl[2] > yl_max) yl_max = yl[2];
	  if(yl[3] > yl_max) yl_max = yl[3];
	  if(yl[0] < yl_min) yl_min = yl[0];
	  if(yl[1] < yl_min) yl_min = yl[1];
	  if(yl[2] < yl_min) yl_min = yl[2];
	  if(yl[3] < yl_min) yl_min = yl[3];
	}
	for(no=0; no<ntile_ocn; no++) {
	  js_ocn[no] = nyo[no]; je_ocn[no] = -1;
	  for(jo = 0; jo <= nyo[no]; jo ++) for(io = 0; io <= nxo[no]; io++) {
	    yy = yocn[no][jo    *(nxo[no]+1) + io];
	    if( yy > yl_min && yy < yl_max ) {
	      if(jo > je_ocn[no] ) je_ocn[no] = jo;
	      if(jo < js_ocn[no] ) js_ocn[no] = jo;
	    }
	  }
	  js_ocn[no] = max(0, js_ocn[no]-1);
	  je_ocn[no] = min(nyo[no]-1, je_ocn[no]+1);
	  is_ocn[no] = 0;
	  ie_ocn[no] = nxo[no] - 1;
	}
      }

      for(ll=is;ll<=ie;ll++) {
	il = ll%nxl[nl];
	jl = ll/nxl[nl];
	if(clip_method == GREAT_CIRCLE_CLIP) { /* clockwise */
	  n0 = jl    *(nxl[nl]+1) + il;
	  n1 = (jl+1)*(nxl[nl]+1) + il;
	  n2 = (jl+1)*(nxl[nl]+1) + il+1;
	  n3 = jl    *(nxl[nl]+1) + il+1;
	  xl[0] = cart_xlnd[nl][n0]; yl[0] = cart_ylnd[nl][n0]; zl[0] = cart_zlnd[nl][n0];
	  xl[1] = cart_xlnd[nl][n1]; yl[1] = cart_ylnd[nl][n1]; zl[1] = cart_zlnd[nl][n1];
	  xl[2] = cart_xlnd[nl][n2]; yl[2] = cart_ylnd[nl][n2]; zl[2] = cart_zlnd[nl][n2];
	  xl[3] = cart_xlnd[nl][n3]; yl[3] = cart_ylnd[nl][n3]; zl[3] = cart_zlnd[nl][n3];
	}
	else {
	  n0 = jl    *(nxl[nl]+1) + il;
	  n1 = jl    *(nxl[nl]+1) + il+1;
	  n2 = (jl+1)*(nxl[nl]+1) + il+1;
	  n3 = (jl+1)*(nxl[nl]+1) + il;
	  xl[0] = xlnd[nl][n0]; yl[0] = ylnd[nl][n0];
	  xl[1] = xlnd[nl][n1]; yl[1] = ylnd[nl][n1];
	  xl[2] = xlnd[nl][n2]; yl[2] = ylnd[nl][n2];
	  xl[3] = xlnd[nl][n3]; yl[3] = ylnd[nl][n3];
	  yl_min  = minval_double(4, yl);
	  yl_max  = maxval_double(4, yl);
	  nl_in   = fix_lon(xl, yl, 4, M_PI);
	  xl_min  = minval_double(nl_in, xl);
	  xl_max  = maxval_double(nl_in, xl);
	  xl_avg  = avgval_double(nl_in, xl);
	}
	for(no=0; no<ntile_ocn; no++) {
	  for(jo = js_ocn[no]; jo <= je_ocn[no]; jo++) for(io = is_ocn[no]; io <= ie_ocn[no]; io++) if(omask[no][jo*nxo[no]+io] > MIN_AREA_FRAC) {
	    double ocn_frac;
	    if(clip_method == GREAT_CIRCLE_CLIP) { /* clockwise */
	      n0 = jo    *(nxo[no]+1) + io;
	      n1 = (jo+1)*(nxo[no]+1) + io;
	      n2 = (jo+1)*(nxo[no]+1) + io+1;
	      n3 = jo    *(nxo[no]+1) + io+1;
	      xo[0] = cart_xocn[no][n0]; yo[0] = cart_yocn[no][n0]; zo[0] = cart_zocn[no][n0];
	      xo[1] = cart_xocn[no][n1]; yo[1] = cart_yocn[no][n1]; zo[1] = cart_zocn[no][n1];
	      xo[2] = cart_xocn[no][n2]; yo[2] = cart_yocn[no][n2]; zo[2] = cart_zocn[no][n2];
	      xo[3] = cart_xocn[no][n3]; yo[3] = cart_yocn[no][n3]; zo[3] = cart_zocn[no][n3];
	    }
	    else {
	      n0 = jo    *(nxo[no]+1) + io;
	      n1 = jo    *(nxo[no]+1) + io+1;
	      n2 = (jo+1)*(nxo[no]+1) + io+1;
	      n3 = (jo+1)*(nxo[no]+1) + io;
	      xo[0] = xocn[no][n0]; yo[0] = yocn[no][n0];
	      xo[1] = xocn[no][n1]; yo[1] = yocn[no][n1];
	      xo[2] = xocn[no][n2]; yo[2] = yocn[no][n2];
	      xo[3] = xocn[no][n3]; yo[3] = yocn[no][n3];
	      yo_min = minval_double(4, yo);
	      yo_max = maxval_double(4, yo);
	      if(yo_min >= yl_max || yo_max <= yl_min ) continue;
	      no_in  = fix_lon(xo, yo, 4, xl_avg);
	      xo_min = minval_double(no_in, xo);
	      xo_max = maxval_double(no_in, xo);
	      /* xo should in the same range as xa after lon_fix, so no need to
		 consider cyclic condition
	      */
	      if(xl_min >= xo_max || xl_max <= xo_min ) continue;
	    }
	    ocn_frac = omask[no][jo*nxo[no]+io];
	    if(clip_method == GREAT_CIRCLE_CLIP) {
	      n_out = clip_2dx2d_great_circle(xl, yl, zl, 4, xo, yo, zo, 4,
						x_out, y_out, z_out);
	    }
	    else
	      n_out = clip_2dx2d( xl, yl, nl_in, xo, yo, no_in, x_out, y_out );

	    if (  n_out > 0 ){
	      if(clip_method == GREAT_CIRCLE_CLIP)
		xarea=great_circle_area ( n_out, x_out, y_out, z_out)*ocn_frac;
	      else
		xarea = poly_area(x_out, y_out, n_out )*ocn_frac;
	      min_area = min(area_ocn[no][jo*nxo[no]+io], area_lnd[nl][ll] );
	      if(xarea/min_area > area_ratio_thresh ) {
		lndxocn_area[nl][no][nlxo[nl][no]] = xarea;
		lndxocn_io[nl][no][nlxo[nl][no]]   = io;
		lndxocn_jo[nl][no][nlxo[nl][no]]   = jo;
		lndxocn_il[nl][no][nlxo[nl][no]]   = il;
		lndxocn_jl[nl][no][nlxo[nl][no]]   = jl;
		if(interp_order == 2) {
		  lndxocn_clon[nl][no][nlxo[nl][no]] = poly_ctrlon ( x_out, y_out, n_out, xl_avg)*ocn_frac;
		  lndxocn_clat[nl][no][nlxo[nl][no]] = poly_ctrlat ( x_out, y_out, n_out )*ocn_frac;
		}
		++(nlxo[nl][no]);
		if(nlxo[nl][no] > MAXXGRID) mpp_error("nlxo is greater than MAXXGRID, increase MAXXGRID");
	      }
	    }
	  } /* end of io, jo loop */
	}
      } /* end of ll( or il, jl) loop */
      mpp_delete_domain2d(&Dom);
      free(js_ocn);
      free(je_ocn);
      free(is_ocn);
      free(ie_ocn);
    }/* for(nl=0; nl<ntile_lnd; nl++) */

    /* calculate the centroid of model grid. */
    if(interp_order == 2) {
      double *l_area, *l_clon, *l_clat;
      double **o_area, **o_clon, **o_clat;

      o_area = (double **)malloc(ntile_ocn*sizeof(double *));
      o_clon = (double **)malloc(ntile_ocn*sizeof(double *));
      o_clat = (double **)malloc(ntile_ocn*sizeof(double *));
      for(no =0; no<ntile_ocn; no++) {
	o_area[no] = (double *)malloc(nxo[no]*nyo[no]*sizeof(double));
	o_clon[no] = (double *)malloc(nxo[no]*nyo[no]*sizeof(double));
	o_clat[no] = (double *)malloc(nxo[no]*nyo[no]*sizeof(double));
	for(lo=0; lo<nxo[no]*nyo[no]; lo++) {
	  o_area[no][lo] = 0;
	  o_clon[no][lo] = 0;
	  o_clat[no][lo] = 0;
	}
      }

      for(nl=0; nl<ntile_lnd; nl++) {
	l_area = (double *)malloc(nxl[nl]*nyl[nl]*sizeof(double));
	l_clon = (double *)malloc(nxl[nl]*nyl[nl]*sizeof(double));
	l_clat = (double *)malloc(nxl[nl]*nyl[nl]*sizeof(double));
	for(ll=0; ll<nxl[nl]*nyl[nl]; ll++) {
	  l_area[ll] = 0;
	  l_clon[ll] = 0;
	  l_clat[ll] = 0;
	}

	for(no=0; no<ntile_ocn; no++) {
	  int nxgrid;
	  nxgrid = nlxo[nl][no];
	  mpp_sum_int(1, &nxgrid);
	  if(nxgrid > 0) {
	    double *g_area, *g_clon, *g_clat;
	    int    *g_il,   *g_jl,   *g_io, *g_jo;
	    int    ii;
	    g_il = (int    *)malloc(nxgrid*sizeof(int   ));
	    g_jl = (int    *)malloc(nxgrid*sizeof(int   ));
	    g_io = (int    *)malloc(nxgrid*sizeof(int   ));
	    g_jo = (int    *)malloc(nxgrid*sizeof(int   ));
	    g_area = (double *)malloc(nxgrid*sizeof(double));
	    g_clon = (double *)malloc(nxgrid*sizeof(double));
	    g_clat = (double *)malloc(nxgrid*sizeof(double));
	    mpp_gather_field_int   (nlxo[nl][no], lndxocn_il[nl][no], g_il);
	    mpp_gather_field_int   (nlxo[nl][no], lndxocn_jl[nl][no], g_jl);
	    mpp_gather_field_int   (nlxo[nl][no], lndxocn_io[nl][no], g_io);
	    mpp_gather_field_int   (nlxo[nl][no], lndxocn_jo[nl][no], g_jo);
	    mpp_gather_field_double(nlxo[nl][no], lndxocn_area[nl][no], g_area);
	    mpp_gather_field_double(nlxo[nl][no], lndxocn_clon[nl][no], g_clon);
	    mpp_gather_field_double(nlxo[nl][no], lndxocn_clat[nl][no], g_clat);
	    for(i=0; i<nxgrid; i++) {
	      ii = g_jl[i]*nxl[nl]+g_il[i];
	      l_area[ii] += g_area[i];
	      l_clon[ii] += g_clon[i];
	      l_clat[ii] += g_clat[i];
	      ii = g_jo[i]*nxo[no]+g_io[i];
	      o_area[no][ii] += g_area[i];
	      o_clon[no][ii] += g_clon[i];
	      o_clat[no][ii] += g_clat[i];
	    }
	    free(g_il);
	    free(g_jl);
	    free(g_io);
	    free(g_jo);
	    free(g_area);
	    free(g_clon);
	    free(g_clat);
	  }
	}
	for(ll=0; ll<nxl[nl]*nyl[nl]; ll++) {
	  if(l_area[ll] > 0) {
	    l_clon[ll] /= l_area[ll];
	    l_clat[ll] /= l_area[ll];
	  }
	}
	/* substract land centroid to get the centroid distance between land grid and exchange grid. */
	for(no=0; no<ntile_ocn; no++) {
	  for(i=0; i<nlxo[nl][no]; i++) {
	    ll = lndxocn_jl[nl][no][i]*nxl[nl] + lndxocn_il[nl][no][i];
	    lndxocn_dil[nl][no][i] = lndxocn_clon[nl][no][i]/lndxocn_area[nl][no][i] - l_clon[ll];
	    lndxocn_djl[nl][no][i] = lndxocn_clat[nl][no][i]/lndxocn_area[nl][no][i] - l_clat[ll];
	  }
	}

	free(l_area);
	free(l_clon);
	free(l_clat);
      }

      /* centroid distance from exchange grid to ocean grid */
      for(no=0; no<ntile_ocn; no++) {
	for(lo=0; lo<nxo[no]*nyo[no]; lo++) {
	  if(o_area[no][lo] > 0) {
	    o_clon[no][lo] /= o_area[no][lo];
	    o_clat[no][lo] /= o_area[no][lo];
	  }
	}
	for(nl=0; nl<ntile_lnd; nl++) {
	  for(i=0; i<nlxo[nl][no]; i++) {
	    lo = lndxocn_jo[nl][no][i]*nxo[no] + lndxocn_io[nl][no][i];
	    lndxocn_dio[nl][no][i] = lndxocn_clon[nl][no][i]/lndxocn_area[nl][no][i] - o_clon[no][lo];
	    lndxocn_djo[nl][no][i] = lndxocn_clat[nl][no][i]/lndxocn_area[nl][no][i] - o_clat[no][lo];
	  }
	}
	free(o_area[no]);
	free(o_clon[no]);
	free(o_clat[no]);
      }

      free(o_area);
      free(o_clon);
      free(o_clat);
    }

    /* write out lndXocn data */
    for(nl = 0; nl < ntile_lnd; nl++) {
       for(no = 0; no < ntile_ocn; no++) {
	 int nxgrid;

	/* get total number of exchange grid on all the pes */
	 nxgrid = nlxo[nl][no];
	mpp_sum_int(1, &nxgrid);

	if(nxgrid >0) {
	  size_t start[4], nwrite[4];
	  int *gdata_int;
	  double *gdata_dbl;

	  char contact[STRING];
	  int fid, dim_string, dim_ncells, dim_two, dims[4];
	  int id_contact, id_xgrid_area, n;
	  int id_tile1_cell, id_tile2_cell, id_tile1_dist, id_tile2_dist;

	  for(i=0; i<4; i++) {
	    start[i] = 0; nwrite[i] = 1;
	  }

	  if(same_mosaic)
	    sprintf(lxo_file[nfile_lxo], "lnd_%s_%sXocn_%s_%s.nc", lmosaic_name, ltile_name[nl], omosaic_name, otile_name[no]);
	  else
	    sprintf(lxo_file[nfile_lxo], "%s_%sX%s_%s.nc", lmosaic_name, ltile_name[nl], omosaic_name, otile_name[no]);
	  sprintf(contact, "%s:%s::%s:%s", lmosaic_name, ltile_name[nl], omosaic_name, otile_name[no]);

	  fid = mpp_open(lxo_file[nfile_lxo], MPP_WRITE);
    print_provenance_gv_gca(fid, history, grid_version, clip_method == GREAT_CIRCLE_CLIP );

	  dim_string = mpp_def_dim(fid, "string", STRING);
	  dim_ncells = mpp_def_dim(fid, "ncells", nxgrid);
	  dim_two    = mpp_def_dim(fid, "two", 2);
	  if(interp_order == 2) {
	    id_contact = mpp_def_var(fid, "contact", MPP_CHAR, 1, &dim_string, 7, "standard_name", "grid_contact_spec",
				     "contact_type", "exchange", "parent1_cell",
				     "tile1_cell", "parent2_cell", "tile2_cell", "xgrid_area_field", "xgrid_area",
				     "distant_to_parent1_centroid", "tile1_distance", "distant_to_parent2_centroid", "tile2_distance");
	  }
	  else {
	    id_contact = mpp_def_var(fid, "contact", MPP_CHAR, 1, &dim_string, 5, "standard_name", "grid_contact_spec",
				     "contact_type", "exchange", "parent1_cell",
				     "tile1_cell", "parent2_cell", "tile2_cell", "xgrid_area_field", "xgrid_area" );
	  }
	  dims[0] = dim_ncells; dims[1] = dim_two;
	  id_tile1_cell = mpp_def_var(fid, "tile1_cell", MPP_INT, 2, dims, 1, "standard_name", "parent1_cell_indices");
	  id_tile2_cell = mpp_def_var(fid, "tile2_cell", MPP_INT, 2, dims, 1, "standard_name", "parent2_cell_indices");
	  id_xgrid_area = mpp_def_var(fid, "xgrid_area", MPP_DOUBLE, 1, &dim_ncells, 2, "standard_name",
				      "exchange_grid_area", "units", "m2");

	  if(interp_order == 2) {
	    id_tile1_dist = mpp_def_var(fid, "tile1_distance", MPP_DOUBLE, 2, dims, 1, "standard_name", "distance_from_parent1_cell_centroid");
	    id_tile2_dist = mpp_def_var(fid, "tile2_distance", MPP_DOUBLE, 2, dims, 1, "standard_name", "distance_from_parent2_cell_centroid");
	  }
	  mpp_end_def(fid);

	  /* the index will start from 1, instead of 0 ( fortran index) */
	  for(i = 0;i < nlxo[nl][no]; i++) {
	    ++(lndxocn_il[nl][no][i]);
	    ++(lndxocn_jl[nl][no][i]);
	    ++(lndxocn_io[nl][no][i]);
	    lndxocn_jo[nl][no][i] += 1 - ocn_south_ext; /* one artificial j-level may be added in the south end */
	  }

	  nwrite[0] = strlen(contact);
	  mpp_put_var_value_block(fid, id_contact, start, nwrite, contact);
	  nwrite[0] = nxgrid;

	  gdata_int = (int *)malloc(nxgrid*sizeof(int));
	  gdata_dbl = (double *)malloc(nxgrid*sizeof(double));

	  mpp_gather_field_double(nlxo[nl][no], lndxocn_area[nl][no], gdata_dbl);

	  mpp_put_var_value(fid, id_xgrid_area, gdata_dbl);
	  mpp_gather_field_int(nlxo[nl][no], lndxocn_il[nl][no], gdata_int);
	  mpp_put_var_value_block(fid, id_tile1_cell, start, nwrite, gdata_int);
	  mpp_gather_field_int(nlxo[nl][no], lndxocn_io[nl][no], gdata_int);
	  mpp_put_var_value_block(fid, id_tile2_cell, start, nwrite, gdata_int);
	  start[1] = 1;
	  mpp_gather_field_int(nlxo[nl][no], lndxocn_jl[nl][no], gdata_int);
	  mpp_put_var_value_block(fid, id_tile1_cell, start, nwrite, gdata_int);
	  mpp_gather_field_int(nlxo[nl][no], lndxocn_jo[nl][no], gdata_int);
	  mpp_put_var_value_block(fid, id_tile2_cell, start, nwrite, gdata_int);
	  if(interp_order == 2) {
	    start[1] = 0;
  	    mpp_gather_field_double(nlxo[nl][no], lndxocn_dil[nl][no], gdata_dbl);
	    mpp_put_var_value_block(fid, id_tile1_dist, start, nwrite, gdata_dbl);
	    mpp_gather_field_double(nlxo[nl][no], lndxocn_dio[nl][no], gdata_dbl);
	    mpp_put_var_value_block(fid, id_tile2_dist, start, nwrite, gdata_dbl);
	    start[1] = 1;
  	    mpp_gather_field_double(nlxo[nl][no], lndxocn_djl[nl][no], gdata_dbl);
	    mpp_put_var_value_block(fid, id_tile1_dist, start, nwrite, gdata_dbl);
	    mpp_gather_field_double(nlxo[nl][no], lndxocn_djo[nl][no], gdata_dbl);
	    mpp_put_var_value_block(fid, id_tile2_dist, start, nwrite, gdata_dbl);
	  }
	  mpp_close(fid);
	  free(gdata_int);
	  free(gdata_dbl);
	  ++nfile_lxo;
	}
       } /* for(no=0; no<ntile_ocn; no++) */
    } /* for(nl=0; nl<ntile_lnd; nl++) */

    /* release the memory */
    for(nl=0; nl<ntile_lnd; nl++) {
      for(no=0; no<ntile_ocn; no++) {
	free(lndxocn_area[nl][no]);
	free(lndxocn_il  [nl][no]);
	free(lndxocn_jl  [nl][no]);
	free(lndxocn_io  [nl][no]);
	free(lndxocn_jo  [nl][no]);
	if(interp_order == 2) {
	  free(lndxocn_clon[nl][no]);
	  free(lndxocn_clat[nl][no]);
	  free(lndxocn_dil [nl][no]);
	  free(lndxocn_djl [nl][no]);
	  free(lndxocn_dio [nl][no]);
	  free(lndxocn_djo [nl][no]);
	}
      }
      free(lndxocn_area[nl]);
      free(lndxocn_il  [nl]);
      free(lndxocn_jl  [nl]);
      free(lndxocn_io  [nl]);
      free(lndxocn_jo  [nl]);
      if(interp_order == 2) {
	free(lndxocn_clon[nl]);
	free(lndxocn_clat[nl]);
	free(lndxocn_dil [nl]);
	free(lndxocn_djl [nl]);
	free(lndxocn_dio [nl]);
	free(lndxocn_djo [nl]);
      }
      free(nlxo[nl]);
    }
    free(nlxo);
    free(lndxocn_area);
    free(lndxocn_il  );
    free(lndxocn_jl  );
    free(lndxocn_io  );
    free(lndxocn_jo  );
    if(interp_order == 2) {
      free(lndxocn_clon);
      free(lndxocn_clat);
      free(lndxocn_dil );
      free(lndxocn_djl );
      free(lndxocn_dio );
      free(lndxocn_djo );
    }

    if(mpp_pe() == mpp_root_pe() && verbose) printf("\nNOTE from make_coupler_mosaic: Complete the process to create exchange grids "
					 "for runoff between land and sea ice.\n" );
  }
  else {
    nfile_lxo = nfile_axo;
    for(i=0; i<nfile_axo; i++) strcpy(lxo_file[i], axo_file[i]);
    if(mpp_pe() == mpp_root_pe() && verbose) printf("\nNOTE from make_coupler_mosaic: Since lmosaic is the same as amosaic, "
					 "no need to compute the exchange grid between lmosaic and omosaic, "
					 "simply use the exchange grid between amosaic and omosaic.\n");
  }

  /***************************************************************************************
     Now generate the exchange grid between wave mosaic and ocean mosaic
     if wave mosaic exist.
  ***************************************************************************************/

  if(wmosaic) {
     nfile_wxo = 0;
     int no, nw, n;
     size_t  **nwxo;

    int     ***wavxocn_iw,   ***wavxocn_jw,   ***wavxocn_io,   ***wavxocn_jo;
    double  ***wavxocn_area, ***wavxocn_diw,  ***wavxocn_djw,  ***wavxocn_dio,  ***wavxocn_djo;
    double  ***wavxocn_clon, ***wavxocn_clat;
    double   min_area;
    time_t time_start, time_end;

    nwxo         = (size_t ** )malloc(ntile_wav*sizeof(size_t *));
    wavxocn_area = (double ***)malloc(ntile_wav*sizeof(double **));
    wavxocn_iw   = (int    ***)malloc(ntile_wav*sizeof(int    **));
    wavxocn_jw   = (int    ***)malloc(ntile_wav*sizeof(int    **));
    wavxocn_io   = (int    ***)malloc(ntile_wav*sizeof(int    **));
    wavxocn_jo   = (int    ***)malloc(ntile_wav*sizeof(int    **));

    if(interp_order == 2 ) {
      wavxocn_diw  = (double ***)malloc(ntile_wav*sizeof(double **));
      wavxocn_djw  = (double ***)malloc(ntile_wav*sizeof(double **));
      wavxocn_dio  = (double ***)malloc(ntile_wav*sizeof(double **));
      wavxocn_djo  = (double ***)malloc(ntile_wav*sizeof(double **));
      wavxocn_clon = (double ***)malloc(ntile_wav*sizeof(double **));
      wavxocn_clat = (double ***)malloc(ntile_wav*sizeof(double **));
    }

    for(nw=0; nw<ntile_wav; nw++) {
      nwxo[nw]         = (size_t * )malloc(ntile_ocn*sizeof(size_t));
      wavxocn_area[nw] = (double **)malloc(ntile_ocn*sizeof(double *));
      wavxocn_iw[nw]   = (int    **)malloc(ntile_ocn*sizeof(int    *));
      wavxocn_jw[nw]   = (int    **)malloc(ntile_ocn*sizeof(int    *));
      wavxocn_io[nw]   = (int    **)malloc(ntile_ocn*sizeof(int    *));
      wavxocn_jo[nw]   = (int    **)malloc(ntile_ocn*sizeof(int    *));

      if(interp_order == 2 ) {
	wavxocn_diw [nw] = (double **)malloc(ntile_ocn*sizeof(double *));
	wavxocn_djw [nw] = (double **)malloc(ntile_ocn*sizeof(double *));
	wavxocn_dio [nw] = (double **)malloc(ntile_ocn*sizeof(double *));
	wavxocn_djo [nw] = (double **)malloc(ntile_ocn*sizeof(double *));
	wavxocn_clon[nw] = (double **)malloc(ntile_ocn*sizeof(double *));
	wavxocn_clat[nw] = (double **)malloc(ntile_ocn*sizeof(double *));
      }

      for(no=0; no<ntile_ocn; no++) {
	wavxocn_area[nw][no] = (double *)malloc(MAXXGRID*sizeof(double));
	wavxocn_iw  [nw][no] = (int    *)malloc(MAXXGRID*sizeof(int   ));
	wavxocn_jw  [nw][no] = (int    *)malloc(MAXXGRID*sizeof(int   ));
	wavxocn_io  [nw][no] = (int    *)malloc(MAXXGRID*sizeof(int   ));
	wavxocn_jo  [nw][no] = (int    *)malloc(MAXXGRID*sizeof(int   ));
	if(interp_order == 2 ) {
	  wavxocn_clon[nw][no] = (double *)malloc(MAXXGRID*sizeof(double));
	  wavxocn_clat[nw][no] = (double *)malloc(MAXXGRID*sizeof(double));
	  wavxocn_diw [nw][no] = (double *)malloc(MAXXGRID*sizeof(double));
	  wavxocn_djw [nw][no] = (double *)malloc(MAXXGRID*sizeof(double));
	  wavxocn_dio [nw][no] = (double *)malloc(MAXXGRID*sizeof(double));
	  wavxocn_djo [nw][no] = (double *)malloc(MAXXGRID*sizeof(double));
	}
      }
    }

    time_start = time(NULL);

    for(nw=0; nw<ntile_wav; nw++) {

      int      l, is, ie, js, je, lw, iw, jw, ia, ja, io, jo, layout[2];
      int      n0, n1, n2, n3, na_in, nw_in, no_in, n_out, n_out2;
      double   xo_min, yo_min, xw_min, yw_min, xw_avg;
      double   xo_max, yo_max, xw_max, yw_max;
      double   xarea;
      double   xw[MV], yw[MV], zw[MV], xo[MV], yo[MV], zo[MV];
      double   x_out[MV], y_out[MV], z_out[MV];
      size_t   count;
      domain2D Dom;
      double   yy;
      int      *js_ocn, *je_ocn;
      int      *is_ocn, *ie_ocn;

      for(no=0; no<ntile_ocn; no++) nwxo[nw][no] = 0;
      layout[0] = mpp_npes();
      layout[1] = 1;

      mpp_define_domain2d(nxw[nw]*nyw[nw], 1, layout, 0, 0, &Dom);
      mpp_get_compute_domain2d(Dom, &is, &ie, &js, &je );
      /* find the js_ocn, je_ocn, and is_ocn, ie_ocn
         In x-direction, cyclic condition will be considered */
      is_ocn = (int *)malloc(ntile_ocn*sizeof(int));
      ie_ocn = (int *)malloc(ntile_ocn*sizeof(int));
      js_ocn = (int *)malloc(ntile_ocn*sizeof(int));
      je_ocn = (int *)malloc(ntile_ocn*sizeof(int));
      yw_min = 9999;
      yw_max = -9999;

      if(clip_method == GREAT_CIRCLE_CLIP) {
	for(no=0; no<ntile_ocn; no++) {
	  is_ocn[no] = 0;
	  ie_ocn[no] = nxo[no]-1;
	  js_ocn[no] = 0;
	  je_ocn[no] = nyo[no]-1;
	}
      }
      else {
	for(lw=is;lw<=ie;lw++) {

	  iw = lw%nxw[nw];
	  jw = lw/nxw[nw];
	  n0 = jw    *(nxw[nw]+1) + iw;
	  n1 = jw    *(nxw[nw]+1) + iw+1;
	  n2 = (jw+1)*(nxw[nw]+1) + iw+1;
	  n3 = (jw+1)*(nxw[nw]+1) + iw;

	  yw[0] = ywav[nw][n0];
	  yw[1] = ywav[nw][n1];
	  yw[2] = ywav[nw][n2];
	  yw[3] = ywav[nw][n3];
	  if(yw[0] > yw_max) yw_max = yw[0];
	  if(yw[1] > yw_max) yw_max = yw[1];
	  if(yw[2] > yw_max) yw_max = yw[2];
	  if(yw[3] > yw_max) yw_max = yw[3];
	  if(yw[0] < yw_min) yw_min = yw[0];
	  if(yw[1] < yw_min) yw_min = yw[1];
	  if(yw[2] < yw_min) yw_min = yw[2];
	  if(yw[3] < yw_min) yw_min = yw[3];
	}
	/*      printf("yw_min=%f, yw_max=%f\n", yw_min, yw_max); */

	for(no=0; no<ntile_ocn; no++) {
	  js_ocn[no] = nyo[no]; je_ocn[no] = -1;
	  for(jo = 0; jo <= nyo[no]; jo ++) for(io = 0; io <= nxo[no]; io++) {
	    yy = yocn[no][jo    *(nxo[no]+1) + io];
	    if( yy > yw_min && yy < yw_max ) {
	      if(jo > je_ocn[no] ) je_ocn[no] = jo;
	      if(jo < js_ocn[no] ) js_ocn[no] = jo;
	    }
	  }
	  js_ocn[no] = max(0, js_ocn[no]-1);
	  je_ocn[no] = min(nyo[no]-1, je_ocn[no]+1);

	  if(no==nw ||  !wav_same_as_ocn ) {
	    is_ocn[no] = 0;
	    ie_ocn[no] = nxo[no]-1;
	  }
	  else {
	    is_ocn[no] = nxo[no]-1;
	    ie_ocn[no] = 0;
	  }
	}
      }
      for(lw=is;lw<=ie;lw++) {
       	if(mpp_pe()==mpp_root_pe() && verbose)printf("nw = %d, lw = %d, is=%d, ie = %d\n", nw, lw, is, ie);
	iw = lw%nxw[nw];
	jw = lw/nxw[nw];
	if(clip_method == GREAT_CIRCLE_CLIP) { /* clockwise */
	  n0 = jw    *(nxw[nw]+1) + iw;
	  n1 = (jw+1)*(nxw[nw]+1) + iw;
	  n2 = (jw+1)*(nxw[nw]+1) + iw+1;
	  n3 = jw    *(nxw[nw]+1) + iw+1;
	  xw[0] = cart_xwav[nw][n0]; yw[0] = cart_ywav[nw][n0]; zw[0] = cart_zwav[nw][n0];
	  xw[1] = cart_xwav[nw][n1]; yw[1] = cart_ywav[nw][n1]; zw[1] = cart_zwav[nw][n1];
	  xw[2] = cart_xwav[nw][n2]; yw[2] = cart_ywav[nw][n2]; zw[2] = cart_zwav[nw][n2];
	  xw[3] = cart_xwav[nw][n3]; yw[3] = cart_ywav[nw][n3]; zw[3] = cart_zwav[nw][n3];
	}
	else {
	  n0 = jw    *(nxw[nw]+1) + iw;
	  n1 = jw    *(nxw[nw]+1) + iw+1;
	  n2 = (jw+1)*(nxw[nw]+1) + iw+1;
	  n3 = (jw+1)*(nxw[nw]+1) + iw;
	  xw[0] = xwav[nw][n0]; yw[0] = ywav[nw][n0];
	  xw[1] = xwav[nw][n1]; yw[1] = ywav[nw][n1];
	  xw[2] = xwav[nw][n2]; yw[2] = ywav[nw][n2];
	  xw[3] = xwav[nw][n3]; yw[3] = ywav[nw][n3];
	  yw_min  = minval_double(4, yw);
	  yw_max  = maxval_double(4, yw);
	  nw_in   = fix_lon(xw, yw, 4, M_PI);
	  xw_min  = minval_double(nw_in, xw);
	  xw_max  = maxval_double(nw_in, xw);
	  xw_avg  = avgval_double(nw_in, xw);
	}

	/* calculate wave/ocean x-cells */
	for(no=0; no<ntile_ocn; no++) {
	  for(jo = js_ocn[no]; jo <= je_ocn[no]; jo++) for(io = is_ocn[no]; io <= ie_ocn[no]; io++) {
	    double ocn_frac;
	    if(clip_method == GREAT_CIRCLE_CLIP) { /* clockwise */
	      n0 = jo    *(nxo[no]+1) + io;
	      n1 = (jo+1)*(nxo[no]+1) + io;
	      n2 = (jo+1)*(nxo[no]+1) + io+1;
	      n3 = jo    *(nxo[no]+1) + io+1;
	      xo[0] = cart_xocn[no][n0]; yo[0] = cart_yocn[no][n0]; zo[0] = cart_zocn[no][n0];
	      xo[1] = cart_xocn[no][n1]; yo[1] = cart_yocn[no][n1]; zo[1] = cart_zocn[no][n1];
	      xo[2] = cart_xocn[no][n2]; yo[2] = cart_yocn[no][n2]; zo[2] = cart_zocn[no][n2];
	      xo[3] = cart_xocn[no][n3]; yo[3] = cart_yocn[no][n3]; zo[3] = cart_zocn[no][n3];
	    }
	    else {
	      n0 = jo    *(nxo[no]+1) + io;
	      n1 = jo    *(nxo[no]+1) + io+1;
	      n2 = (jo+1)*(nxo[no]+1) + io+1;
	      n3 = (jo+1)*(nxo[no]+1) + io;
	      xo[0] = xocn[no][n0]; yo[0] = yocn[no][n0];
	      xo[1] = xocn[no][n1]; yo[1] = yocn[no][n1];
	      xo[2] = xocn[no][n2]; yo[2] = yocn[no][n2];
	      xo[3] = xocn[no][n3]; yo[3] = yocn[no][n3];
	      yo_min = minval_double(4, yo);
	      yo_max = maxval_double(4, yo);
	      no_in  = fix_lon(xo, yo, 4, xw_avg);
	      xo_min = minval_double(no_in, xo);
	      xo_max = maxval_double(no_in, xo);
	    }
	    ocn_frac = omask[no][jo*nxo[no]+io];
	    if(ocn_frac > MIN_AREA_FRAC) { /* over sea/ice */
	      /* xo should in the same range as xa after lon_fix, so no need to
		 consider cyclic condition
	      */
	      if(clip_method == GREAT_CIRCLE_CLIP) {
		n_out = clip_2dx2d_great_circle(xw, yw, zw, 4, xo, yo, zo, 4,
						x_out, y_out, z_out);
	      }
	      else {
		if(xw_min >= xo_max || xw_max <= xo_min || yw_min >= yo_max || yw_max <= yo_min ) continue;
		n_out = clip_2dx2d( xw, yw, nw_in, xo, yo, no_in, x_out, y_out);
	      }
	      if (  n_out > 0) {
		xarea = poly_area(x_out, y_out, n_out )*ocn_frac;
		min_area = min(area_ocn[no][jo*nxo[no]+io], area_wav[nw][lw]);
		if(xarea/min_area > area_ratio_thresh) {

		  wavxocn_area[nw][no][nwxo[nw][no]] = xarea;
		  wavxocn_io[nw][no][nwxo[nw][no]]   = io;
		  wavxocn_jo[nw][no][nwxo[nw][no]]   = jo;
		  wavxocn_iw[nw][no][nwxo[nw][no]]   = iw;
		  wavxocn_jw[nw][no][nwxo[nw][no]]   = jw;
		  if(interp_order == 2) {
		    wavxocn_clon[nw][no][nwxo[nw][no]] = poly_ctrlon ( x_out, y_out, n_out, xw_avg)*ocn_frac;
		    wavxocn_clat[nw][no][nwxo[nw][no]] = poly_ctrlat ( x_out, y_out, n_out )*ocn_frac;
		  }
		  ++(nwxo[nw][no]);
		  if(nwxo[nw][no] > MAXXGRID) mpp_error("nwxo is greater than MAXXGRID, increase MAXXGRID");
		}
	      }
	    }
	  }
	}
      }/* end of la loop */

      mpp_delete_domain2d(&Dom);
      free(is_ocn);
      free(ie_ocn);
      free(js_ocn);
      free(je_ocn);

    } /* end of nw loop */


    time_end = time(NULL);
    if(verbose) printf("on pe %d, calculating waveXocean used %f seconds.\n", mpp_pe(), difftime(time_end, time_start));
    /* calculate the centroid of model grid, as well as land_mask and ocean_mask */
    {
      double **w_area, **w_area2;
      int    nw, lw;
      w_area = (double **)malloc(ntile_wav*sizeof(double *));
      for(nw =0; nw<ntile_wav; nw++) {
	w_area[nw] = (double *)malloc(nxw[nw]*nyw[nw]*sizeof(double));
	for(lw=0; lw<nxw[nw]*nyw[nw]; lw++) {
	  w_area[nw][lw] = 0;
	}
      }

      if(interp_order == 1) {
	for(nw=0; nw<ntile_wav; nw++) {
	  for(no=0; no<ntile_ocn; no++) {
	    int nxgrid;
	    nxgrid = nwxo[nw][no];
	    mpp_sum_int(1, &nxgrid);
	    if(nxgrid > 0) {
	      double *g_area;
	      int    *g_iw, *g_jw;
	      int    ii;

	      g_iw = (int    *)malloc(nxgrid*sizeof(int   ));
	      g_jw = (int    *)malloc(nxgrid*sizeof(int   ));
	      g_area = (double *)malloc(nxgrid*sizeof(double));
	      mpp_gather_field_int   (nwxo[nw][no], wavxocn_io[nw][no], g_iw);
	      mpp_gather_field_int   (nwxo[nw][no], wavxocn_jo[nw][no], g_jw);
	      mpp_gather_field_double(nwxo[nw][no], wavxocn_area[nw][no], g_area);
	      for(i=0; i<nxgrid; i++) {
		ii = g_jw[i]*nxw[nw]+g_iw[i];
		w_area2[nw][ii] += g_area[i];
	      }
	      free(g_iw);
	      free(g_jw);
	      free(g_area);
	    }
	  }
	}
      }
      else { /* interp_order == 2 */
	double **o_area, **o_clon, **o_clat;
	double  *w_clon,  *w_clat;
	int      lo;

	o_area = (double **)malloc(ntile_ocn*sizeof(double *));
	o_clon = (double **)malloc(ntile_ocn*sizeof(double *));
	o_clat = (double **)malloc(ntile_ocn*sizeof(double *));
	for(no =0; no<ntile_ocn; no++) {
          o_area[no] = (double *)malloc(nxo[no]*nyo[no]*sizeof(double));
	  o_clon[no] = (double *)malloc(nxo[no]*nyo[no]*sizeof(double));
	  o_clat[no] = (double *)malloc(nxo[no]*nyo[no]*sizeof(double));
	  for(lo=0; lo<nxo[no]*nyo[no]; lo++) {
            o_area[no][lo] = 0;
	    o_clon[no][lo] = 0;
	    o_clat[no][lo] = 0;
	  }
	}
	for(nw=0; nw<ntile_wav; nw++) {
      	  w_clon = (double *)malloc(nxw[nw]*nyw[nw]*sizeof(double));
	  w_clat = (double *)malloc(nxw[nw]*nyw[nw]*sizeof(double));
	  for(lw=0; lw<nxw[nw]*nyw[nw]; lw++) {
	    w_clon [lw] = 0;
	    w_clat [lw] = 0;
	  }

	  for(no=0; no<ntile_ocn; no++) {
	    int nxgrid;
	    nxgrid = nwxo[nw][no];
	    mpp_sum_int(1, &nxgrid);
	    if(nxgrid > 0) {
	      double *g_area, *g_clon, *g_clat;
	      int    *g_iw,   *g_jw,   *g_io, *g_jo;
	      int    ii;
	      g_iw = (int    *)malloc(nxgrid*sizeof(int   ));
	      g_jw = (int    *)malloc(nxgrid*sizeof(int   ));
	      g_io = (int    *)malloc(nxgrid*sizeof(int   ));
	      g_jo = (int    *)malloc(nxgrid*sizeof(int   ));
	      g_area = (double *)malloc(nxgrid*sizeof(double));
	      g_clon = (double *)malloc(nxgrid*sizeof(double));
	      g_clat = (double *)malloc(nxgrid*sizeof(double));
	      mpp_gather_field_int   (nwxo[nw][no], wavxocn_iw[nw][no], g_iw);
	      mpp_gather_field_int   (nwxo[nw][no], wavxocn_jw[nw][no], g_jw);
	      mpp_gather_field_int   (nwxo[nw][no], wavxocn_io[nw][no], g_io);
	      mpp_gather_field_int   (nwxo[nw][no], wavxocn_jo[nw][no], g_jo);
	      mpp_gather_field_double(nwxo[nw][no], wavxocn_area[nw][no], g_area);
	      mpp_gather_field_double(nwxo[nw][no], wavxocn_clon[nw][no], g_clon);
	      mpp_gather_field_double(nwxo[nw][no], wavxocn_clat[nw][no], g_clat);
	      for(i=0; i<nxgrid; i++) {
		ii = g_jw[i]*nxw[nw]+g_iw[i];
		w_area[nw][ii] += g_area[i];
		w_clon[ii] += g_clon[i];
		w_clat[ii] += g_clat[i];
		ii = g_jo[i]*nxo[no]+g_io[i];
		o_area[no][ii] += g_area[i];
		o_clon[no][ii] += g_clon[i];
		o_clat[no][ii] += g_clat[i];
	      }
	      free(g_iw);
	      free(g_jw);
	      free(g_io);
	      free(g_jo);
	      free(g_area);
	      free(g_clon);
	      free(g_clat);
	    }
	  }

	  for(lw=0; lw<nxw[nw]*nyw[nw]; lw++) {
	    if(w_area[nw][lw] > 0) {
	      w_clon[lw] /= w_area[nw][lw];
	      w_clat[lw] /= w_area[nw][lw];
	    }
	  }

	  for(no=0; no<ntile_ocn; no++) {
	    for(i=0; i<nwxo[nw][no]; i++) {
	      lw = wavxocn_jw[nw][no][i]*nxw[nw] + wavxocn_iw[nw][no][i];
	      wavxocn_diw[nw][no][i] = wavxocn_clon[nw][no][i]/wavxocn_area[nw][no][i] - w_clon[lw];
	      wavxocn_djw[nw][no][i] = wavxocn_clat[nw][no][i]/wavxocn_area[nw][no][i] - w_clat[lw];
	    }
	  }

	  free(w_clon);
	  free(w_clat);
	}


	/* centroid distance from exchange grid to ocean grid */
	for(no=0; no<ntile_ocn; no++) {
	  for(lo=0; lo<nxo[no]*nyo[no]; lo++) {
	    if(o_area[no][lo] > 0) {
	      o_clon[no][lo] /= o_area[no][lo];
	      o_clat[no][lo] /= o_area[no][lo];
	    }
	  }
	  for(nw=0; nw<ntile_wav; nw++) {
	    for(i=0; i<nwxo[nw][no]; i++) {
	      lo = wavxocn_jo[nw][no][i]*nxo[no] + wavxocn_io[nw][no][i];
	      wavxocn_dio[nw][no][i] = wavxocn_clon[nw][no][i]/wavxocn_area[nw][no][i] - o_clon[no][lo];
	      wavxocn_djo[nw][no][i] = wavxocn_clat[nw][no][i]/wavxocn_area[nw][no][i] - o_clat[no][lo];
	    }
	  }
	  free(o_area[no]);
	  free(o_clon[no]);
	  free(o_clat[no]);
	}
	free(o_clon);
	free(o_clat);
	free(o_area);
      }

      /* calculate land/sea fraction for wave grid from wavxocn */
  {
	int    iw, jw;
	double ocn_frac;
	int    id_mask, fid, dims[2];
	char   wav_mask_file[STRING];
	double *mask;
	int ny;

	for(nw=0; nw<ntile_wav; nw++) {
	  mask  = (double *)malloc(nxw[nw]*nyw[nw]*sizeof(double));
	  for(jw=0; jw<nyw[nw]; jw++) for(iw=0; iw<nxw[nw]; iw++) {
	    i = jw*nxw[nw]+iw;
	    mask[i]  = w_area[nw][i]/area_wav[nw][i];
	  }


	  if(ntile_wav > 1)
	    sprintf(wav_mask_file, "wave_mask_tile%d.nc", nw+1);
	  else
	    strcpy(wav_mask_file, "wave_mask.nc");

	  fid = mpp_open(wav_mask_file, MPP_WRITE);
    print_provenance_gv_gca(fid, history, grid_version, clip_method == GREAT_CIRCLE_CLIP );

	  dims[1] = mpp_def_dim(fid, "nx", nxw[nw]);
	  dims[0] = mpp_def_dim(fid, "ny", nyw[nw]);
	  id_mask = mpp_def_var(fid, "mask", MPP_DOUBLE, 2, dims,  2, "standard_name",
				"land/sea fraction at T-cell centers", "units", "none");
	  mpp_end_def(fid);
	  mpp_put_var_value(fid, id_mask, mask);
	  mpp_close(fid);
	  free(mask);
	}


  }
    }


    /* write out the wavxocn */

    for(nw = 0; nw < ntile_wav; nw++) {
       for(no = 0; no < ntile_ocn; no++) {
	 int nxgrid;

	 /* get total number of exchange grid on all the pes */
	 nxgrid = nwxo[nw][no];
	 mpp_sum_int(1, &nxgrid);

	if(nxgrid >0) {
	  size_t start[4], nwrite[4];
	  int *gdata_int;
	  double *gdata_dbl;

	  char contact[STRING];
	  int fid, dim_string, dim_ncells, dim_two, dims[4];
	  int id_contact, id_xgrid_area, n;
	  int id_tile1_cell, id_tile2_cell, id_tile1_dist, id_tile2_dist;

	  for(i=0; i<4; i++) {
	    start[i] = 0; nwrite[i] = 1;
	  }

	  if(same_mosaic)
	    sprintf(wxo_file[nfile_wxo], "wav_%s_%sXocn_%s_%s.nc", wmosaic_name, wtile_name[nw], omosaic_name, otile_name[no]);
	  else
	    sprintf(wxo_file[nfile_wxo], "%s_%sX%s_%s.nc", wmosaic_name, wtile_name[nw], omosaic_name, otile_name[no]);
	  sprintf(contact, "%s:%s::%s:%s", wmosaic_name, wtile_name[nw], omosaic_name, otile_name[no]);

	  fid = mpp_open(wxo_file[nfile_wxo], MPP_WRITE);
    print_provenance_gv_gca(fid, history, grid_version, clip_method == GREAT_CIRCLE_CLIP );

	  dim_string = mpp_def_dim(fid, "string", STRING);
	  dim_ncells = mpp_def_dim(fid, "ncells", nxgrid);
	  dim_two    = mpp_def_dim(fid, "two", 2);
	  if(interp_order == 2) {
	    id_contact = mpp_def_var(fid, "contact", MPP_CHAR, 1, &dim_string, 7, "standard_name", "grid_contact_spec",
				     "contact_type", "exchange", "parent1_cell",
				     "tile1_cell", "parent2_cell", "tile2_cell", "xgrid_area_field", "xgrid_area",
				     "distant_to_parent1_centroid", "tile1_distance", "distant_to_parent2_centroid", "tile2_distance");
	  }
	  else {
	    id_contact = mpp_def_var(fid, "contact", MPP_CHAR, 1, &dim_string, 5, "standard_name", "grid_contact_spec",
				     "contact_type", "exchange", "parent1_cell",
				     "tile1_cell", "parent2_cell", "tile2_cell", "xgrid_area_field", "xgrid_area" );
	  }
	  dims[0] = dim_ncells; dims[1] = dim_two;
	  id_tile1_cell = mpp_def_var(fid, "tile1_cell", MPP_INT, 2, dims, 1, "standard_name", "parent1_cell_indices");
	  id_tile2_cell = mpp_def_var(fid, "tile2_cell", MPP_INT, 2, dims, 1, "standard_name", "parent2_cell_indices");
	  id_xgrid_area = mpp_def_var(fid, "xgrid_area", MPP_DOUBLE, 1, &dim_ncells, 2, "standard_name",
				      "exchange_grid_area", "units", "m2");

	  if(interp_order == 2) {
	    id_tile1_dist = mpp_def_var(fid, "tile1_distance", MPP_DOUBLE, 2, dims, 1, "standard_name", "distance_from_parent1_cell_centroid");
	    id_tile2_dist = mpp_def_var(fid, "tile2_distance", MPP_DOUBLE, 2, dims, 1, "standard_name", "distance_from_parent2_cell_centroid");
	  }
	  mpp_end_def(fid);

	  /* the index will start from 1, instead of 0 ( fortran index) */
	  for(i = 0;i < nwxo[nw][no]; i++) {
	    ++(wavxocn_iw[nw][no][i]);
	    ++(wavxocn_jw[nw][no][i]);
	    ++(wavxocn_io[nw][no][i]);
	    wavxocn_jo[nw][no][i] += 1 - ocn_south_ext; /* one artificial j-level may be added in the south end */
	  }

	  nwrite[0] = strlen(contact);
	  mpp_put_var_value_block(fid, id_contact, start, nwrite, contact);
	  nwrite[0] = nxgrid;

	  gdata_int = (int *)malloc(nxgrid*sizeof(int));
	  gdata_dbl = (double *)malloc(nxgrid*sizeof(double));

	  mpp_gather_field_double(nwxo[nw][no], wavxocn_area[nw][no], gdata_dbl);
	  if(check) {
	    for(n=0; n<nxgrid; n++) wxo_area_sum += gdata_dbl[n];
	  }
	  mpp_put_var_value(fid, id_xgrid_area, gdata_dbl);
	  mpp_gather_field_int(nwxo[nw][no], wavxocn_iw[nw][no], gdata_int);
	  mpp_put_var_value_block(fid, id_tile1_cell, start, nwrite, gdata_int);
	  mpp_gather_field_int(nwxo[nw][no], wavxocn_io[nw][no], gdata_int);
	  mpp_put_var_value_block(fid, id_tile2_cell, start, nwrite, gdata_int);
	  start[1] = 1;
	  mpp_gather_field_int(nwxo[nw][no], wavxocn_jw[nw][no], gdata_int);
	  mpp_put_var_value_block(fid, id_tile1_cell, start, nwrite, gdata_int);
	  mpp_gather_field_int(nwxo[nw][no], wavxocn_jo[nw][no], gdata_int);
	  mpp_put_var_value_block(fid, id_tile2_cell, start, nwrite, gdata_int);
	  if(interp_order == 2) {
	    start[1] = 0;
  	    mpp_gather_field_double(nwxo[nw][no], wavxocn_diw[nw][no], gdata_dbl);
	    mpp_put_var_value_block(fid, id_tile1_dist, start, nwrite, gdata_dbl);
	    mpp_gather_field_double(nwxo[nw][no], wavxocn_dio[nw][no], gdata_dbl);
	    mpp_put_var_value_block(fid, id_tile2_dist, start, nwrite, gdata_dbl);
	    start[1] = 1;
  	    mpp_gather_field_double(nwxo[nw][no], wavxocn_djw[nw][no], gdata_dbl);
	    mpp_put_var_value_block(fid, id_tile1_dist, start, nwrite, gdata_dbl);
	    mpp_gather_field_double(nwxo[nw][no], wavxocn_djo[nw][no], gdata_dbl);
	    mpp_put_var_value_block(fid, id_tile2_dist, start, nwrite, gdata_dbl);
	  }
	  mpp_close(fid);
	  free(gdata_int);
	  free(gdata_dbl);
	  ++nfile_wxo;
	}
       } /* for(no=0; no<ntile_ocn; no++) */

    } /* for(nw=0; nl<ntile_wav; nw++) */

    /* release the memory */
    for(nw=0; nw<ntile_wav; nw++) {
      /* memory used for wave X ocean */
      for(no=0; no<ntile_ocn; no++) {
	free(wavxocn_area[nw][no]);
	free(wavxocn_iw  [nw][no]);
	free(wavxocn_jw  [nw][no]);
	free(wavxocn_io  [nw][no]);
	free(wavxocn_jo  [nw][no]);
	if(interp_order == 2) {
	  free(wavxocn_clon[nw][no]);
	  free(wavxocn_clat[nw][no]);
	  free(wavxocn_diw [nw][no]);
	  free(wavxocn_djw [nw][no]);
	  free(wavxocn_dio [nw][no]);
	  free(wavxocn_djo [nw][no]);
	}
      }
      free(wavxocn_area[nw]);
      free(wavxocn_iw  [nw]);
      free(wavxocn_jw  [nw]);
      free(wavxocn_io  [nw]);
      free(wavxocn_jo  [nw]);
      if(interp_order == 2) {
	free(wavxocn_clon[nw]);
	free(wavxocn_clat[nw]);
	free(wavxocn_diw [nw]);
	free(wavxocn_djw [nw]);
	free(wavxocn_dio [nw]);
	free(wavxocn_djo [nw]);
      }

      free(nwxo[nw]);
    }
    free(nwxo);
    free(wavxocn_area);
    free(wavxocn_iw  );
    free(wavxocn_jw  );
    free(wavxocn_io  );
    free(wavxocn_jo  );
    if(interp_order == 2) {
      free(wavxocn_clon);
      free(wavxocn_clat);
      free(wavxocn_diw );
      free(wavxocn_djw );
      free(wavxocn_dio );
      free(wavxocn_djo );
    }

    if(mpp_pe() == mpp_root_pe() && verbose) printf("\nNOTE from make_coupler_mosaic: Complete the process to create exchange grids "
					 "between wave and ocean\n" );

  }

  /* tiling check */
  if(mpp_pe() == mpp_root_pe()) {
    int n, i;
    double axo_area_frac;
    double axl_area_frac;
    double axo_area_frac_nest;
    double axl_area_frac_nest;
    double tiling_area;
    double tiling_area_nest;
    double atm_area_sum, atm_area_sum_nest;

    if(check) {
      double max_ratio, cur_ratio;
      int  max_ia, max_ja, max_ta;

      /* make sure each atmosphere grid area equal sum of the axo area and axl area (atm_xarea) */
      max_ta = -1;
      max_ia = -1;
      max_ja = -1;
      max_ratio = -1;

      for(n=0; n<ntile_atm; n++) {
	for(i=0; i<nxa[n]*nya[n]; i++) {
	  cur_ratio = fabs(atm_xarea[n][i] - area_atm[n][i])/area_atm[n][i];
	  if(cur_ratio > 1.e-5) {
	    printf("at tile =%d, i=%d, j=%d, ratio=%g, area1=%g, area2=%g\n",
		   n+1, i%nxa[n], i/nxa[n], cur_ratio,  area_atm[n][i], atm_xarea[n][i] );
	  }
	  if(cur_ratio > max_ratio) {
	    max_ratio = cur_ratio;
	    max_ta = n;
	    max_ja = i/nxa[n];
	    max_ia = i%nxa[n];
	  }
	}
      }

      printf("The maximum area mismatch is at tile=%d, i=%d, j=%d, ratio=%g, area_atm_grid=%g, area_atm_xgrid=%g\n",
	     max_ta+1, max_ia, max_ja, max_ratio, area_atm[max_ta][max_ja*nxa[max_ta]+max_ia],
	     atm_xarea[max_ta][max_ja*nxa[max_ta]+max_ia]);
      /* for cubic grid, when number of model points is odd, some grid cell area will be negative,
	 need to think about how to solve this issue in the future */
      /*      atm_area_sum = 4*M_PI*RADIUS*RADIUS; */
      atm_area_sum = 0;
      atm_area_sum_nest = 0;
      for(n=0; n<ntile_atm; n++) {
	if( n== tile_nest )
	  for(i=0; i<nxa[n]*nya[n]; i++) atm_area_sum_nest += area_atm[n][i];
	else
	  for(i=0; i<nxa[n]*nya[n]; i++) atm_area_sum += area_atm[n][i];
      }
      axo_area_frac = axo_area_sum/atm_area_sum*100;
      axl_area_frac = axl_area_sum/atm_area_sum*100;
      tiling_area   = (atm_area_sum-axo_area_sum-axl_area_sum)/atm_area_sum*100;
      printf("\nNOTE: axo_area_sum is %f and ocean fraction is %f%%\n", axo_area_sum, axo_area_frac);
      printf("NOTE: axl_area_sum is %f and land  fraction is %f%%\n", axl_area_sum, axl_area_frac);
      printf("NOTE: tiling error is %f%%\n", tiling_area );
      if(wmosaic) {
	printf("\nNOTE: axo_area_sum is %f\n", axo_area_sum);
        printf("\nNOTE: wxo_area_sum is %f\n", wxo_area_sum);
      }
      if(tile_nest >=0) {
	axo_area_frac_nest = axo_area_sum_nest/atm_area_sum_nest*100;
	axl_area_frac_nest = axl_area_sum_nest/atm_area_sum_nest*100;
	tiling_area_nest   = (atm_area_sum_nest-axo_area_sum_nest-axl_area_sum_nest)/atm_area_sum_nest*100;
	printf("\nNOTE: axo_area_sum_nest is %f and ocean fraction is %f%%\n", axo_area_sum_nest, axo_area_frac_nest);
	printf("NOTE: axl_area_sum_nest is %f and land  fraction is %f%%\n", axl_area_sum_nest, axl_area_frac_nest);
	printf("NOTE: tiling error for nest region is %f%%\n", tiling_area_nest );
      }

    }
  }

  /*Fianlly create the coupler mosaic file mosaic_name.nc */
  {
    int fid, dim_string, dim_axo, dim_lxo, dim_axl, dim_axw, dim_wxo, dims[4], n;
    size_t start[4], nwrite[4];
    int id_lmosaic_dir, id_lmosaic_file, id_omosaic_dir, id_omosaic_file, id_wmosaic_dir;
    int id_amosaic_dir, id_amosaic_file, id_otopog_dir, id_otopog_file, id_wmosaic_file;
    int id_xgrids_dir, id_axo_file, id_lxo_file, id_axl_file, id_wxo_file;
    int id_amosaic, id_lmosaic, id_omosaic, id_wmosaic;

    fid = mpp_open(mosaic_file, MPP_WRITE);
    print_provenance_gv_gca(fid, history, grid_version, clip_method == GREAT_CIRCLE_CLIP );

    dim_string = mpp_def_dim(fid, "string", STRING);
    if(nfile_axo >0) dim_axo = mpp_def_dim(fid, "nfile_aXo", nfile_axo);
    if(nfile_axl >0) dim_axl = mpp_def_dim(fid, "nfile_aXl", nfile_axl);
    if(nfile_lxo >0) dim_lxo = mpp_def_dim(fid, "nfile_lXo", nfile_lxo);
    if(nfile_wxo >0) dim_wxo = mpp_def_dim(fid, "nfile_wXo", nfile_wxo);

    id_amosaic_dir  = mpp_def_var(fid, "atm_mosaic_dir", MPP_CHAR, 1, &dim_string,
				  1, "standard_name", "directory_storing_atmosphere_mosaic");
    id_amosaic_file = mpp_def_var(fid, "atm_mosaic_file", MPP_CHAR, 1, &dim_string,
				  1, "standard_name", "atmosphere_mosaic_file_name");
    id_amosaic      = mpp_def_var(fid, "atm_mosaic", MPP_CHAR, 1, &dim_string,
				  1, "standard_name", "atmosphere_mosaic_name");
    id_lmosaic_dir  = mpp_def_var(fid, "lnd_mosaic_dir", MPP_CHAR, 1, &dim_string,
				  1, "standard_name", "directory_storing_land_mosaic");
    id_lmosaic_file = mpp_def_var(fid, "lnd_mosaic_file", MPP_CHAR, 1, &dim_string,
				  1, "standard_name", "land_mosaic_file_name");
    id_lmosaic      = mpp_def_var(fid, "lnd_mosaic", MPP_CHAR, 1, &dim_string,
				  1, "standard_name", "land_mosaic_name");
    id_omosaic_dir  = mpp_def_var(fid, "ocn_mosaic_dir", MPP_CHAR, 1, &dim_string,
				  1, "standard_name", "directory_storing_ocean_mosaic");
    id_omosaic_file = mpp_def_var(fid, "ocn_mosaic_file", MPP_CHAR, 1, &dim_string,
				  1, "standard_name", "ocean_mosaic_file_name");
    id_omosaic      = mpp_def_var(fid, "ocn_mosaic", MPP_CHAR, 1, &dim_string,
				  1, "standard_name", "ocean_mosaic_name");
    if(wmosaic) {
      id_wmosaic_dir  = mpp_def_var(fid, "wav_mosaic_dir", MPP_CHAR, 1, &dim_string,
				    1, "standard_name", "directory_storing_wave_mosaic");
      id_wmosaic_file = mpp_def_var(fid, "wav_mosaic_file", MPP_CHAR, 1, &dim_string,
				    1, "standard_name", "wave_mosaic_file_name");
      id_wmosaic      = mpp_def_var(fid, "wav_mosaic", MPP_CHAR, 1, &dim_string,
				  1, "standard_name", "wave_mosaic_name");
    }

    id_otopog_dir   = mpp_def_var(fid, "ocn_topog_dir", MPP_CHAR, 1, &dim_string,
				  1, "standard_name", "directory_storing_ocean_topog");
    id_otopog_file  = mpp_def_var(fid, "ocn_topog_file", MPP_CHAR, 1, &dim_string,
				  1, "standard_name", "ocean_topog_file_name");
    /* since exchange grid is created in this tool, we may add command line option to specify where to store the output */
    /*    id_xgrids_dir = mpp_def_var(fid, "xgrids_dir", MPP_CHAR, 1, &dim_string,
	  1, "standard_name", "directory_storing_xgrids"); */
    if(nfile_axo >0) {
      dims[0] = dim_axo; dims[1] = dim_string;
      id_axo_file = mpp_def_var(fid, "aXo_file", MPP_CHAR, 2, dims, 1, "standard_name", "atmXocn_exchange_grid_file");
    }
    if(nfile_axl >0) {
      dims[0] = dim_axl; dims[1] = dim_string;
      id_axl_file = mpp_def_var(fid, "aXl_file", MPP_CHAR, 2, dims, 1, "standard_name", "atmXlnd_exchange_grid_file");
    }
    if(nfile_lxo >0) {
      dims[0] = dim_lxo; dims[1] = dim_string;
      id_lxo_file = mpp_def_var(fid, "lXo_file", MPP_CHAR, 2, dims, 1, "standard_name", "lndXocn_exchange_grid_file");
    }
    if(nfile_wxo >0) {
      dims[0] = dim_wxo; dims[1] = dim_string;
      id_wxo_file = mpp_def_var(fid, "wXo_file", MPP_CHAR, 2, dims, 1, "standard_name", "wavXocn_exchange_grid_file");
    }
    mpp_end_def(fid);
    for(i=0; i<4; i++) { start[i] = 0; nwrite[i] = 1; }

    nwrite[0] = strlen(amosaic_dir);
    mpp_put_var_value_block(fid, id_amosaic_dir, start, nwrite, amosaic_dir);
    nwrite[0] = strlen(amosaic_file);
    mpp_put_var_value_block(fid, id_amosaic_file, start, nwrite, amosaic_file);
    nwrite[0] = strlen(amosaic_name);
    mpp_put_var_value_block(fid, id_amosaic, start, nwrite, amosaic_name);
    nwrite[0] = strlen(lmosaic_dir);
    mpp_put_var_value_block(fid, id_lmosaic_dir, start, nwrite, lmosaic_dir);
    nwrite[0] = strlen(lmosaic_file);
    mpp_put_var_value_block(fid, id_lmosaic_file, start, nwrite, lmosaic_file);
    nwrite[0] = strlen(lmosaic_name);
    mpp_put_var_value_block(fid, id_lmosaic, start, nwrite, lmosaic_name);
    nwrite[0] = strlen(omosaic_dir);
    mpp_put_var_value_block(fid, id_omosaic_dir, start, nwrite, omosaic_dir);
    nwrite[0] = strlen(omosaic_file);
    mpp_put_var_value_block(fid, id_omosaic_file, start, nwrite, omosaic_file);
    nwrite[0] = strlen(omosaic_name);
    mpp_put_var_value_block(fid, id_omosaic, start, nwrite, omosaic_name);
    if(wmosaic) {
      nwrite[0] = strlen(wmosaic_dir);
      mpp_put_var_value_block(fid, id_wmosaic_dir, start, nwrite, wmosaic_dir);
      nwrite[0] = strlen(wmosaic_file);
      mpp_put_var_value_block(fid, id_wmosaic_file, start, nwrite, wmosaic_file);
      nwrite[0] = strlen(wmosaic_name);
      mpp_put_var_value_block(fid, id_wmosaic, start, nwrite, wmosaic_name);
    }

    nwrite[0] = strlen(otopog_dir);
    mpp_put_var_value_block(fid, id_otopog_dir, start, nwrite, otopog_dir);
    nwrite[0] = strlen(otopog_file);
    mpp_put_var_value_block(fid, id_otopog_file, start, nwrite, otopog_file);
    nwrite[0] = 1;

    for(n=0; n<nfile_axo; n++) {
      start[0] = n; nwrite[1] = strlen(axo_file[n]);
      mpp_put_var_value_block(fid, id_axo_file, start, nwrite, axo_file[n]);
    }
    for(n=0; n<nfile_axl; n++) {
      start[0] = n; nwrite[1] = strlen(axl_file[n]);
      mpp_put_var_value_block(fid, id_axl_file, start, nwrite, axl_file[n]);
    }
    for(n=0; n<nfile_lxo; n++) {
      start[0] = n; nwrite[1] = strlen(lxo_file[n]);
      mpp_put_var_value_block(fid, id_lxo_file, start, nwrite, lxo_file[n]);
    }

    for(n=0; n<nfile_wxo; n++) {
      start[0] = n; nwrite[1] = strlen(wxo_file[n]);
      mpp_put_var_value_block(fid, id_wxo_file, start, nwrite, wxo_file[n]);
    }

    mpp_close(fid);
  }

  if(mpp_pe()== mpp_root_pe() && verbose){
     if(nbad == 0){
       printf("\n***** Congratulation! You have successfully run make_coupler_mosaic\n");}
     else{
       printf("\n***** You have run make_coupler_mosaic, but there were at least %d grid cells with issues.\n",nbad);
       mpp_error("make_coupler_mosaic: omask is not equal ocn_frac");}
  }
  mpp_end();

  return 0;

} /* main */
