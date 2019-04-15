/*
  The program create remap file from Cubic sphere grid to Cubic sphere grid.
  The ratio of grid size of input and output grid must be integer. The purpose
  of this tool is to create remap file on single processor for very high resoluton
  grid ( It will large processor count and long time to create such file ). Currently
  this tool can only remap file for data on grid cell center.

 AUTHOR: Zhi Liang (Zhi.Liang@noaa.gov)
          NOAA Geophysical Fluid Dynamics Lab, Princeton, NJ
 
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
 
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 
  For the full text of the GNU General Public License,
  write to: Free Software Foundation, Inc.,
            675 Mass Ave, Cambridge, MA 02139, USA.  
-----------------------------------------------------------------------
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
#include <time.h>
#include "constant.h"
#include "read_mosaic.h"
#include "mpp_io.h"
#include "mpp.h"
#include "mosaic_util.h"
#include "tool_util.h"
#include "create_xgrid.h"

#define EPSLN (1.e-10)
#define D2R (M_PI/180)
#define CONSERVE_ORDER1 1
#define CONSERVE_ORDER2 2

char *usage[] = {
  ""
  " make_remap_file  --input_mosaic input_mosaic --output_mosaic output_mosaic ",
  "                  --remap_file remap_file [--interp_method method]          ",
  "                                                                            ",
  " make_remap_file create remap file from Cubic sphere grid to Cubic sphere   ",
  " grid. The ratio of grid size of input and output grid must be integer.     ",
  " The purpose of this tool is to create remap file on single processor for   ",
  " very high resoluton grid ( It will large processor count and long time     ",
  " to create such file ).  Currently this tool can only remap file for data   ",
  " on grid cell center.                                                       ",
  "                                                                            ",
  "make_remap_file takes the following flags:                                                    ",
  "                                                                                      ",
  "REQUIRED:                                                                             ",
  "                                                                                      ",
  "--input_mosaic  input_mosaic  specify the input mosaic information. This file         ",
  "                              contains list of tile files which specify the grid      ",
  "                              information for each tile.                              ",
  "                                                                                      ",
  "--output_mosaic output_mosaic specify the output mosaic information. This file        ",
  "                              contains list of tile files which specify the grid      ",
  "                              information for each tile. If output_mosaic is not      ",
  "                              specified, nlon and nlat must be specified.             ",
  "                                                                                      ",
  "--remap_file   remap_file     specify the file name that saves remapping information. ",
  "                              If remap_file is specified and the file does not exist, ",
  "                              remapping information will be calculated ans stored in  ",
  "                              remap_file. If remap_file is specified and the file     ",
  "                              exists, remapping information will be read from         ",
  "                              remap_file.                                             ",
  "                                                                                      ",
  "OPTIONAL FLAGS                                                                        ",
  "                                                                                      ",  
  "--interp_method interp_method specify the remapping algorithm to be used. Default is  ",
  "                              'conserve_order1'. Currently only 'conserve_order1',    ",
  "                              'conserve_order2' remapping scheme are implemented in   ",
  "                              this tool.                                              "
  "                                                                                      ",
  " Example:                                                                             ",
  "                                                                                      ",
  " make_remap_file --input_mosaic C3072_mosaic.nc --output_mosaic C6144_mosaic.nc       ",
  "                 --remap_file remap_C3072_to_C6144.nc --interp_method conserve_order2 ",
  "                                                                                      ",
   NULL};

int main(int argc, char* argv[])
{
  char    *mosaic_in=NULL;            /* input mosaic file name */
  char    *mosaic_out=NULL;           /* input mosaic file name */
  char    dir_in[STRING];               /* input file location */
  char    dir_out[STRING];              /* output file location */
  int     ntiles_in = 0;              /* number of tiles in input mosaic */
  int     ntiles_out = 0;             /* number of tiles in output mosaic */
  char    *remap_file = NULL;
  char    interp_method[STRING] = "conserve_order1";

  char    history[512];
  char remap_file_base[STRING];
  int remap_method;
  int c, option_index, n;
  int mid_in, mid_out, len, nxgrid;
  double *x_in=NULL, *y_in=NULL, *tmp_in=NULL;
  double *x_out=NULL, *y_out=NULL, *tmp_out=NULL;
  double *xarea=NULL, *xclon=NULL, *xclat=NULL;
  double *area_in=NULL, *clon_in=NULL, *clat_in=NULL;
  double *di_in=NULL, *dj_in;
  int *i_in=NULL, *j_in=NULL, *t_in=NULL;
  int *i_out=NULL, *j_out=NULL;
  int nx_in, ny_in, ni_in, nj_in;
  int nx_out, ny_out, ni_out, nj_out; 
  int errflg = (argc == 1);

  static struct option long_options[] = {
    {"input_mosaic",     required_argument, NULL, 'i'},
    {"output_mosaic",    required_argument, NULL, 'o'},
    {"remap_file",       required_argument, NULL, 'r'},
    {"interp_method",    required_argument, NULL, 'm'},
    {"help",             no_argument,       NULL, 'h'},
    {0, 0, 0, 0},
  };

  mpp_init(&argc, &argv);
  if(mpp_npes() > 1) mpp_error("make_remap_file: this tool must be run on single MPI rank");

  while ((c = getopt_long(argc, argv, "", long_options, &option_index)) != -1) {
    switch (c) {
    case 'i':
      mosaic_in  = optarg;
      break;
    case 'o':
      mosaic_out = optarg;
      break;
    case 'r':
      remap_file = optarg;
      break;    
    case 'm':
      strcpy(interp_method, optarg);
      break;    
    case '?':
      errflg++;
      break;
    }
  }

  if (errflg) {
    char **u = usage;
    while (*u) { fprintf(stderr, "%s\n", *u); u++; }
    exit(2);
  }      
  /* check the arguments */
  if( !mosaic_in  ) mpp_error("make_remap_file: input_mosaic is not specified");
  if( !mosaic_out ) mpp_error("make_remap_file: output_mosaic is not specified");
  if( !remap_file ) mpp_error("make_remap_file: remap_file is not specified");
  
  remap_method = CONSERVE_ORDER1;
  if( strcmp(interp_method, "conserve_order1") && strcmp(interp_method, "conserve_order2"))
    mpp_error("make_remap_file: interp_method should be conserve_order1 or conserve_order2");
  if(!strcmp(interp_method, "conserve_order2") ) {
    printf("****make_remap_file: second order conservative remap information will be created\n");
    remap_method = CONSERVE_ORDER2;
  }
  else
    printf("****make_remap_file: first order conservative remap information will be created\n");

  /* define history to be the history in the grid file */
  {
    int i;
    strcpy(history,argv[0]);

    for(i=1;i<argc;i++) {
      strcat(history, " ");
      strcat(history, argv[i]);
    }
  }
  /* get the mosaic information of input and output mosaic*/
  mid_in = mpp_open(mosaic_in, MPP_READ);
  ntiles_in = mpp_get_dimlen(mid_in, "ntiles");
  mid_out = mpp_open(mosaic_out, MPP_READ);
  ntiles_out = mpp_get_dimlen(mid_out, "ntiles");
  if( ntiles_out != ntiles_in)  mpp_error("make_remap_file: ntiles_out not equal to ntiles_in"); 

  /*get the path of input_mosaic and output_mosaic, assume the grid file are in the same directory */
  get_file_path(mosaic_in, dir_in);
  get_file_path(mosaic_out, dir_out);

  /* get the remap_file base */
  {
    len = strlen(remap_file);
    if(len >= STRING) mpp_error("make_remap_file: length of remap_file should be less than STRING");
    if( strcmp(remap_file+len-3, ".nc")==0 ) {
      strncpy(remap_file_base, remap_file, len-3);
      remap_file_base[len-3] = 0;
    }
    else
      strcpy(remap_file_base, remap_file);
  }
    
  /* loop through ntiles to get the input and output grid, also create the remap file */
  for(n=0; n<ntiles_in; n++) {
    char   filename[STRING];
    char   grid_file[STRING];
    char   my_remap_file[STRING];
    double lon_in_avg;
    double px_out[20], py_out[20];
    double px_in[20], py_in[20];
    size_t start[4], nread[4], nwrite[4];
    int nxp_in, nyp_in, nip_in, njp_in;
    int nxp_out, nyp_out, nip_out, njp_out;
    int i, j, ratio, pos, i1, j1;
    int vid, fid, n_out, n_in;
    int dim_string, dim_ncells, dim_two, dims[2];
    int id_tile1, id_tile1_cell, id_tile2_cell, id_xgrid_area;
    int id_tile1_dist;

    
    /**********************************************************************
                         Get input grid
    ***********************************************************************/
    start[0] = n; start[1] = 0; nread[0] = 1; nread[1] = STRING;
    vid = mpp_get_varid(mid_in, "gridfiles");
    mpp_get_var_value_block(mid_in, vid, start, nread, filename);
    sprintf(grid_file, "%s/%s", dir_in, filename);
    fid = mpp_open(grid_file, MPP_READ);

    nx_in = mpp_get_dimlen(fid, "nx");
    ny_in = mpp_get_dimlen(fid, "ny");
    if(nx_in%2) mpp_error("make_remap_file: the size of dimension nx in input grid should be even (on supergrid)");
    if(ny_in%2) mpp_error("make_remap_file: the size of dimension ny in input grid should be even (on supergrid)");
    /* Assume all the tiles have the same size */
    if(n>0 && (nx_in != 2*ni_in || ny_in != 2*nj_in)) mpp_error("make_remap_file: all the tiles should have same grid size");
    ni_in = nx_in/2;
    nj_in = ny_in/2;
    nxp_in = nx_in+1;
    nyp_in = ny_in+1;
    nip_in = ni_in+1;
    njp_in = nj_in+1;
    if(n ==0) {
      x_in = (double *)malloc(nip_in*njp_in*sizeof(double));
      y_in = (double *)malloc(nip_in*njp_in*sizeof(double));
      tmp_in = (double*)malloc(nxp_in*nyp_in*sizeof(double));
    }
    vid = mpp_get_varid(fid, "x");
    mpp_get_var_value(fid, vid, tmp_in);
    for(j=0; j<njp_in; j++) for(i=0; i<nip_in; i++) {
      x_in[j*nip_in+i] = tmp_in[2*j*nxp_in+2*i]*D2R;
    }
    vid = mpp_get_varid(fid, "y");
    mpp_get_var_value(fid, vid, tmp_in);
    for(j=0; j<njp_in; j++) for(i=0; i<nip_in; i++) {
      y_in[j*nip_in+i] = tmp_in[2*j*nxp_in+2*i]*D2R;
    }    
    mpp_close(fid);

    /**********************************************************************
                         Get output grid
    ***********************************************************************/
    start[0] = n; start[1] = 0; nread[0] = 1; nread[1] = STRING;
    vid = mpp_get_varid(mid_out, "gridfiles");
    mpp_get_var_value_block(mid_out, vid, start, nread, filename);
    sprintf(grid_file, "%s/%s", dir_out, filename);
    fid = mpp_open(grid_file, MPP_READ);

    nx_out = mpp_get_dimlen(fid, "nx");
    ny_out = mpp_get_dimlen(fid, "ny");
    if(nx_out%2) mpp_error("make_remap_file: the size of dimension nx in output grid should be even (on supergrid)");
    if(ny_out%2) mpp_error("make_remap_file: the size of dimension ny in output grid should be even (on supergrid)");
    /* Assume all the tiles have the same size */
    if(n>1 && (nx_out != 2*ni_out || ny_out != 2*nj_out)) mpp_error("make_remap_file: all the tiles should have same grid size");
    ni_out = nx_out/2;
    nj_out = ny_out/2;
    nxp_out = nx_out+1;
    nyp_out = ny_out+1;
    nip_out = ni_out+1;
    njp_out = nj_out+1;
    if(n ==0) {
      x_out = (double *)malloc(nip_out*njp_out*sizeof(double));
      y_out = (double *)malloc(nip_out*njp_out*sizeof(double));
      tmp_out = (double*)malloc(nxp_out*nyp_out*sizeof(double));
    }
    vid = mpp_get_varid(fid, "x");
    mpp_get_var_value(fid, vid, tmp_out);
    for(j=0; j<njp_out; j++) for(i=0; i<nip_out; i++) {
      x_out[j*nip_out+i] = tmp_out[2*j*nxp_out+2*i]*D2R;
    }
    vid = mpp_get_varid(fid, "y");
    mpp_get_var_value(fid, vid, tmp_out);
    for(j=0; j<njp_out; j++) for(i=0; i<nip_out; i++) {
      y_out[j*nip_out+i] = tmp_out[2*j*nxp_out+2*i]*D2R;
    }    
    mpp_close(fid);

    /**********************************************************************
            check the ratio of input and output grid size is integer
            Currently we assume nx_out > nx_in. Also check the x_in/y_in
            match some of the points of x_out/y_out
    **********************************************************************/
    if(nx_out%nx_in !=0) mpp_error("make_remap_file: nx_out/nx_in is not integer");
    ratio = nx_out/nx_in;
    for(j=0; j<njp_in; j++) for(i=0; i<nip_in; i++) {
      if(fabs(x_in[j*nip_in+i]-x_out[j*ratio*nip_out+i*ratio]) > EPSLN)
	mpp_error("make_remap_file: input grid x does not match output grid");
      if(fabs(y_in[j*nip_in+i]-y_out[j*ratio*nip_out+i*ratio]) > EPSLN)
	mpp_error("make_remap_file: input grid y does not match output grid");
    }
    
    /**********************************************************************
                  Create remap information, number of exchange grid will
                  equal number of output grid.
    **********************************************************************/
    if(n==0) { /* all the tiles have the same size, only do memory allocation at the first tile */
      nxgrid = ni_out*nj_out;
      i_in   = (int    *)malloc(nxgrid*sizeof(int   ));
      j_in   = (int    *)malloc(nxgrid*sizeof(int   ));
      i_out  = (int    *)malloc(nxgrid*sizeof(int   ));
      j_out  = (int    *)malloc(nxgrid*sizeof(int   ));
      xarea   = (double *)malloc(nxgrid*sizeof(double));
      t_in   = (int    *)malloc(nxgrid*sizeof(int   ));
      area_in = (double *)malloc(ni_in*nj_in*sizeof(double));
      for(i=0; i<ni_in*nj_in; i++) area_in[i] = 0.0;
      if(remap_method == CONSERVE_ORDER2) {
	xclon = (double *)malloc(nxgrid*sizeof(double));
	xclat =  (double *)malloc(nxgrid*sizeof(double));
	di_in = (double *)malloc(nxgrid*sizeof(double));
	dj_in = (double *)malloc(nxgrid*sizeof(double));
	clon_in = (double *)malloc(ni_in*nj_in*sizeof(double));
	clat_in = (double *)malloc(ni_in*nj_in*sizeof(double));
      }
    }
    for(i=0; i<ni_in*nj_in; i++) area_in[i] = 0.0;
    if(remap_method == CONSERVE_ORDER2) {
      for(i=0; i<ni_in*nj_in; i++) {
	clon_in[i] = 0.0;
	clat_in[i] = 0.0;
      }
    }
    pos = 0;
    for(j=0; j<nj_in; j++) for(i=0; i<ni_in; i++) {
	px_in[0] = x_in[j*nip_in+i];
	py_in[0] = y_in[j*nip_in+i];
	px_in[1] = x_in[j*nip_in+i+1];
	py_in[1] = y_in[j*nip_in+i+1];
	px_in[2] = x_in[(j+1)*nip_in+i+1];
	py_in[2] = y_in[(j+1)*nip_in+i+1];
	px_in[3] = x_in[(j+1)*nip_in+i];
	py_in[3] = y_in[(j+1)*nip_in+i];
	n_in = fix_lon(px_in, py_in, 4, M_PI);
	lon_in_avg = avgval_double(n_in, px_in);
      for(j1=0; j1< ratio; j1++) for(i1=0; i1<ratio; i1++) {
	i_in[pos] = i;
	j_in[pos] = j;
	i_out[pos] = i*ratio + i1;
	j_out[pos] = j*ratio + j1;
	t_in[pos]  = n;
	/* calculate exchange grid area */
	px_out[0] = x_out[j_out[pos]*nip_out+i_out[pos]];
	py_out[0] = y_out[j_out[pos]*nip_out+i_out[pos]];
	px_out[1] = x_out[j_out[pos]*nip_out+i_out[pos]+1];
	py_out[1] = y_out[j_out[pos]*nip_out+i_out[pos]+1];
	px_out[2] = x_out[(j_out[pos]+1)*nip_out+i_out[pos]+1];
	py_out[2] = y_out[(j_out[pos]+1)*nip_out+i_out[pos]+1];
	px_out[3] = x_out[(j_out[pos]+1)*nip_out+i_out[pos]];
	py_out[3] = y_out[(j_out[pos]+1)*nip_out+i_out[pos]];
	n_out = fix_lon(px_out, py_out, 4, M_PI);
	xarea[pos] = poly_area (px_out, py_out, n_out );
	area_in[j*ni_in+i] += xarea[pos];
	if(remap_method == CONSERVE_ORDER2) {
	  xclon[pos] = poly_ctrlon(px_out, py_out, n_out, lon_in_avg);
	  xclat[pos] = poly_ctrlat (px_out, py_out, n_out );
	  clon_in[j*ni_in+i] += xclon[pos];
	  clat_in[j*ni_in+i] += xclat[pos];
	}
	pos++;
      }
    }

    if(pos != nxgrid) {
      printf("pos=%d, nxgrid=%d\n", pos, nxgrid);
      mpp_error("make_remap_file: pos != nxgrid");
    }

    if(remap_method == CONSERVE_ORDER2) {
      for(i=0; i<ni_in*nj_in; i++) {
	clon_in[i] /= area_in[i];
	clat_in[i] /= area_in[i];
      }
      for(i=0; i<nxgrid; i++) {
	xclon[i] /= xarea[i];
	xclat[i] /= xarea[i];
	di_in[i] = xclon[i] - clon_in[j_in[i]*ni_in+i_in[i]];
	dj_in[i] = xclat[i] - clat_in[j_in[i]*ni_in+i_in[i]];
      }
    }
    /* write out the remap information */
    if(ntiles_in > 1) 
      sprintf(my_remap_file, "%s.tile%d.nc", remap_file_base, n+1);
    else
      sprintf(my_remap_file, "%s.nc", remap_file_base);

    /*convert from c-index to fortran index */
    for(i=0; i<nxgrid; i++) {
      t_in[i]++;  
      i_in[i]++;
      j_in[i]++;
      i_out[i]++;
      j_out[i]++;
    }
    
    fid = mpp_open(my_remap_file, MPP_WRITE);
    mpp_def_global_att(fid, "history", history);
    dim_string = mpp_def_dim(fid, "string", STRING);
    dim_ncells = mpp_def_dim(fid, "ncells", nxgrid);
    dim_two    = mpp_def_dim(fid, "two", 2);
    dims[0] = dim_ncells; dims[1] = dim_two;
    id_tile1      = mpp_def_var(fid, "tile1",      NC_INT, 1, &dim_ncells, 1,
				"standard_name", "tile_number_in_mosaic1");
    id_tile1_cell = mpp_def_var(fid, "tile1_cell", NC_INT, 2, dims, 1,
				"standard_name", "parent_cell_indices_in_mosaic1");
    id_tile2_cell = mpp_def_var(fid, "tile2_cell", NC_INT, 2, dims, 1,
				"standard_name", "parent_cell_indices_in_mosaic2");
    id_xgrid_area = mpp_def_var(fid, "xgrid_area", NC_DOUBLE, 1, &dim_ncells, 2,
				"standard_name", "exchange_grid_area", "units", "m2");
    if(remap_method == CONSERVE_ORDER2) id_tile1_dist = mpp_def_var(fid, "tile1_distance", NC_DOUBLE, 2, dims, 1,
								    "standard_name", "distance_from_parent1_cell_centroid");
    mpp_end_def(fid);
    for(i=0; i<4; i++) {
      start[i] = 0; nwrite[i] = 1;
    }
    nwrite[0] = nxgrid;
    mpp_put_var_value(fid, id_tile1, t_in);

    mpp_put_var_value_block(fid, id_tile1_cell, start, nwrite, i_in);
    mpp_put_var_value_block(fid, id_tile2_cell, start, nwrite, i_out);

    start[1] = 1;
    mpp_put_var_value_block(fid, id_tile1_cell, start, nwrite, j_in);
    mpp_put_var_value_block(fid, id_tile2_cell, start, nwrite, j_out);

    mpp_put_var_value(fid, id_xgrid_area, xarea);

    if(remap_method == CONSERVE_ORDER2) {
      start[1] = 0;
      mpp_put_var_value_block(fid, id_tile1_dist, start, nwrite, di_in);
      start[1] = 1;
      mpp_put_var_value_block(fid, id_tile1_dist, start, nwrite, dj_in);
    }
	  
    mpp_close(fid);
  }

  /*release the memory */
  free(x_in);
  free(y_in);
  free(tmp_in);
  free(x_out);
  free(y_out);
  free(tmp_out);
  free(i_in);
  free(j_in);
  free(i_out);
  free(j_out);
  free(xarea);
  free(t_in);
  if(remap_method == CONSERVE_ORDER2) {
    free(xclon);
    free(xclat);
    free(di_in);
    free(dj_in);
    free(clon_in);
    free(clat_in);
  }
    
  printf("Successfully running make_remap_file and the following output file are generated.\n");
      
  mpp_end();
  return 0;
  
} /* end of main */
  

  
