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
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include "mpp.h"
#include "mpp_io.h"
#include "constant.h"
#include "tool_util.h"
#include "read_mosaic.h"

char *usage[] = {
  "",
  "                                                                                ",
  "                                                                                ",
  "                    Usage of make_land_domain                                   ",
  "                                                                                ",
  "   make_land_domain --grid_file grid_file --land_restart land_restart           ",
  "                    [--output_file output_file]                                 ",
  "                                                                                ",
  " REQUIRED:                                                                      ",
  "                                                                                ",
  " --grid_file grid_file       coupled mosaic file name                           ",
  " --land_restart land_restart specify the source file of land restart file. The  ",
  "                                                                                ",
  "                             file name does not tile# and normally in the for   ",
  "                             'land.res'                                         ",
  "                                                                                ",
  " OPTIONAL FLAGS                                                                 ",
  "                                                                                ",
  " --output_file output_file   Specify the output file name. The default value is ",
  "                             'land_domain.nc'.                                  ",
  "                                                                                ",
  NULL };


int main(int argc, char* argv[])
{
  int c,n;
  int errflg = (argc == 1);
  int    option_index = 0;
  char *grid_file=NULL;
  char *land_restart=NULL;
  char outfile[STRING] = "land_domain.nc";
  char griddir[STRING], solo_mosaic[STRING], land_mosaic[STRING];
  char mosaic_name[STRING];
  int  **land_mask=NULL, **land_ntile=NULL;
  int  *nland_face=NULL, *grid_index=NULL, *grid_ntile=NULL;
  int  fid, vid;
  int  nfaces,face,nx,ny,nlands;

    /*
   * process command line
   */

  static struct option long_options[] = {
    {"grid_file",    required_argument, NULL, 'g'},
    {"land_restart", required_argument, NULL, 'l'},
    {"output_file",  required_argument, NULL, 'o'},
    {0, 0, 0, 0},
  };

  mpp_init(&argc, &argv);

  if(mpp_npes() != 1) mpp_error("make_land_domain: this tool must be run on single procesor");
  
    while ((c = getopt_long(argc, argv, "", long_options, &option_index)) != -1)
    switch (c) {
    case 'g':
      grid_file = optarg;
      break;      
    case 'l':
      land_restart = optarg;
      break;
    case 'o':
      strcpy(outfile, optarg);
      break;      
    case '?':
      errflg++;
      break;
    }

  if( !grid_file ) errflg++;
  if( !land_restart ) errflg++;
  if (errflg) {
    char **u = usage;
    while (*u) { fprintf(stderr, "%s\n", *u); u++; }
    mpp_error("make_land: check the command line arguments");
  }


    /* write out arguments */
  if(mpp_pe() == mpp_root_pe()) {
    printf("grid_file    is %s\n", grid_file);
    printf("land_restart is %s\n", land_restart);
  }

  get_file_path(grid_file, griddir);
  fid = mpp_open(grid_file, MPP_READ);
  vid = mpp_get_varid(fid, "lnd_mosaic_file");
  mpp_get_var_value(fid, vid, solo_mosaic);
  sprintf(land_mosaic, "%s/%s", griddir, solo_mosaic);

  /* get the grid resolution */
  {
    int  *nxl=NULL, *nyl=NULL;
    nfaces = read_mosaic_ntiles(land_mosaic);
    nxl = (int *)malloc(nfaces*sizeof(int));
    nyl = (int *)malloc(nfaces*sizeof(int));
  
    read_mosaic_grid_sizes(land_mosaic, nxl, nyl);

    /* nx, ny should have the same value on each face */
    for(n=1; n<nfaces; n++) {
      if(nxl[n] != nxl[0] || nyl[n] != nyl[0])
	mpp_error("remap_land: all the faces of source grid should have the same number of grid points");
    }  
    nx = nxl[0];
    ny = nyl[0];
    free(nxl);
    free(nyl);
  }
  
  /* read the exchange grid information and get the land/sea mask of land model*/
  {
    int nfile_aXl;    
    double *land_area;    
    land_area = (double *)malloc(nx*ny*sizeof(double));
    for(n=0; n<nx*ny; n++) land_area[n] = 0;


    nfile_aXl = mpp_get_dimlen( fid, "nfile_aXl");
    vid = mpp_get_varid(fid, "lnd_mosaic");
    mpp_get_var_value(fid, vid, mosaic_name);
  
    vid = mpp_get_varid(fid, "aXl_file");

    land_mask = (int **)malloc(nfaces*sizeof(int *));
    nland_face = (int *)malloc(nfaces*sizeof(int));
    nlands = 0;
    for(face=0; face<nfaces; face++) {
      char face_name[STRING];

    
      nland_face[face] = 0;
      land_mask[face] = (int *)malloc(nx*ny*sizeof(int));
      for(n=0; n<nx*ny; n++) land_mask[face][n] = 0;
      sprintf(face_name, "%s_tile%d", mosaic_name, face+1);
      for(n=0; n<nx*ny; n++) land_area[n] = 0.0;  
      for(n=0; n<nfile_aXl; n++) {
	size_t start[4], nread[4];
	int nxgrid;
	char aXl_file[STRING], filepath[STRING];
	start[0] = n;
	start[1] = 0;
	nread[0] = 1;
	nread[1] = STRING;
	mpp_get_var_value_block(fid, vid, start, nread, aXl_file);
	if( !strstr(aXl_file, face_name) ) continue;
	sprintf(filepath, "%s/%s", griddir, aXl_file);
	nxgrid = read_mosaic_xgrid_size(filepath);
	if(nxgrid>0) {
	  int l;
	  int *i1, *j1, *i2, *j2;
	  double *area;

	  i1 = (int *)malloc(nxgrid*sizeof(int));
	  j1 = (int *)malloc(nxgrid*sizeof(int));
	  i2 = (int *)malloc(nxgrid*sizeof(int));
	  j2 = (int *)malloc(nxgrid*sizeof(int));
	  area = (double *)malloc(nxgrid*sizeof(double));
	  read_mosaic_xgrid_order1(filepath, i1, j1, i2, j2, area);
	  for(l=0; l<nxgrid; l++) land_area[j2[l]*nx+i2[l]] += (area[l]*4*M_PI*RADIUS*RADIUS);
	  free(i1);
	  free(j1);
	  free(i2);
	  free(j2);
	  free(area);
	}
      }
      for(n=0; n<nx*ny; n++) {
	if(land_area[n] >0) {
	  land_mask[face][n] = 1;
	  nland_face[face]++;
	  nlands++;
	}
      }
    }
    free(land_area);
  }
  
  mpp_close(fid);
  /* compute total number of land points */

  printf("total number of land points is %d\n", nlands);
  land_ntile = (int **)malloc(nfaces*sizeof(int *));
  for(face=0; face<nfaces; face++) {
    land_ntile[face] = (int *)malloc(nx*nx*sizeof(int));
    for(n=0; n<nx*ny; n++) land_ntile[face][n] = 0;
  }
  
  /* open the land restart file and get tile_index */
  {
    int fid2, vid2;
    int npts;

    npts = nx*ny;
    
    for(face=0; face<nfaces; face++) {
      char land_restart_file[STRING];
      int num_idx, i;
      int *idx;
      if(nfaces==1)
	strcpy(land_restart, land_restart_file);
      else
	sprintf(land_restart_file, "%s.tile%d.nc", land_restart, face+1);
      fid2 = mpp_open(land_restart_file, MPP_READ);
      vid2 = mpp_get_varid(fid2, "tile_index");
      num_idx = mpp_get_dimlen(fid2, "tile_index");
      idx = (int *)malloc(num_idx*sizeof(int));
      mpp_get_var_value(fid2, vid2, idx);
      for(n=0; n<num_idx; n++) {
	if(idx[n]<0) mpp_error("make_land: idx should be non-negatve");
	i = idx[n]%npts;
	(land_ntile[face][i])++;
      }
      
      free(idx);
      mpp_close(fid2);
    }
  }

  /*make sure land_mask and land_ntile are consistent */
  for(face=0; face<nfaces; face++) {
    for(n=0; n<nx*ny; n++) {
      if(land_mask[face][n] > 0 && land_ntile[face][n] == 0)
	mpp_error("make_land: land_mask > 0 but land_ntile = 0");
      else if(land_mask[face][n] == 0 && land_ntile[face][n] > 0)
	mpp_error("make_land: land_mask = 0 but land_ntile > 0");
    }
  }
    
  
  /* pack the grid tile information into unstructured grid */
  grid_index = (int *)malloc(nlands*sizeof(int));
  grid_ntile = (int *)malloc(nlands*sizeof(int));
  for(n=0; n<nlands; n++) {
    grid_index[n] = -1;
    grid_ntile[n] = 0;
  }
  
  {
    int pos;
    
    pos = 0;
    for(face=0; face<nfaces; face++) {
      for(n=0; n<nx*ny; n++) {
	if(land_mask[face][n] > 0) {
	  grid_index[pos] = n;
	  grid_ntile[pos] = land_ntile[face][n];
	  pos++;
	}
      }
    }

  }
  
  /* write out the data */
  {
    char history[1280];
    int dim_nfaces, dim_nlands;
    int id_nland_face, id_grid_index, id_grid_ntile;
    int n;
    
    strcpy(history,argv[0]);
    for(n=1;n<argc;n++) {
      strcat(history, " ");
      strcat(history, argv[n]);
    }

    
    fid = mpp_open(outfile, MPP_WRITE);
    dim_nfaces = mpp_def_dim(fid, "nfaces", nfaces);
    dim_nlands =  mpp_def_dim(fid, "nlands", nlands);

    id_nland_face = mpp_def_var(fid, "nland_face", MPP_INT, 1, &dim_nfaces, 2, "standard_name",
			    "number of land points at each face", "units", "none");
    id_grid_index = mpp_def_var(fid, "grid_index", MPP_INT, 1, &dim_nlands, 2, "standard_name",
			    "grid index in structured grid", "units", "none");
    id_grid_ntile = mpp_def_var(fid, "grid_ntile", MPP_INT, 1, &dim_nlands, 2, "standard_name",
			    "number of tiles in each grid", "units", "none");

    mpp_def_global_att(fid, "history", history);
    mpp_end_def(fid);

    mpp_put_var_value(fid, id_nland_face, nland_face);
    mpp_put_var_value(fid, id_grid_index, grid_index);
    mpp_put_var_value(fid, id_grid_ntile, grid_ntile);
    mpp_close(fid);
  }

  /* free the memory */
  for(face=0; face<nfaces; face++) {
    free(land_mask[face]);
    free(land_ntile[face]);
  }
  free(land_mask);
  free(land_ntile);
  free(nland_face);
  free(grid_index);
  free(grid_ntile);
  
}
