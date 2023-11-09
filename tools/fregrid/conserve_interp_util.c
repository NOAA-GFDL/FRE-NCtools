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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <netcdf.h>
#include "constant.h"
#include "globals.h"
#include "mpp.h"
#include "mpp_io.h"
#include "read_mosaic.h"
#include "mosaic_util.h"
#include "conserve_interp_util.h"

/*******************************************************************************
  void read_remap_file( int ntiles_in, int ntiles_out, Grid_config *grid_out,
                        Interp_config *interp, unsigned int opcode)
  Reads in the weight/remap file if provided
*******************************************************************************/
void read_remap_file( int ntiles_in, int ntiles_out, Grid_config *grid_out,
          Interp_config *interp, unsigned int opcode){

  int *i_in=NULL, *j_in=NULL, *i_out=NULL, *j_out=NULL;
  double *xgrid_area=NULL, *tmp_area=NULL, *xgrid_clon=NULL, *xgrid_clat=NULL;
  int n, i;
  int zero=0;
  size_t nxgrid;
  double garea;

  garea = 4*M_PI*RADIUS*RADIUS;

  for(n=0; n<ntiles_out; n++) {
    if( interp[n].file_exist ) { /* reading from file */
      int *t_in, *ind;
      int fid, vid;

      nxgrid     = read_mosaic_xgrid_size(interp[n].remap_file);
      malloc_xgrid_arrays(nxgrid, &i_in, &j_in, &i_out, &j_out, &xgrid_area, &xgrid_clon, &xgrid_clat);
      t_in       = (int    *)malloc(nxgrid*sizeof(int   ));
      ind        = (int    *)malloc(nxgrid*sizeof(int   ));
      if(opcode & CONSERVE_ORDER1)
        read_mosaic_xgrid_order1(interp[n].remap_file, i_in, j_in, i_out, j_out, xgrid_area);
      else
        read_mosaic_xgrid_order2(interp[n].remap_file, i_in, j_in, i_out, j_out, xgrid_area, xgrid_clon, xgrid_clat);

      /*--- rescale the xgrid area */
      for(i=0; i<nxgrid; i++) xgrid_area[i] *= garea;
      fid = mpp_open(interp[n].remap_file, MPP_READ);
      vid = mpp_get_varid(fid, "tile1");
      mpp_get_var_value(fid, vid, t_in);
      mpp_close(fid);
      /*distribute the exchange grid on each pe according to target grid index*/
      interp[n].nxgrid = 0;
      for(i=0; i<nxgrid; i++) {
        if( i_out[i] <= grid_out[n].iec && i_out[i] >= grid_out[n].isc &&
            j_out[i] <= grid_out[n].jec && j_out[i] >= grid_out[n].jsc )
          ind[interp[n].nxgrid++] = i;
      }
      interp[n].i_in   = (int    *)malloc(interp[n].nxgrid*sizeof(int   ));
      interp[n].j_in   = (int    *)malloc(interp[n].nxgrid*sizeof(int   ));
      interp[n].i_out  = (int    *)malloc(interp[n].nxgrid*sizeof(int   ));
      interp[n].j_out  = (int    *)malloc(interp[n].nxgrid*sizeof(int   ));
      interp[n].area   = (double *)malloc(interp[n].nxgrid*sizeof(double));
      interp[n].t_in   = (int    *)malloc(interp[n].nxgrid*sizeof(int   ));

      for(i=0; i< interp[n].nxgrid; i++) {
        interp[n].i_in [i] = i_in [ind[i]];
        interp[n].j_in [i] = j_in [ind[i]];
        interp[n].t_in [i] = t_in [ind[i]] - 1;
        interp[n].i_out[i] = i_out[ind[i]] - grid_out[n].isc;
        interp[n].j_out[i] = j_out[ind[i]] - grid_out[n].jsc;
        interp[n].area [i] = xgrid_area[ind[i]];
      }
      if(opcode & CONSERVE_ORDER2) {
        interp[n].di_in   = (double *)malloc(interp[n].nxgrid*sizeof(double));
        interp[n].dj_in   = (double *)malloc(interp[n].nxgrid*sizeof(double));
        for(i=0; i< interp[n].nxgrid; i++) {
          interp[n].di_in[i] = xgrid_clon[ind[i]];
          interp[n].dj_in[i] = xgrid_clat[ind[i]];
        }
      }
      free(t_in);
      free(ind);
      malloc_xgrid_arrays(zero, &i_in, &j_in, &i_out, &j_out, &xgrid_area, &xgrid_clon, &xgrid_clat);
    }//if read from file
  } // ntiles
  if(mpp_pe() == mpp_root_pe())printf("NOTE: Finish reading index and weight for conservative interpolation from file.\n");

}//end read_regrid_weights

/*******************************************************************************
  void malloc_xgrid_arrays( int nsize, int **i_in, int **j_in, int **i_out, int **j_out,
                            double **xgrid_area, double **xgrid_clon, double **xgrid_clat )
  allocates arrays that will hold exchange grid information
*******************************************************************************/
void malloc_xgrid_arrays( int nsize, int **i_in, int **j_in, int **i_out, int **j_out,
                          double **xgrid_area, double **xgrid_clon, double **xgrid_clat )
{

  // free if malloc-ed
  if(*i_in!=NULL) {
    free(*i_in);
    *i_in=NULL;
  }
  if(*j_in!=NULL) {
    free(*j_in);
    *j_in=NULL;
  }
  if(*i_out!=NULL) {
    free(*i_out);
    *i_out=NULL;
  }
  if(*j_out!=NULL) {
    free(*j_out);
    *j_out=NULL;
  }
  if(*xgrid_area!=NULL) {
    free(*xgrid_area);
    *xgrid_area=NULL;
  }
  if(*xgrid_clon!=NULL) {
    free(*xgrid_clon);
    *xgrid_clon=NULL;
  }
  if(*xgrid_clat!=NULL) {
    free(*xgrid_clat);
    *xgrid_clat=NULL;
  }

  if(nsize>0) {
    *i_in       = (int *) malloc(nsize * sizeof(int   ));
    *j_in       = (int *) malloc(nsize * sizeof(int   ));
    *i_out      = (int *) malloc(nsize * sizeof(int   ));
    *j_out      = (int *) malloc(nsize * sizeof(int   ));
    *xgrid_area = (double *) malloc(nsize * sizeof(double));
    *xgrid_clon = (double *) malloc(nsize * sizeof(double));
    *xgrid_clat = (double *) malloc(nsize * sizeof(double));
  }

}

/*******************************************************************************
  void get_CellStruct
  Gathers exchange grid information from all ranks and
  stores information in cell_in structure.
  Cell_in stores exchange grid information corresponding to each input parent cell
*******************************************************************************/
void get_CellStruct(const int tile_in, const int nx_in, const int nxgrid, int *i_in, int *j_in,
                    double *xgrid_area, double *xgrid_clon, double *xgrid_clat,
                    CellStruct *cell_in)
{

  int g_nxgrid;
  int *g_i_in, *g_j_in;
  double *g_area, *g_clon, *g_clat;

  g_nxgrid = nxgrid;
  mpp_sum_int(1, &g_nxgrid);
  if(g_nxgrid > 0) {
    g_i_in = (int    *)malloc(g_nxgrid*sizeof(int   ));
    g_j_in = (int    *)malloc(g_nxgrid*sizeof(int   ));
    g_area = (double *)malloc(g_nxgrid*sizeof(double));
    g_clon = (double *)malloc(g_nxgrid*sizeof(double));
    g_clat = (double *)malloc(g_nxgrid*sizeof(double));
    mpp_gather_field_int   (nxgrid, i_in,       g_i_in);
    mpp_gather_field_int   (nxgrid, j_in,       g_j_in);
    mpp_gather_field_double(nxgrid, xgrid_area, g_area);
    mpp_gather_field_double(nxgrid, xgrid_clon, g_clon);
    mpp_gather_field_double(nxgrid, xgrid_clat, g_clat);
    for(int i=0; i<g_nxgrid; i++) {
      int ii;
      ii = g_j_in[i]*nx_in+g_i_in[i];
      cell_in[tile_in].area[ii] += g_area[i];
      cell_in[tile_in].clon[ii] += g_clon[i];
      cell_in[tile_in].clat[ii] += g_clat[i];
    }
    free(g_i_in);
    free(g_j_in);
    free(g_area);
    free(g_clon);
    free(g_clat);
  }

}
/*******************************************************************************
  void get_interp
  stores exchange grid information to the interp structure
********************************************************************************/
void get_interp( const int opcode, const int nxgrid, Interp_config *interp, const int m, const int n,
                 const int *i_in, const int *j_in, const int *i_out, const int *j_out,
                 const double *xgrid_clon, const double *xgrid_clat, const double *xgrid_area )
{

  int nxgrid_prev;
  int i;

  if(nxgrid > 0) {
    nxgrid_prev = interp[n].nxgrid;
    interp[n].nxgrid += nxgrid;
    if(nxgrid_prev == 0 ) {
      interp[n].i_in   = (int    *)malloc(interp[n].nxgrid*sizeof(int   ));
      interp[n].j_in   = (int    *)malloc(interp[n].nxgrid*sizeof(int   ));
      interp[n].i_out  = (int    *)malloc(interp[n].nxgrid*sizeof(int   ));
      interp[n].j_out  = (int    *)malloc(interp[n].nxgrid*sizeof(int   ));
      interp[n].area   = (double *)malloc(interp[n].nxgrid*sizeof(double));
      interp[n].t_in   = (int    *)malloc(interp[n].nxgrid*sizeof(int   ));
      for(i=0; i<interp[n].nxgrid; i++) {
        interp[n].t_in [i] = m;
        interp[n].i_in [i] = i_in [i];
        interp[n].j_in [i] = j_in [i];
        interp[n].i_out[i] = i_out[i];
        interp[n].j_out[i] = j_out[i];
        interp[n].area[i]  = xgrid_area[i];
      }
      if(opcode & CONSERVE_ORDER2) {
        interp[n].di_in   = (double *)malloc(interp[n].nxgrid*sizeof(double));
        interp[n].dj_in   = (double *)malloc(interp[n].nxgrid*sizeof(double));
        for(i=0; i<interp[n].nxgrid; i++) {
          interp[n].di_in [i] = xgrid_clon[i]/xgrid_area[i];
          interp[n].dj_in [i] = xgrid_clat[i]/xgrid_area[i];
        }
      }
    }
    else {
      interp[n].i_in   = (int    *)realloc(interp[n].i_in,  interp[n].nxgrid*sizeof(int   ));
      interp[n].j_in   = (int    *)realloc(interp[n].j_in,  interp[n].nxgrid*sizeof(int   ));
      interp[n].i_out  = (int    *)realloc(interp[n].i_out, interp[n].nxgrid*sizeof(int   ));
      interp[n].j_out  = (int    *)realloc(interp[n].j_out, interp[n].nxgrid*sizeof(int   ));
      interp[n].area   = (double *)realloc(interp[n].area,  interp[n].nxgrid*sizeof(double));
      interp[n].t_in   = (int    *)realloc(interp[n].t_in,  interp[n].nxgrid*sizeof(int   ));
      for(i=0; i<nxgrid; i++) {
        interp[n].t_in [nxgrid_prev+i] = m;
        interp[n].i_in [nxgrid_prev+i] = i_in [i];
        interp[n].j_in [nxgrid_prev+i] = j_in [i];
        interp[n].i_out[nxgrid_prev+i] = i_out[i];
        interp[n].j_out[nxgrid_prev+i] = j_out[i];
        interp[n].area [nxgrid_prev+i] = xgrid_area[i];
      }
      if(opcode & CONSERVE_ORDER2) {
        interp[n].di_in   = (double *)realloc(interp[n].di_in, interp[n].nxgrid*sizeof(double));
        interp[n].dj_in   = (double *)realloc(interp[n].dj_in, interp[n].nxgrid*sizeof(double));
        for(i=0; i<nxgrid; i++) {
          interp[n].di_in [i+nxgrid_prev] = xgrid_clon[i]/xgrid_area[i];
          interp[n].dj_in [i+nxgrid_prev] = xgrid_clat[i]/xgrid_area[i];
        }
      }
    }
  }  /* if(nxgrid>0) */

}
/*******************************************************************************
  void get_jstart_jend
  get the starting and ending indices of the input grid that
  overlaps with the output grid
********************************************************************************/
void get_jstart_jend( const int nx_out, const int ny_out, const int nx_in, const int ny_in,
                      const double *lat_out, const double *lat_in,
                      int *jstart, int *jend, int *ny_now )
{

  double y_min, y_max, yy ;
  int i, j;

  y_min = minval_double((nx_out+1)*(ny_out+1), lat_out);
  y_max = maxval_double((nx_out+1)*(ny_out+1), lat_out);
  *jstart = ny_in; *jend = -1;
  for(j=0; j<=ny_in; j++) for(i=0; i<=nx_in; i++) {
      yy = lat_in[j*(nx_in+1)+i];
      if( yy > y_min ) {
        if(j < *jstart ) *jstart = j;
      }
      if( yy < y_max ) {
        if(j > *jend ) *jend = j;
      }

    }
  *jstart = max(0, *jstart-1);
  *jend   = min(ny_in-1, *jend+1);
  *ny_now = *jend-*jstart+1;

}
