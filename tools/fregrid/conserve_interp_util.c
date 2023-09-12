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
#include "conserve_interp_util.h"

void read_remap_file( int ntiles_in, int ntiles_out, Grid_config *grid_out,
          Interp_config *interp, unsigned int opcode){

  int *i_in=NULL, *j_in=NULL, *i_out=NULL, *j_out=NULL;
  double *xgrid_area=NULL, *tmp_area=NULL, *xgrid_clon=NULL, *xgrid_clat=NULL;
  int n, i;
  size_t nxgrid;
  double garea;

  garea = 4*M_PI*RADIUS*RADIUS;

  for(n=0; n<ntiles_out; n++) {
    if( interp[n].file_exist ) { /* reading from file */
      int *t_in, *ind;
      int fid, vid;

      nxgrid     = read_mosaic_xgrid_size(interp[n].remap_file);
      i_in       = (int    *)malloc(nxgrid   * sizeof(int   ));
      j_in       = (int    *)malloc(nxgrid   * sizeof(int   ));
      i_out      = (int    *)malloc(nxgrid   * sizeof(int   ));
      j_out      = (int    *)malloc(nxgrid   * sizeof(int   ));
      xgrid_area = (double *)malloc(nxgrid   * sizeof(double));
      if(opcode & CONSERVE_ORDER2) {
        xgrid_clon = (double *)malloc(nxgrid   * sizeof(double));
        xgrid_clat = (double *)malloc(nxgrid   * sizeof(double));
      }
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
      free(i_in);
      free(j_in);
      free(i_out);
      free(j_out);
      free(xgrid_area);
      if(opcode & CONSERVE_ORDER2) {
        free(xgrid_clon);
        free(xgrid_clat);
      }
    }//if read from file
  } // ntiles
  if(mpp_pe() == mpp_root_pe())printf("NOTE: Finish reading index and weight for conservative interpolation from file.\n");

}//end read_regrid_weights


void malloc_arrays( int nsize, int **i_in, int **j_in, int **i_out, int **j_out,
                    double **xgrid_area, double **xgrid_clon, double **xgrid_clat )
{
  *i_in       = (int *) malloc(nsize * sizeof(int   ));
  *j_in       = (int *) malloc(nsize * sizeof(int   ));
  *i_out      = (int *) malloc(nsize * sizeof(int   ));
  *j_out      = (int *) malloc(nsize * sizeof(int   ));
  *xgrid_area = (double *) malloc(nsize * sizeof(double));
  *xgrid_clon = (double *) malloc(nsize * sizeof(double));
  *xgrid_clat = (double *) malloc(nsize * sizeof(double));
}
