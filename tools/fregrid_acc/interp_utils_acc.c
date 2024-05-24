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
#include <openacc.h>
#include "globals_acc.h"
#include "general_utils_acc.h"

/*******************************************************************************
void copy_grid_to_device( const int itile, Grid_config *grid )
Copies lat lon coordinates to device
*******************************************************************************/
void copy_grid_to_device_acc( const int npoints, const double *lat, const double *lon )
{

#pragma acc enter data copyin(lon[:npoints], lat[:npoints])

}

/*******************************************************************************
void copy_interp_to_device( Xgrid_config *interp )
Copies the interp struct to device
*******************************************************************************/
void copy_xgrid_to_device_acc( const int ntiles_in, const int ntiles_out, const Xgrid_config *xgrid,
                               const unsigned int opcode )
{

#pragma acc enter data copyin(xgrid[:ntiles_out])
  for(int n=0 ; n<ntiles_out; n++) {

#pragma acc enter data copyin( xgrid[n].per_intile[:ntiles_in] )

    for(int m=0 ; m<ntiles_in ; m++) {

      int nxcells = xgrid[n].per_intile[m].nxcells;
#pragma acc enter data copyin( xgrid[n].per_intile[m].input_parent_lon_indices[:nxcells], \
                               xgrid[n].per_intile[m].input_parent_lat_indices[:nxcells], \
                               xgrid[n].per_intile[m].output_parent_lon_indices[:nxcells], \
                               xgrid[n].per_intile[m].output_parent_lat_indices[:nxcells], \
                               xgrid[n].per_intile[m].xcell_area[:nxcells] )
        if( opcode & CONSERVE_ORDER2) {
#pragma acc enter data copyin( xgrid[n].per_intile[m].dcentroid_lon[:nxcells], \
                               xgrid[n].per_intile[m].dcentroid_lat[:nxcells])
        }

    } // ntiles_in
  }//ntiles_out

} //end copy_interp_to_device

/*******************************************************************************
void get_bounding_indices
gets indices for kat that overlap with the ref_lat
TODO: THIS FUNCTION NEEDS A UNIT TEST
*******************************************************************************/
void get_bounding_indices(const int ref_nx, const int ref_ny, const int nx, const int ny,
                          const double *ref_lat, const double *lat, int *jstart, int *jend, int *new_ny)
{

  int ref_min, ref_max, yy;
  int jstart2, jend2;

  ref_min = minval_double((ref_ny+1)*(ref_nx+1), ref_lat);
  ref_max = maxval_double((ref_ny+1)*(ref_nx+1), ref_lat);

  jstart2 = ny;
  jend2 = -1;

  for(int j=0; j<ny+1; j++) {
    for(int i=0; i<nx+1; i++) {
      yy = lat[j*(nx+1)+i];
      if( yy > ref_min ) jstart2 = min(jstart2, j);
      if( yy < ref_max ) jend2   = max(jend2, j);
    }
  }

  *jstart = max(0, jstart2-1);
  *jend   = min(ny-1, jend2+1);
  *new_ny = jend2-jstart2+1;

}
/*******************************************************************************
void get_mask
get mask between input and output grid.
*******************************************************************************/
void create_mask_on_device(const int mask_size, double **mask)
{

  double *pmask;
  pmask = mask;

  pmask = (double *)malloc(mask_size*sizeof(double));

#pragma acc enter data create(pmask[:mask_size])
#pragma acc parallel loop independent present(pmask[:mask_size])
  for( int i=0 ; i<mask_size; i++) pmask[i]=1.0;

}
