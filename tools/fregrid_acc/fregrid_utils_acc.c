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
#include "globals.h"

/*******************************************************************************
void copy_grid_to_device( const int itile, Grid_config *grid )
Copies lat lon coordinates to device
*******************************************************************************/
void copy_grid_to_device( const int itile, const Grid_config *grid )
{

  int nxp, nyp;

  nxp = grid[itile].nxc +1;
  nyp = grid[itile].nyc +1;

#pragma acc enter data copyin(grid[itile])
#pragma acc enter data copyin(grid[itile].lonc[0:nxp*nyp], grid[itile].latc[0:nxp*nyp])

}

/*******************************************************************************
void copy_interp_to_device( Interp_config *interp )
Copies the interp struct to device
*******************************************************************************/
void copy_interp_to_device( const int ntiles_in, const int ntiles_out, const Interp_config *interp,
                            const unsigned int opcode )
{

#pragma acc enter data copyin(interp[:ntiles_out])
  for(int n=0 ; n<ntiles_out; n++) {

#pragma acc enter data copyin( interp[n].interp_mini[:ntiles_in] )

    for(int m=0 ; m<ntiles_in ; m++) {

      int nxgrid_acc = interp[n].interp_mini[m].nxgrid;
#pragma acc enter data copyin( interp[n].interp_mini[m].i_in[:nxgrid_acc], \
                               interp[n].interp_mini[m].j_in[:nxgrid_acc], \
                               interp[n].interp_mini[m].i_out[:nxgrid_acc], \
                               interp[n].interp_mini[m].j_out[:nxgrid_acc], \
                               interp[n].interp_mini[m].area[:nxgrid_acc] )
        if( opcode & CONSERVE_ORDER2) {
#pragma acc enter data copyin( interp[n].interp_mini[m].di_in[:nxgrid_acc], \
                               interp[n].interp_mini[m].dj_in[:nxgrid_acc])
        }

    } // ntiles_in
  }//ntiles_out

} //end copy_interp_to_device
