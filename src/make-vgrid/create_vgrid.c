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
#include "create_vgrid.h"
#include "tool_util.h"

/***********************************************************************
   void create_vgrid(int nbnds, double *bnds, int *nz, char center, double *zeta)
   This routine is used to create vertical grid. The created grid is on super grid.
   the refinement is assumed to be 2.

   INPUT:
      nbnds:  number of vertical regions for varying resolution.      
      bnds:   boundaries for defining vertical regions of varying resolution.
              The size of bnds is nz.
      nz:     Number of model grid points for each vertical regions of varying resolution.
              The size of nz will be nz-1
      center: grid cell center location. Its value can be 'N', 'T' or 'C' with default
              value 'N'. When the value is 'N', supergrid location will be calculated.
              When the value is 'T', model grid corner location will be calculated,
              other grid location of the supergrid will be derived through T-cell
              centered condition. When the value is 'C', model grid corner location will
              be calculated, other grid location of the supergrid will be derived through
              C-cell centered condition.

   OUTPUT:
      zeta:   created vertical grid location.

***********************************************************************/

void create_vgrid(int nbnds, double *bnds, int *nz, double *dbnds, double stretch, int use_legacy, double *zeta, int *np, const char *center)
{
  int i;
  double *zb=NULL;


  if(use_legacy)
    zb = compute_grid_bound_legacy(nbnds, bnds, dbnds, stretch, np, center);
  else
    zb = compute_grid_bound(nbnds, bnds, nz, np, center);
  
  for(i=0; i<(*np+1); i++) zeta[i] =zb[i];

  free(zb);
    
  
}; /* create_vgrid */
