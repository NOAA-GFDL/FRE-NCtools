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
#include <math.h>
#include <string.h>
#include "fregrid_util_acc.h"
#include "mpp.h"
#include "mpp_io.h"
#include "mosaic_util.h"
#include "globals.h"
#include "interp.h"
#include "create_xgrid.h"

#define D2R (M_PI/180)
#define EPSLN10 (1.e-10)

void get_input_output_cell_area_acc(int ntiles_in, Grid_config *grid_in, int ntiles_out, Grid_config *grid_out, unsigned int opcode)
{
  double *lonc=NULL, *latc=NULL;
  int n, nx, ny, ind1, ind2, i, j, halo;

  halo = grid_in->halo;
  for(n=0; n<ntiles_in; n++) {
    nx = grid_in[n].nx;
    ny = grid_in[n].ny;

    grid_in[n].cell_area = (double *)malloc(nx*ny*sizeof(double));
    lonc = (double *)malloc((nx+1)*(ny+1)*sizeof(double));
    latc = (double *)malloc((nx+1)*(ny+1)*sizeof(double));

    for(j=0; j<=ny; j++) for(i=0; i<=nx; i++) {
        ind1 = j*(nx+1)+i;
        ind2 = (j+halo)*(nx+1+2*halo)+i+halo;
        lonc[ind1] = grid_in[n].lonc[ind2];
        latc[ind1] = grid_in[n].latc[ind2];
      }

    /*calculate grid_in cell area */
    if( opcode & GREAT_CIRCLE )
      get_grid_great_circle_area(&nx, &ny, lonc, latc, grid_in[n].cell_area);
    else
      get_grid_area(&nx, &ny, lonc, latc, grid_in[n].cell_area);

    free(lonc);
    free(latc);
  }

  for(n=0; n<ntiles_out; n++) {
    nx = grid_out[n].nxc;
    ny = grid_out[n].nyc;

    grid_out[n].cell_area = (double *)malloc(nx*ny*sizeof(double));

    /*calculate grid_in cell area */
    if( opcode & GREAT_CIRCLE )
      get_grid_great_circle_area(&nx, &ny, grid_out[n].lonc, grid_out[n].latc, grid_out[n].cell_area);
    else
      get_grid_area(&nx, &ny, grid_out[n].lonc, grid_out[n].latc, grid_out[n].cell_area);

  }

}

/*******************************************************************************
  void get_output_grid_by_size(Mosaic_config *mosaic, int nlon, int nlat, int finer_steps, unsigned int opcode)
  calculate output grid based on nlon, nlat and finer steps.

*******************************************************************************/
void get_output_grid_by_size_acc(int ntiles, Grid_config *grid, double lonbegin, double lonend, double latbegin, double latend,
                                 int nlon, int nlat, int finer_steps, int center_y, unsigned int opcode)
{
  double      dlon, dlat, lon_fine, lat_fine, lon_range, lat_range;
  int         nx_fine, ny_fine, i, j, layout[2];
  int nxc, nyc, ii, jj;

  if(ntiles !=1) mpp_error("fregrid_utils: ntiles of output mosaic must be 1 for bilinear interpolation");
  if(finer_steps && !(opcode&BILINEAR)) mpp_error("fregrid_util: finer_steps must be 0 when interp_method is not bilinear");

  grid->nx      = nlon;
  grid->ny      = nlat;
  grid->nx_fine = pow(2,finer_steps)*nlon;
  grid->ny_fine = pow(2,finer_steps)*(nlat-1)+1;
  nx_fine       = grid->nx_fine;
  ny_fine       = grid->ny_fine;
  lon_range     = lonend - lonbegin;
  lat_range     = latend - latbegin;
  grid->is_tripolar = 0;
  grid->lont1D = (double *)malloc(nlon*sizeof(double));
  grid->latt1D = (double *)malloc(nlat*sizeof(double));
  grid->lonc1D = (double *)malloc((nlon+1)*sizeof(double));
  grid->latc1D = (double *)malloc((nlat+1)*sizeof(double));

  dlon=lon_range/nlon;
  for(i=0; i<nlon; i++) grid->lont1D[i]  = (lonbegin + (i + 0.5)*dlon)*D2R;
  for(i=0; i<=nlon; i++) grid->lonc1D[i] = (lonbegin + i*dlon)*D2R;

  layout[0] = 1;
  layout[1] = mpp_npes();
  mpp_define_domain2d(grid->nx, grid->ny, layout, 0, 0, &(grid->domain));
  mpp_get_compute_domain2d(grid->domain, &(grid->isc), &(grid->iec), &(grid->jsc), &(grid->jec));
  grid->nxc = grid->iec - grid->isc + 1;
  grid->nyc = grid->jec - grid->jsc + 1;
  nxc       = grid->nxc;
  nyc       = grid->nyc;
  if(center_y) {
    dlat=lat_range/nlat;
    for(j=0; j<nlat; j++) grid->latt1D[j] = (latbegin+(j+0.5)*dlat)*D2R;
    for(j=0; j<=nlat; j++) grid->latc1D[j] = (latbegin+j*dlat)*D2R;
  }
  else {
    dlat=lat_range/(nlat-1);
    for(j=0; j<nlat; j++) grid->latt1D[j] = (latbegin+j*dlat)*D2R;
    for(j=0; j<=nlat; j++) grid->latc1D[j] = (latbegin+(j-0.5)*dlat)*D2R;
  }

  if(opcode & BILINEAR) {
    grid->latt1D_fine = (double *)malloc(ny_fine*sizeof(double));
    grid->lont   = (double *)malloc(nx_fine*ny_fine*sizeof(double));
    grid->latt   = (double *)malloc(nx_fine*ny_fine*sizeof(double));
    grid->xt     = (double *)malloc(nx_fine*ny_fine*sizeof(double));
    grid->yt     = (double *)malloc(nx_fine*ny_fine*sizeof(double));
    grid->zt     = (double *)malloc(nx_fine*ny_fine*sizeof(double));
    grid->vlon_t = (double *)malloc(3*nx_fine*ny_fine*sizeof(double));
    grid->vlat_t = (double *)malloc(3*nx_fine*ny_fine*sizeof(double));

    dlon = lon_range/nx_fine;
    for(i=0; i<nx_fine; i++) {
      lon_fine = (lonbegin + (i + 0.5)*dlon)*D2R;
      for(j=0; j<ny_fine; j++) grid->lont[j*nx_fine+i] = lon_fine;
    }
    if(center_y) {
      dlat=lat_range/ny_fine;
      for(j=0; j<ny_fine; j++)grid->latt1D_fine[j] = (latbegin+(j+0.5)*dlat)*D2R;
    }
    else {
      dlat = lat_range/(ny_fine-1);
      for(j=0; j<ny_fine; j++) grid->latt1D_fine[j] = (latbegin+j*dlat)*D2R;

    }
    for(j=0; j<ny_fine; j++) for(i=0; i<nx_fine; i++) {
        grid->latt[j*nx_fine+i] = grid->latt1D_fine[j];
      }
    /* get the cartesian coordinates */
    latlon2xyz(nx_fine*ny_fine, grid->lont, grid->latt, grid->xt, grid->yt, grid->zt);

    unit_vect_latlon(nx_fine*ny_fine, grid->lont, grid->latt, grid->vlon_t, grid->vlat_t);

  }

  grid->lonc  = (double *) malloc((nxc+1)*(nyc+1)*sizeof(double));
  grid->latc  = (double *) malloc((nxc+1)*(nyc+1)*sizeof(double));
  for(j=0; j<=nyc; j++) {
    jj = j + grid->jsc;
    for(i=0; i<=nxc; i++) {
      ii = i + grid->isc;
      grid->lonc[j*(nxc+1)+i] = grid->lonc1D[ii];
      grid->latc[j*(nxc+1)+i] = grid->latc1D[jj];
    }
  }
  if(opcode & VECTOR) { /* no rotation is needed for regular lat-lon grid. */
    grid->rotate = 0;
  }

}; /* get_output_grid_by_size */


void setup_vertical_interp_acc(VGrid_config *vgrid_in, VGrid_config *vgrid_out)
{
  int nk1, nk2, kstart, kend, k;

  nk1 = vgrid_in->nz;
  nk2 = vgrid_out->nz;


  for(kstart=0; kstart<nk2; kstart++) {
    if(vgrid_out->z[kstart] >= vgrid_in->z[0]) break;
  }
  for(kend=nk2-1; kend>=0; kend--) {
    if(vgrid_out->z[kend] <= vgrid_in->z[nk1-1]) break;
  }


  if(kstart >0 && mpp_pe()==mpp_root_pe()) {
    printf("NOTE from fregrid_util: the value from level 0 to level %d will be set to the value at the shallowest source levle.\n", kstart-1);
  }
  if(kend <nk2-1 && mpp_pe()==mpp_root_pe()) {
    printf("NOTE from fregrid_util: the value from level %d to level %d will be set to the value at the deepest source levle.\n", kend+1, nk2-1);
  }
  vgrid_out->kstart = kstart;
  vgrid_out->kend   = kend;
  vgrid_out->need_interp = 1;
  if(nk1 == nk2 ){
    for(k=0; k<nk1; k++) {
      if(fabs(vgrid_out->z[k]-vgrid_in->z[k]) > EPSLN10 ) break;
    }
    if(k==nk1) vgrid_out->need_interp = 0;
  }
}

void do_vertical_interp_acc(VGrid_config *vgrid_in, VGrid_config *vgrid_out, Grid_config *grid_out, Field_config *field, int varid)
{
  int nk1, nk2, nx, ny, kstart, kend, i, k;
  double *tmp;

  if(vgrid_out->need_interp && field->var[varid].has_zaxis ) {
    nk1 = vgrid_in->nz;
    nk2 = vgrid_out->nz;
    nx  = grid_out->nx;
    ny  = grid_out->ny;
    tmp = (double *)malloc(nx*ny*nk1*sizeof(double));
    for(i=0; i<nx*ny*nk1; i++) tmp[i] = field->data[i];
    if(nk1 != nk2 ) {
      free(field->data);
      field->data =  (double *)malloc(nx*ny*nk2*sizeof(double));
    }

    kstart = vgrid_out->kstart;
    kend   = vgrid_out->kend;
    for(k=0; k<kstart; k++) {
      for(i=0; i<nx*ny; i++) field->data[k*nx*ny+i] = tmp[i];
    }
    for(k=kend; k<nk2; k++) {
      for(i=0; i<nx*ny; i++) field->data[k*nx*ny+i] = tmp[(nk1-1)*nx*ny+i];
    }
    nk2 = kend - kstart + 1;
    linear_vertical_interp(nx, ny, nk1, nk2, vgrid_in->z, vgrid_out->z+kstart, tmp, field->data+kstart*nx*ny);
    free(tmp);
  }

}
