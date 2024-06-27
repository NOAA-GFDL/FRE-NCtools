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
#include <openacc.h>
#include "general_utils_acc.h"
#include "create_xgrid_acc.h"
#include "create_xgrid_utils_acc.h"
#include "globals_acc.h"

/*******************************************************************************
void get_upbound_nxcells_2dx2d_acc
This function computes the upperbound to nxgrid.  This upper bound will be used
to malloc arrays used in create_xgrid
*******************************************************************************/
int get_upbound_nxcells_2dx2d_acc(const int nlon_input_cells,  const int nlat_input_cells,
                                  const int nlon_output_cells, const int nlat_output_cells,
                                  const int jlat_overlap_starts, const int jlat_overlap_ends,
                                  const double *input_grid_lon, const double *input_grid_lat,
                                  const double *output_grid_lon, const double *output_grid_lat,
                                  const double *skip_input_cells,
                                  const Grid_cells_struct_config *output_grid_cells,
                                  int *approx_xcells_per_ij1, int *ij2_start, int *ij2_end)
{

  int input_grid_ncells  = nlon_input_cells*nlat_input_cells;
  int output_grid_ncells = nlon_output_cells*nlat_output_cells;
  int input_grid_npts    = (nlon_input_cells+1)*(nlat_input_cells+1);
  int output_grid_npts   = (nlon_output_cells+1)*(nlat_output_cells+1);

  int ij1_start = jlat_overlap_starts*nlon_input_cells;
  int ij1_end = (jlat_overlap_ends+1)*nlon_input_cells;
  int upbound_nxcells=0;

#pragma acc data present(output_grid_lon[:output_grid_npts],  \
                         output_grid_lat[:output_grid_npts],  \
                         input_grid_lon[:input_grid_npts],   \
                         input_grid_lat[:input_grid_npts],          \
                         output_grid_cells[:1],                     \
                         approx_xcells_per_ij1[:input_grid_ncells], \
                         ij2_start[:input_grid_ncells],             \
                         ij2_end[:input_grid_ncells],               \
                         skip_input_cells[:input_grid_ncells])
#pragma acc data copyin(input_grid_ncells, \
                        output_grid_ncells)
#pragma acc data copy(upbound_nxcells)
#pragma acc parallel loop independent reduction(+:upbound_nxcells)
  for( int ij1=ij1_start ; ij1<ij1_end ; ij1++) {
    if( skip_input_cells[ij1] > MASK_THRESH ) {

      int i_approx_xcells_per_ij1=0;
      int ij2_min=output_grid_ncells, ij2_max=0;
      double input_cell_lon_vertices[MV], input_cell_lat_vertices[MV];

      get_cell_vertices_acc(ij1, nlon_input_cells, input_grid_lon, input_grid_lat,
                            input_cell_lon_vertices, input_cell_lat_vertices);

      double input_cell_lat_min = minval_double_acc(4, input_cell_lat_vertices);
      double input_cell_lat_max = maxval_double_acc(4, input_cell_lat_vertices);
      int nvertices = fix_lon_acc(input_cell_lon_vertices, input_cell_lat_vertices, 4, M_PI);
      double input_cell_lon_min = minval_double_acc(nvertices, input_cell_lon_vertices);
      double input_cell_lon_max = maxval_double_acc(nvertices, input_cell_lon_vertices);
      double input_cell_lon_cent = avgval_double_acc(nvertices, input_cell_lon_vertices);

      approx_xcells_per_ij1[ij1]=0;

#pragma acc loop independent reduction(+:upbound_nxcells) reduction(+:i_approx_xcells_per_ij1) \
                             reduction(min:ij2_min) reduction(max:ij2_max)

      for(int ij2=0; ij2<output_grid_ncells; ij2++) {

        double dlon_cent, output_cell_lon_min, output_cell_lon_max;
        double rotate=0.0;

        if(output_grid_cells->lat_min[ij2] >= input_cell_lat_max) continue;
        if(output_grid_cells->lat_max[ij2] <= input_cell_lat_min) continue;

        dlon_cent = output_grid_cells->lon_cent[ij2] - input_cell_lon_cent;
        if(dlon_cent < -M_PI) rotate = +TPI;
        if(dlon_cent > M_PI)  rotate = -TPI;

        // adjust according to input_grid_lon_cent
        // TODO: breakup grid into quadrants to avoid if statements?
        output_cell_lon_min = output_grid_cells->lon_min[ij2] + rotate;
        output_cell_lon_max = output_grid_cells->lon_max[ij2] + rotate;

        //output_cell_lon should in the same range as input_cell_lon after lon_fix,
        // so no need to consider cyclic condition
        if(output_cell_lon_min >= input_cell_lon_max ) continue;
        if(output_cell_lon_max <= input_cell_lon_min ) continue;

        //Note, the check for AREA_RATIO_THRESH has been removed
        //Thus, the computed value of upbound_nxcells will be equal to or greater than nxgrid
        i_approx_xcells_per_ij1++;
        upbound_nxcells++;
        ij2_min = min(ij2_min, ij2);
        ij2_max = max(ij2_max, ij2);

      } //ij2
      approx_xcells_per_ij1[ij1] = i_approx_xcells_per_ij1;
      ij2_start[ij1] = ij2_min ;
      ij2_end[ij1]   = ij2_max;

    } //mask
  } //ij1

  return upbound_nxcells;

}

/*******************************************************************************
  void create_xgrid_2DX2D_order1
  This routine generate exchange grids between two grids for the first order
  conservative interpolation. nlon_input_cells,ninput_grid_lat,nlon_output_cells,nlat_output_cells are the size of the grid cell
  and input_grid_lon,input_grid_lat, output_grid_lon,output_grid_lat are geographic grid location of grid cell bounds.
*******************************************************************************/
int create_xgrid_2dx2d_order1_acc(const int nlon_input_cells,  const int nlat_input_cells,
                                  const int nlon_output_cells, const int nlat_output_cells,
                                  const int jlat_overlap_starts, const int jlat_overlap_ends,
                                  const double *input_grid_lon,  const double *input_grid_lat,
                                  const double *output_grid_lon, const double *output_grid_lat,
                                  const int upbound_nxcells, const double *mask_input_grid,
                                  const Grid_cells_struct_config *output_grid_cells,
                                  int *approx_nxcells_per_ij1, int *ij2_start, int *ij2_end,
                                  Interp_per_input_tile *interp_for_input_tile)
{

  if(upbound_nxcells<1) return 0;

  int nxcells=0;

  int input_grid_ncells  = nlon_input_cells*nlat_input_cells;
  int output_grid_ncells = nlon_output_cells*nlat_output_cells;
  int input_grid_npts    = (nlon_input_cells+1)*(nlat_input_cells+1);
  int output_grid_npts   = (nlon_output_cells+1)*(nlat_output_cells+1);

  int ij1_start = jlat_overlap_starts*nlon_input_cells;
  int ij1_end   = (jlat_overlap_ends+1)*nlon_input_cells;

  double xcell_dclon=-99.99, xcell_dclat=-99.99;

  int *parent_input_index  = NULL ; parent_input_index  = (int *)malloc(upbound_nxcells*sizeof(int));
  int *parent_output_index = NULL ; parent_output_index = (int *)malloc(upbound_nxcells*sizeof(int));
  int *nxcells_per_ij1     = NULL ; nxcells_per_ij1 = (int *)malloc(input_grid_ncells*sizeof(int));
  double *store_xcell_area = NULL ; store_xcell_area =  (double *)malloc(upbound_nxcells*sizeof(double));

#pragma acc enter data create(parent_input_index[:upbound_nxcells],  \
                              parent_output_index[:upbound_nxcells], \
                              store_xcell_area[:upbound_nxcells],    \
                              nxcells_per_ij1[:input_grid_ncells])

#pragma acc data present(output_grid_lon[:output_grid_npts],         \
                         output_grid_lat[:output_grid_npts],         \
                         input_grid_lon[:input_grid_npts],           \
                         input_grid_lat[:input_grid_npts],           \
                         output_grid_cells[:1],                      \
                         approx_nxcells_per_ij1[:input_grid_ncells], \
                         ij2_start[:input_grid_ncells],              \
                         ij2_end[:input_grid_ncells],                \
                         mask_input_grid[:input_grid_ncells],        \
                         nxcells_per_ij1[:input_grid_ncells],        \
                         parent_input_index[:upbound_nxcells],       \
                         parent_output_index[:upbound_nxcells],      \
                         store_xcell_area[:upbound_nxcells])         \
  copyin(input_grid_ncells, output_grid_ncells)                      \
  copy(nxcells)
#pragma acc parallel loop reduction(+:nxcells)
  for(int ij1=ij1_start; ij1<ij1_end; ij1++) {
    if(mask_input_grid[ij1] > MASK_THRESH)  {

      double input_cell_lon_vertices[MV], input_cell_lat_vertices[MV];
      int approx_nxcells_b4_ij1=0, ixcell=0;

      get_cell_vertices_acc(ij1, nlon_input_cells, input_grid_lon, input_grid_lat,
                            input_cell_lon_vertices, input_cell_lat_vertices);
      double input_cell_lat_min = minval_double_acc(4, input_cell_lat_vertices);
      double input_cell_lat_max = maxval_double_acc(4, input_cell_lat_vertices);
      int nvertices1 = fix_lon_acc(input_cell_lon_vertices, input_cell_lat_vertices, 4, M_PI);
      double input_cell_lon_min = minval_double_acc(nvertices1, input_cell_lon_vertices);
      double input_cell_lon_max = maxval_double_acc(nvertices1, input_cell_lon_vertices);
      double input_cell_lon_cent = avgval_double_acc(nvertices1, input_cell_lon_vertices);
      double input_cell_area = poly_area_acc(input_cell_lon_vertices, input_cell_lat_vertices, nvertices1);

#pragma acc loop seq
      for(int i=1; i<=ij1 ; i++) approx_nxcells_b4_ij1 += approx_nxcells_per_ij1[i-1];
      nxcells_per_ij1[ij1]=0;

#pragma acc loop seq reduction(+:ixcell)
      for(int ij2=ij2_start[ij1]; ij2<=ij2_end[ij1]; ij2++) {

        int nvertices2, xvertices=1;
        double dlon_cent, output_cell_lon_min, output_cell_lon_max, output_cell_area;
        double output_cell_lon_vertices[MAX_V], output_cell_lat_vertices[MAX_V];
        double xcell_lon_vertices[MV], xcell_lat_vertices[MV];

        double rotate=0.0;

        if(output_grid_cells->lat_min[ij2] >= input_cell_lat_max) continue;
        if(output_grid_cells->lat_max[ij2] <= input_cell_lat_min) continue;

        /* adjust according to input_grid_lon_cent*/
        output_cell_lon_min = output_grid_cells->lon_min[ij2];
        output_cell_lon_max = output_grid_cells->lon_max[ij2];
        nvertices2 = output_grid_cells->nvertices[ij2];
        output_cell_area = output_grid_cells->area[ij2];

        dlon_cent = output_grid_cells->lon_cent[ij2] - input_cell_lon_cent;

        if(dlon_cent < -M_PI) rotate = TPI;
        if(dlon_cent > M_PI)  rotate = -TPI;

        output_cell_lon_min += rotate;
        output_cell_lon_max += rotate;
        for (int l=0; l<nvertices2; l++) {
          output_cell_lon_vertices[l] = output_grid_cells->lon_vertices[ij2][l] + rotate;
          output_cell_lat_vertices[l] = output_grid_cells->lat_vertices[ij2][l];
        }

        //output_cell_lon should in the same range as input_cell_lon after lon_fix,
        // so no need to consider cyclic condition
        if(output_cell_lon_min >= input_cell_lon_max ) continue;
        if(output_cell_lon_max <= input_cell_lon_min ) continue;

        if ( (xvertices = clip_2dx2d_acc( input_cell_lon_vertices, input_cell_lat_vertices, nvertices1,
                                          output_cell_lon_vertices, output_cell_lat_vertices, nvertices2,
                                          xcell_lon_vertices, xcell_lat_vertices)) > 0 ){
          double xcell_area = poly_area_acc(xcell_lon_vertices, xcell_lat_vertices, xvertices);
          if( xcell_area/min(input_cell_area, output_cell_area) > AREA_RATIO_THRESH ) {
            store_xcell_area[approx_nxcells_b4_ij1+ixcell] = xcell_area;
            parent_input_index[approx_nxcells_b4_ij1+ixcell]  = ij1;
            parent_output_index[approx_nxcells_b4_ij1+ixcell] = ij2;
            ixcell++;
          }
        }
      }
      nxcells+=ixcell;
      nxcells_per_ij1[ij1]=ixcell;
    }
  }

  copy_data_to_interp_on_device_acc(nxcells, input_grid_ncells, upbound_nxcells, nxcells_per_ij1,
                                    &xcell_dclon, &xcell_dclat,
                                    approx_nxcells_per_ij1, parent_input_index, parent_output_index,
                                    store_xcell_area, interp_for_input_tile);

#pragma acc exit data delete( parent_input_index[:upbound_nxcells],  \
                              parent_output_index[:upbound_nxcells], \
                              store_xcell_area[:upbound_nxcells],    \
                              nxcells_per_ij1[:input_grid_ncells])

  free(parent_input_index) ; parent_input_index = NULL;
  free(parent_output_index); parent_output_index = NULL;
  free(nxcells_per_ij1)    ; nxcells_per_ij1 = NULL;
  free(store_xcell_area)   ; store_xcell_area = NULL;

  return nxcells;

};/* get_xgrid_2Dx2D_order1 */


/*******************************************************************************
  void create_xgrid_2DX2D_order2
  This routine generate exchange grids between two grids for the second order
*******************************************************************************/
int create_xgrid_2dx2d_order2_acc(const int nlon_input_cells,  const int nlat_input_cells,
                                  const int nlon_output_cells, const int nlat_output_cells,
                                  const int jlat_overlap_starts, const int jlat_overlap_ends,
                                  const double *input_grid_lon,  const double *input_grid_lat,
                                  const double *output_grid_lon, const double *output_grid_lat,
                                  const int upbound_nxcells, const double *mask_input_grid,
                                  const Grid_cells_struct_config *output_grid_cells,
                                  int *approx_nxcells_per_ij1, int *ij2_start, int *ij2_end,
                                  Interp_per_input_tile *interp_for_input_tile, double *readin_input_area)
{

  if(upbound_nxcells<1) return 0;

  int nxcells=0;

  int input_grid_ncells  = nlon_input_cells*nlat_input_cells;
  int output_grid_ncells = nlon_output_cells*nlat_output_cells;
  int input_grid_npts    = (nlon_input_cells+1)*(nlat_input_cells+1);
  int output_grid_npts   = (nlon_output_cells+1)*(nlat_output_cells+1);

  int ij1_start = jlat_overlap_starts*nlon_input_cells;
  int ij1_end = (jlat_overlap_ends+1)*nlon_input_cells;

  int *parent_input_index=NULL   ; parent_input_index = (int *)malloc(upbound_nxcells*sizeof(int));
  int *parent_output_index=NULL  ; parent_output_index = (int *)malloc(upbound_nxcells*sizeof(int));
  double *store_xcell_area=NULL  ; store_xcell_area  = (double *)malloc(upbound_nxcells*sizeof(double));
  double *store_xcell_dclon=NULL ; store_xcell_dclon = (double *)malloc(upbound_nxcells*sizeof(double));
  double *store_xcell_dclat=NULL ; store_xcell_dclat = (double *)malloc(upbound_nxcells*sizeof(double));

  int *nxcells_per_ij1=NULL ; nxcells_per_ij1 = (int *)malloc(input_grid_ncells*sizeof(int));
  double *summed_input_area=NULL; summed_input_area = (double *)malloc(input_grid_ncells*sizeof(double));
  double *summed_input_clat=NULL; summed_input_clat = (double *)malloc(input_grid_ncells*sizeof(double));
  double *summed_input_clon=NULL; summed_input_clon = (double *)malloc(input_grid_ncells*sizeof(double));

#pragma acc enter data create(parent_input_index[:upbound_nxcells],  \
                              parent_output_index[:upbound_nxcells], \
                              store_xcell_area[:upbound_nxcells],    \
                              nxcells_per_ij1[:input_grid_ncells],   \
                              store_xcell_dclon[:upbound_nxcells],   \
                              store_xcell_dclat[:upbound_nxcells],   \
                              summed_input_area[:input_grid_ncells], \
                              summed_input_clon[:input_grid_ncells], \
                              summed_input_clat[:input_grid_ncells])

#pragma acc data present(output_grid_lon[:output_grid_npts],         \
                         output_grid_lat[:output_grid_npts],         \
                         input_grid_lon[:input_grid_npts],           \
                         input_grid_lat[:input_grid_npts],           \
                         output_grid_cells[:1],                      \
                         approx_nxcells_per_ij1[:input_grid_ncells], \
                         ij2_start[:input_grid_ncells],              \
                         ij2_end[:input_grid_ncells],                \
                         mask_input_grid[:input_grid_ncells],        \
                         nxcells_per_ij1[:input_grid_ncells],        \
                         parent_input_index[:upbound_nxcells],       \
                         parent_output_index[:upbound_nxcells],      \
                         store_xcell_area[:upbound_nxcells],         \
                         store_xcell_dclon[:upbound_nxcells],        \
                         store_xcell_dclat[:upbound_nxcells])        \
  copyin(input_grid_ncells, output_grid_ncells) copy(nxcells)
#pragma acc parallel loop reduction(+:nxcells)
  for(int ij1=ij1_start; ij1<ij1_end; ij1++) {
    if(mask_input_grid[ij1] > MASK_THRESH)  {

      double input_cell_lon_vertices[MV], input_cell_lat_vertices[MV];
      double summed_input_area_ij1=0.;
      double summed_input_clon_ij1=0.;
      double summed_input_clat_ij1=0.;
      int approx_nxcells_b4_ij1=0, ixcell=0;

      get_cell_vertices_acc(ij1, nlon_input_cells, input_grid_lon, input_grid_lat,
                            input_cell_lon_vertices, input_cell_lat_vertices);
      double input_cell_lat_min = minval_double_acc(4, input_cell_lat_vertices);
      double input_cell_lat_max = maxval_double_acc(4, input_cell_lat_vertices);
      int nvertices1 = fix_lon_acc(input_cell_lon_vertices, input_cell_lat_vertices, 4, M_PI);
      double input_cell_lon_min = minval_double_acc(nvertices1, input_cell_lon_vertices);
      double input_cell_lon_max = maxval_double_acc(nvertices1, input_cell_lon_vertices);
      double input_cell_lon_cent = avgval_double_acc(nvertices1, input_cell_lon_vertices);
      double input_cell_area = poly_area_acc(input_cell_lon_vertices, input_cell_lat_vertices, nvertices1);

#pragma acc loop seq
      for(int i=1; i<=ij1 ; i++) approx_nxcells_b4_ij1 += approx_nxcells_per_ij1[i-1];

#pragma acc loop seq reduction(+:ixcell)                \
                     reduction(+:summed_input_area_ij1) \
                     reduction(+:summed_input_clon_ij1) \
                     reduction(+:summed_input_clat_ij1)
      for(int ij2=ij2_start[ij1]; ij2<=ij2_end[ij1]; ij2++) {

        int nvertices2, xvertices=1;
        double dlon_cent, output_cell_lon_min, output_cell_lon_max;
        double output_cell_lon_vertices[MAX_V], output_cell_lat_vertices[MAX_V];
        double output_cell_area;
        double xcell_lon_vertices[MV], xcell_lat_vertices[MV];

        double rotate=0.0;

        if(output_grid_cells->lat_min[ij2] >= input_cell_lat_max) continue;
        if(output_grid_cells->lat_max[ij2] <= input_cell_lat_min) continue;

        dlon_cent = output_grid_cells->lon_cent[ij2] - input_cell_lon_cent;
        if(dlon_cent < -M_PI) rotate = TPI;
        if(dlon_cent > M_PI)  rotate = -TPI;

        /* adjust according to input_grid_lon_cent*/
        output_cell_lon_min = output_grid_cells->lon_min[ij2] + rotate;
        output_cell_lon_max = output_grid_cells->lon_max[ij2] + rotate;
        nvertices2 = output_grid_cells->nvertices[ij2];
        output_cell_area = output_grid_cells->area[ij2];

        for (int l=0; l<nvertices2; l++) {
          output_cell_lon_vertices[l] = output_grid_cells->lon_vertices[ij2][l] + rotate;
          output_cell_lat_vertices[l] = output_grid_cells->lat_vertices[ij2][l];
        }

        //output_cell_lon should in the same range as input_cell_lon after lon_fix,
        // so no need to consider cyclic condition
        if(output_cell_lon_min >= input_cell_lon_max ) continue;
        if(output_cell_lon_max <= input_cell_lon_min ) continue;

        if ( (xvertices = clip_2dx2d_acc( input_cell_lon_vertices, input_cell_lat_vertices, nvertices1,
                                          output_cell_lon_vertices, output_cell_lat_vertices, nvertices2,
                                          xcell_lon_vertices, xcell_lat_vertices)) > 0 ){
          double xcell_area = poly_area_acc(xcell_lon_vertices, xcell_lat_vertices, xvertices);
          if( xcell_area/min(input_cell_area, output_cell_area) > AREA_RATIO_THRESH ) {
            double xcell_clon, xcell_clat;
            store_xcell_area[approx_nxcells_b4_ij1+ixcell] = xcell_area;
            parent_input_index[approx_nxcells_b4_ij1+ixcell]  = ij1;
            parent_output_index[approx_nxcells_b4_ij1+ixcell] = ij2;
            poly_ctrlon_acc(xcell_lon_vertices, xcell_lat_vertices, xvertices, input_cell_lon_cent, &xcell_clon);
            poly_ctrlat_acc(xcell_lon_vertices, xcell_lat_vertices, xvertices, &xcell_clat);
            store_xcell_dclon[approx_nxcells_b4_ij1+ixcell] = xcell_clon/xcell_area;
            store_xcell_dclat[approx_nxcells_b4_ij1+ixcell] = xcell_clat/xcell_area;
            summed_input_area_ij1 += xcell_area;
            summed_input_clon_ij1 += xcell_clon;
            summed_input_clat_ij1 += xcell_clat;
            ixcell++;
          }
        }
      }

      nxcells+=ixcell;
      nxcells_per_ij1[ij1]=ixcell;
      summed_input_area[ij1] = summed_input_area_ij1;
      summed_input_clon[ij1] = summed_input_clon_ij1;
      summed_input_clat[ij1] = summed_input_clat_ij1;

    }
  }

  copy_data_to_interp_on_device_acc(nxcells, input_grid_ncells, upbound_nxcells, nxcells_per_ij1,
                                    store_xcell_dclon, store_xcell_dclat, approx_nxcells_per_ij1, parent_input_index,
                                    parent_output_index, store_xcell_area, interp_for_input_tile);

#pragma acc parallel loop present(interp_for_input_tile->dcentroid_lat[:nxcells], \
                                  interp_for_input_tile->dcentroid_lat[:nxcells], \
                                  input_grid_lon[:input_grid_ncells],   \
                                  input_grid_lat[:input_grid_ncells],   \
                                  summed_input_area[:input_grid_ncells], \
                                  summed_input_clon[:input_grid_ncells], \
                                  summed_input_clat[:input_grid_ncells], \
                                  interp_for_input_tile->input_parent_cell_index[:nxcells]) \
                           copyin(readin_input_area[:input_grid_ncells])
    for(int ix=0 ; ix<nxcells ; ix++){
      int ij1 = interp_for_input_tile->input_parent_cell_index[ix];
      double input_area = summed_input_area[ij1];
      double input_clon = summed_input_clon[ij1];
      double input_clat = summed_input_clat[ij1];
      double readin_area = readin_input_area[ij1];
      if(fabs(input_area - readin_area)/readin_area > AREA_RATIO) {
        double x[4], y[4], input_cell_lon_cent;
        get_cell_vertices_acc(ij1, nlon_input_cells, input_grid_lon, input_grid_lat, x, y);
        int n = fix_lon_acc(x, y, 4, M_PI);
        input_cell_lon_cent = avgval_double_acc(n, x);
        poly_ctrlon_acc(x, y, n, input_cell_lon_cent, &input_clon);
        poly_ctrlat_acc(x, y, n, &input_clat);
        input_area = readin_area;
      }
      interp_for_input_tile->dcentroid_lon[ix] -= input_clon/input_area;
      interp_for_input_tile->dcentroid_lat[ix] -= input_clat/input_area;
    }

#pragma acc exit data delete( parent_input_index[:upbound_nxcells],     \
                              parent_output_index[:upbound_nxcells],    \
                              store_xcell_area[:upbound_nxcells],       \
                              nxcells_per_ij1[:input_grid_ncells],      \
                              store_xcell_dclon[:upbound_nxcells],      \
                              store_xcell_dclat[:upbound_nxcells],      \
                              summed_input_area[:input_grid_ncells],    \
                              summed_input_clon[:input_grid_ncells],    \
                              summed_input_clat[:input_grid_ncells])

    free(parent_input_index)    ; parent_input_index = NULL;
    free(parent_output_index)   ; parent_output_index = NULL;
    free(store_xcell_area)      ; store_xcell_area = NULL;
    free(nxcells_per_ij1)       ; nxcells_per_ij1 = NULL;
    free(store_xcell_dclon)     ; store_xcell_dclon = NULL;
    free(store_xcell_dclat)     ; store_xcell_dclat = NULL;
    free(summed_input_area); summed_input_area = NULL;
    free(summed_input_clon); summed_input_clon = NULL;
    free(summed_input_clat); summed_input_clat = NULL;

  return nxcells;

};/* get_xgrid_2Dx2D_order2 */

int create_xgrid_great_circle_acc(const int *nlon_input_cells, const int *nlat_input_cells,
                                  const int *nlon_output_cells, const int *nlat_output_cells,
                                  const double *input_grid_lon, const double *input_grid_lat,
                                  const double *output_grid_lon, const double *output_grid_lat,
                                  const double *mask_input_grid,
                                  int *i_in, int *j_in, int *i_out, int *j_out,
                                  double *xgrid_area, double *xgrid_clon, double *xgrid_clat)
{

  int nx1, nx2, ny1, ny2, nx1p, nx2p, ny1p, ny2p, nxgrid, n1_in, n2_in;
  int n0, n1, n2, n3, i1, j1, i2, j2, l, n;
  double x1_in[MV], y1_in[MV], z1_in[MV];
  double x2_in[MV], y2_in[MV], z2_in[MV];
  double x_out[MV], y_out[MV], z_out[MV];
  double *x1=NULL, *y1=NULL, *z1=NULL;
  double *x2=NULL, *y2=NULL, *z2=NULL;

  double xctrlon, xctrlat;
  double *area1, *area2, min_area;

  nx1 = *nlon_input_cells;
  ny1 = *nlat_input_cells;
  nx2 = *nlon_output_cells;
  ny2 = *nlat_output_cells;
  nxgrid = 0;
  nx1p = nx1 + 1;
  nx2p = nx2 + 1;
  ny1p = ny1 + 1;
  ny2p = ny2 + 1;

  /* first convert lon-lat to cartesian coordinates */
  x1 = (double *)malloc(nx1p*ny1p*sizeof(double));
  y1 = (double *)malloc(nx1p*ny1p*sizeof(double));
  z1 = (double *)malloc(nx1p*ny1p*sizeof(double));
  x2 = (double *)malloc(nx2p*ny2p*sizeof(double));
  y2 = (double *)malloc(nx2p*ny2p*sizeof(double));
  z2 = (double *)malloc(nx2p*ny2p*sizeof(double));

  latlon2xyz_acc(nx1p*ny1p, input_grid_lon, input_grid_lat, x1, y1, z1);
  latlon2xyz_acc(nx2p*ny2p, output_grid_lon, output_grid_lat, x2, y2, z2);

  area1  = (double *)malloc(nx1*ny1*sizeof(double));
  area2 = (double *)malloc(nx2*ny2*sizeof(double));
  get_grid_great_circle_area_acc(nlon_input_cells, nlat_input_cells, input_grid_lon, input_grid_lat, area1);
  get_grid_great_circle_area_acc(nlon_output_cells, nlat_output_cells, output_grid_lon, output_grid_lat, area2);
  n1_in = 4;
  n2_in = 4;

  for(int j1=0; j1<ny1; j1++) for(int i1=0; i1<nx1; i1++) {
      if( mask_input_grid[j1*nx1+i1] > MASK_THRESH ) {
          /* clockwise */
          n0 = j1*nx1p+i1;       n1 = (j1+1)*nx1p+i1;
          n2 = (j1+1)*nx1p+i1+1; n3 = j1*nx1p+i1+1;
          x1_in[0] = x1[n0]; y1_in[0] = y1[n0]; z1_in[0] = z1[n0];
          x1_in[1] = x1[n1]; y1_in[1] = y1[n1]; z1_in[1] = z1[n1];
          x1_in[2] = x1[n2]; y1_in[2] = y1[n2]; z1_in[2] = z1[n2];
          x1_in[3] = x1[n3]; y1_in[3] = y1[n3]; z1_in[3] = z1[n3];

          for(j2=0; j2<ny2; j2++) for(i2=0; i2<nx2; i2++) {
              int n_in, n_out;
              double xarea;

              n0 = j2*nx2p+i2;       n1 = (j2+1)*nx2p+i2;
              n2 = (j2+1)*nx2p+i2+1; n3 = j2*nx2p+i2+1;
              x2_in[0] = x2[n0]; y2_in[0] = y2[n0]; z2_in[0] = z2[n0];
              x2_in[1] = x2[n1]; y2_in[1] = y2[n1]; z2_in[1] = z2[n1];
              x2_in[2] = x2[n2]; y2_in[2] = y2[n2]; z2_in[2] = z2[n2];
              x2_in[3] = x2[n3]; y2_in[3] = y2[n3]; z2_in[3] = z2[n3];

              if (  (n_out = clip_2dx2d_great_circle_acc( x1_in, y1_in, z1_in, n1_in, x2_in, y2_in, z2_in, n2_in,
                                                          x_out, y_out, z_out)) > 0) {
                xarea = great_circle_area_acc( n_out, x_out, y_out, z_out ) ;
                min_area = min(area1[j1*nx1+i1], area2[j2*nx2+i2]);
                if( xarea/min_area > AREA_RATIO_THRESH ) {
#ifdef debug_test_create_xgrid
                  printf("(i2,j2)=(%d,%d), (i1,j1)=(%d,%d), xarea=%g\n", i2, j2, i1, j1, xarea);
#endif
                  xgrid_area[nxgrid] = xarea;
                  xgrid_clon[nxgrid] = 0; /*z1l: will be developed very soon */
                  xgrid_clat[nxgrid] = 0;
                  i_in[nxgrid]       = i1;
                  j_in[nxgrid]       = j1;
                  i_out[nxgrid]      = i2;
                  j_out[nxgrid]      = j2;
                  ++nxgrid;
            }
          }
        }
     }
  }


  free(area1);
  free(area2);

  free(x1);
  free(y1);
  free(z1);
  free(x2);
  free(y2);
  free(z2);

  return nxgrid;

};/* create_xgrid_great_circle */
