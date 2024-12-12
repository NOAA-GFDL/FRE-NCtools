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

// This test tests the function create_xgrid_2dx2d_order2_gpu for a simple
// case where the input grid is identical to the output grid.  The exchange
// grid is identical to the input/output grid.  Only the parent indices corresponding
// to each exchange cell are checked.

#include <stdlib.h>
#include <stdio.h>
#include <openacc.h>
#include "globals_gpu.h"
#include "interp_utils_gpu.h"
#include "conserve_interp_gpu.h"
#include "create_xgrid_utils_gpu.h"
#include "create_xgrid_gpu.h"

void generate_input_is_same_as_output_grids(Grid_config *input_grid, Grid_config *output_grid,
                                            Interp_per_input_tile *interp_answers);
int run_tests(Grid_config *input_grid,  Grid_config *output_grid,
              Interp_per_input_tile *interp);
void check_answers(const Interp_per_input_tile *interp_answers, const Interp_per_input_tile *interp);
void check_ianswers(const int n, const int *answers, const int *checkme);

//TODO:  add more complicated tests
int main()
{

  Grid_config input_grid, output_grid;
  Interp_per_input_tile interp_answers, interp;

  printf("CHECKING CASE 1:  Input grid = Output_grid\n");

  // get input_grid and output_grid
  generate_input_is_same_as_output_grids(&input_grid, &output_grid, &interp_answers);

  // run create_xgrid order 1
  run_tests(&input_grid, &output_grid, &interp);

  // check answers
  check_answers(&interp_answers, &interp);

  //test successful
  return 0;

}

void generate_input_is_same_as_output_grids(Grid_config *input_grid, Grid_config *output_grid,
                                            Interp_per_input_tile *interp_answers)
{

  int nlon_cells = 360-1; //[ 0, 1, ... 360]
  int nlat_cells = 3-1;   //[ -30, 0, 30]
  int ncells = nlon_cells * nlat_cells;
  int ngridpts = (nlon_cells+1) * (nlat_cells+1);
  double dlat = 30.0;
  double dlon = 1.0;

  int i=0;

  input_grid->nxc  = nlon_cells  ; input_grid->nyc  = nlat_cells;
  output_grid->nxc = nlon_cells  ; output_grid->nyc = nlat_cells;

  input_grid->latc = (double *)malloc(ngridpts*sizeof(double));
  input_grid->lonc = (double *)malloc(ngridpts*sizeof(double));
  output_grid->latc = (double *)malloc(ngridpts*sizeof(double));
  output_grid->lonc = (double *)malloc(ngridpts*sizeof(double));

  for( int ilat=0 ; ilat<nlat_cells+1 ; ilat++){
    double latitude = ilat * dlat * D2R;
    for( int ilon=0 ; ilon<nlon_cells+1 ; ilon++) {
      input_grid->latc[i] = latitude;
      output_grid->latc[i] = latitude;
      input_grid->lonc[i] = (ilon * dlon)*D2R;
      output_grid->lonc[i] = (ilon * dlon)*D2R;
      i++;
    }
  }

  //answers
  interp_answers->nxcells = ncells;
  interp_answers->input_parent_cell_index = (int *)malloc(ncells*sizeof(int));
  interp_answers->output_parent_cell_index = (int *)malloc(ncells*sizeof(int));
  //interp_answers
  for(int ij1=0 ; ij1<ncells ; ij1++) {
    interp_answers->input_parent_cell_index[ij1] = ij1;
    interp_answers->output_parent_cell_index[ij1] = ij1;
  }

}


void check_answers(const Interp_per_input_tile *interp_answers, const Interp_per_input_tile *interp)
{

  //reset interp on host to reassure device data is correct
  for(int i=0 ; i<interp_answers->nxcells; i++) {
    interp->input_parent_cell_index[i] = -99;
    interp->output_parent_cell_index[i] = -99;
  }

#pragma acc update host( interp->nxcells )
  int nxcells = interp->nxcells;

#pragma acc exit data copyout( interp->input_parent_cell_index[:nxcells],\
                               interp->output_parent_cell_index[:nxcells])

  printf("checking nxcells\n");
  check_ianswers(1, &(interp_answers->nxcells), &interp->nxcells);

  printf("checking input_parent_cell_index\n");
  check_ianswers(nxcells, interp_answers->input_parent_cell_index, interp->input_parent_cell_index);

  printf("checking output_parent_cell_index\n");
  check_ianswers(nxcells, interp_answers->output_parent_cell_index, interp->output_parent_cell_index);

}


int run_tests(Grid_config *input_grid, Grid_config *output_grid, Interp_per_input_tile *interp)
{

  int nlon_input_cells  = input_grid->nxc;
  int nlat_input_cells  = input_grid->nyc;
  int nlon_output_cells = output_grid->nxc;
  int nlat_output_cells = output_grid->nyc;
  int ncells_input  = nlon_input_cells * nlat_input_cells;
  int ncells_output = nlon_output_cells * nlat_output_cells;
  int ngridpts_input = (nlon_input_cells+1)*(nlat_input_cells+1);
  int ngridpts_output = (nlon_output_cells+1)*(nlat_output_cells+1);

  int upbound_nxcells, nxcells, jlat_overlap_starts, jlat_overlap_ends;
  int *approx_nxcells_per_ij1, *ij2_start, *ij2_end;
  double *input_grid_mask;
  Grid_cells_struct_config output_grid_cells;

  //copy grid to device
  copy_grid_to_device_gpu(ngridpts_input, input_grid->latc, input_grid->lonc);
  copy_grid_to_device_gpu(ngridpts_output, output_grid->latc, output_grid->lonc);

  //get mask to skip input cells in creating interp
  get_input_grid_mask_gpu(ncells_input, &input_grid_mask);

  //get output grid cell info
  get_grid_cell_struct_gpu(nlon_output_cells, nlat_output_cells, output_grid, &output_grid_cells);

  //get bounding index
  get_bounding_indices_gpu(nlon_output_cells, nlat_output_cells, nlon_input_cells, nlat_input_cells,
                           output_grid->latc, input_grid->latc, &jlat_overlap_starts, &jlat_overlap_ends);

  //malloc and create arrays
  create_upbound_nxcells_arrays_on_device_gpu(ncells_input, &approx_nxcells_per_ij1, &ij2_start, &ij2_end);

  upbound_nxcells = get_upbound_nxcells_2dx2d_gpu(input_grid->nxc, input_grid->nyc, output_grid->nxc, output_grid->nyc,
                                                  jlat_overlap_starts, jlat_overlap_ends,
                                                  input_grid->lonc, input_grid->latc,
                                                  output_grid->lonc, output_grid->latc,
                                                  input_grid_mask, &output_grid_cells,
                                                  approx_nxcells_per_ij1, ij2_start, ij2_end);

  nxcells = create_xgrid_2dx2d_order1_gpu( nlon_input_cells, nlat_input_cells,
                                           nlon_output_cells, nlat_output_cells,
                                           jlat_overlap_starts, jlat_overlap_ends,
                                           input_grid->lonc, input_grid->latc,
                                           output_grid->lonc, output_grid->latc,
                                           upbound_nxcells, input_grid_mask, &output_grid_cells,
                                           approx_nxcells_per_ij1, ij2_start, ij2_end,
                                           interp);

  free_grid_cell_struct_gpu(ncells_output, &output_grid_cells);
  free_upbound_nxcells_arrays_gpu(ncells_input, &approx_nxcells_per_ij1, &ij2_start, &ij2_end);
  free_input_grid_mask_gpu(ncells_input, &input_grid_mask);

  return 0;

}


void check_ianswers( const int n, const int *answers, const int *checkme)
{

  for(int i=0 ; i<n ; i++){
    if( answers[i]!=checkme[i] ) {
      printf("EXPECTED %d but got %d\n for element %d\n", answers[i], checkme[i], i);
      exit(1);
    }
  }

}
