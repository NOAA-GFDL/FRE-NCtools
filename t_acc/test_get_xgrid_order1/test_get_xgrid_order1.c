#include <stdlib.h>
#include <stdio.h>
#include <openacc.h>
#include "globals_acc.h"
#include "interp_utils_acc.h"
#include "create_xgrid_utils_acc.h"
#include "create_xgrid_acc.h"

typedef struct {
  int ncells_input;
  int *approx_nxcells_per_ij1;
  int *ij2_start;
  int *ij2_end;
  int upbound_nxcells;
} Answers;

void generate_input_is_same_as_output_grids(Grid_config *input_grid, Grid_config *output_grid,
                                            Answers *answers, Xinfo_per_input_tile *xgrid_answers);
int run_tests(Grid_config *input_grid,  Grid_config *output_grid, Grid_cells_struct_config *input_grid_cells,
              Grid_cells_struct_config *output_grid_cells, int **approx_nxcells_per_ij1, int **ij2_start, int **ij2_end,
              int *upbound_nxcells, Xinfo_per_input_tile *xgrid_answers);
void cleanup_test(Answers *answers, Grid_config *input_grid, Grid_config *output_grid,
                  Grid_cells_struct_config *input_grid_cells, Grid_cells_struct_config *output_grid_cells,
                  int *approx_nxcells_per_ij1, int *ij2_start, int *ij2_end);
void check_answers(const Answers *answers, int *approx_nxcells_per_ij1, int *ij2_start,
                   int *ij2_end, const int upbound_nxcells);
void check_ianswers(const int n, const int *answers, const int *checkme);

//TODO:  add more complicated tests
int main()
{

  Grid_config input_grid, output_grid;
  Grid_cells_struct_config input_grid_cells, output_grid_cells;
  int upbound_nxcells;
  int *approx_nxcells_per_ij1, *ij2_start, *ij2_end;

  Answers answers;
  Xinfo_per_input_tile xgrid_answers;

  printf("CHECKING CASE 1:  Input grid = Output_grid\n");
  generate_input_is_same_as_output_grids(&input_grid, &output_grid, &answers, &xgrid_answers);
  run_tests(&input_grid, &output_grid, &input_grid_cells, &output_grid_cells, &approx_nxcells_per_ij1, &ij2_start, &ij2_end,
            &upbound_nxcells, &xgrid_answers);
  check_answers(&answers, approx_nxcells_per_ij1, ij2_start, ij2_end, upbound_nxcells);
  cleanup_test(&answers, &input_grid, &output_grid, &input_grid_cells, &output_grid_cells, approx_nxcells_per_ij1,
               ij2_start, ij2_end);

  return 0;

}

void generate_input_is_same_as_output_grids(Grid_config *input_grid, Grid_config *output_grid,
                                            Answers *answers, Xinfo_per_input_tile *xgrid_answers)
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
  answers->ncells_input = ncells;
  answers->upbound_nxcells = ncells;
  answers->approx_nxcells_per_ij1 = (int *)malloc(ncells*sizeof(int));
  answers->ij2_start = (int *)malloc(ncells*sizeof(int));
  answers->ij2_end = (int *)malloc(ncells*sizeof(int));
  for( int ij1=0 ; ij1<ncells ; ij1++ ){
    answers->approx_nxcells_per_ij1[ij1] = 1;
    answers->ij2_start[ij1] = ij1;
    answers->ij2_end[ij1] = ij1;
  }

  //xgrid_answers
  xgrid_answers->input_parent_cell_indices = (int *)malloc(ncells*sizeof(int));
  xgrid_answers->output_parent_cell_indices = (int*)malloc(ncells*sizeof(int));
  xgrid_answers->xcell_area = (double *)malloc(ncells*sizeof(int));
  for(int ij1=0 ; ij1<ncells ; ij1++) {
    xgrid_answers->input_parent_cell_indices[ij1] = ij1;
    xgrid_answers->output_parent_cell_indices[ij1] = ij1;
  }

}


void check_answers(const Answers *answers, int *approx_nxcells_per_ij1, int *ij2_start,
                   int *ij2_end, const int upbound_nxcells)
{

  int ncells_input = answers->ncells_input;

  //double checking answers got copied out
  for(int icell=0 ; icell<ncells_input; icell++) {
    ij2_start[icell] = -99;
    ij2_end[icell] = -99;
    approx_nxcells_per_ij1[icell] = -99;
  }

#pragma acc exit data copyout( approx_nxcells_per_ij1[:ncells_input], \
                               ij2_start[:ncells_input],              \
                               ij2_end[:ncells_input])

  printf("checking upbound_nxcells\n");
  check_ianswers(1, &(answers->upbound_nxcells), &upbound_nxcells);

  printf("checking approx_nxcells_per_ij1\n");
  check_ianswers(ncells_input, answers->approx_nxcells_per_ij1, approx_nxcells_per_ij1);

  printf("checking ij2_start\n");
  check_ianswers(ncells_input, answers->ij2_start, ij2_start);

  printf("checking ij2_end\n");
  check_ianswers(ncells_input, answers->ij2_end, ij2_end);

}

int run_tests(Grid_config *input_grid, Grid_config *output_grid, Grid_cells_struct_config *input_grid_cells,
              Grid_cells_struct_config *output_grid_cells, int **approx_nxcells_per_ij1, int **ij2_start, int **ij2_end,
              int *upbound_nxcells, Xinfo_per_input_tile *xgrid_answers)
{

  int nlon_input_cells  = input_grid->nxc;
  int nlat_input_cells  = input_grid->nyc;
  int nlon_output_cells = output_grid->nxc;
  int nlat_output_cells = output_grid->nyc;
  int ncells_input  = nlon_input_cells * nlat_input_cells;
  int ncells_output = nlon_output_cells * nlat_output_cells;
  int ngridpts_input = (nlon_input_cells+1)*(nlat_input_cells+1);
  int ngridpts_output = (nlon_output_cells+1)*(nlat_output_cells+1);

  int nxcells;
  int jlat_overlap_starts, jlat_overlap_ends, nlat_overlapping_cells;

  int *p_approx_nxcells_per_ij1, *p_ij2_start, *p_ij2_end;

  Xinfo_per_input_tile xgrid;

  //copy grid to device
  copy_grid_to_device_acc(ngridpts_input, input_grid->latc, input_grid->lonc);
  copy_grid_to_device_acc(ngridpts_output, output_grid->latc, output_grid->lonc);

  //get mask to skip input cells in creating xgrid
  get_skip_cells_acc(ncells_input, &(input_grid_cells->skip_cells));

  //get output grid cell info
  get_grid_cells_struct_acc(nlon_output_cells, nlat_output_cells, output_grid->lonc, output_grid->latc,
                            output_grid_cells);

  //get bounding indices
  get_bounding_indices_acc(nlon_output_cells, nlat_output_cells, nlon_input_cells, nlat_input_cells,
                           output_grid->latc, input_grid->latc, &jlat_overlap_starts, &jlat_overlap_ends,
                           &nlat_overlapping_cells);
  if(jlat_overlap_starts != 0) printf("SMETHING IS WRONG WITH JLAT_OVERLAP_STARTS %d\n", jlat_overlap_starts);
  if(jlat_overlap_ends != nlat_input_cells-1) printf("SOMETHING IS WRONG WITH JLAT_OVERLAP_ENDS %d %d\n", jlat_overlap_ends, nlat_input_cells);
  if(nlat_overlapping_cells != nlat_input_cells) printf("SOMETHING IS WRONG WITH NLAT_OVERLAPPING_CELLS\n");

  //malloc and create arrays
  create_upbound_nxcells_arrays_on_device_acc(ncells_input, &p_approx_nxcells_per_ij1, &p_ij2_start, &p_ij2_end);

  *upbound_nxcells = get_upbound_nxcells_2dx2d_acc(input_grid->nxc, input_grid->nyc, output_grid->nxc, output_grid->nyc,
                                                   jlat_overlap_starts, jlat_overlap_ends,
                                                   input_grid->lonc, input_grid->latc, output_grid->lonc, output_grid->latc,
                                                   input_grid_cells->skip_cells, output_grid_cells,
                                                   p_approx_nxcells_per_ij1, p_ij2_start, p_ij2_end);

  create_xgrid_per_intile_arrays_on_device_acc(*upbound_nxcells, 0, &xgrid);
  nxcells = create_xgrid_2dx2d_order1_acc( nlon_input_cells, nlat_input_cells,
                                           nlon_output_cells, nlat_output_cells,
                                           jlat_overlap_starts, jlat_overlap_ends,
                                           input_grid->lonc, input_grid->latc,
                                           output_grid->lonc, output_grid->latc,
                                           *upbound_nxcells,
                                           input_grid_cells->skip_cells,
                                           output_grid_cells,
                                           p_approx_nxcells_per_ij1, p_ij2_start, p_ij2_end,
                                           &xgrid);

#pragma acc exit data copyout(xgrid.input_parent_cell_indices[:ncells_input], \
                              xgrid.output_parent_cell_indices[:ncells_input], \
                              xgrid.xcell_area[:ncells_input])

  for(int i=0 ; i<ncells_input ; i++) {
    if(xgrid_answers->input_parent_cell_indices[i] != xgrid.input_parent_cell_indices[i]) {
      printf("SOMETHING WRONG INPUT_PARENT_CELL_INDICES\n");
      exit(1);
    }
    if(xgrid_answers->output_parent_cell_indices[i] != xgrid.input_parent_cell_indices[i]) {
      printf("SOMETHING WRONG INPUT_PARENT_CELL_INDICES\n");
      exit(1);
    }


  }

  *approx_nxcells_per_ij1 = p_approx_nxcells_per_ij1;
  *ij2_start = p_ij2_start;
  *ij2_end = p_ij2_end;

  return 0;

}

void cleanup_test(Answers *answers, Grid_config *input_grid, Grid_config *output_grid,
                Grid_cells_struct_config *input_grid_cells, Grid_cells_struct_config *output_grid_cells,
                int *approx_nxcells_per_ij1, int *ij2_start, int *ij2_end)
{

  int ncells_output = output_grid->nxc * output_grid->nyc;
  int ncells_input  = input_grid->nxc * input_grid->nyc;

  free_output_grid_cell_struct_from_all_acc(ncells_output, output_grid_cells);
  free_upbound_xcells_array_from_all_acc(ncells_input, approx_nxcells_per_ij1, ij2_start, ij2_end);
  free_skip_cells_on_all_acc(ncells_input, input_grid_cells->skip_cells);

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

// void offset by half
