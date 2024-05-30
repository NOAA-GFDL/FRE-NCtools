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

void generate_input_is_same_as_output_grids(Grid_config *input_grid, Grid_config *output_grid, Answers *answers);
void generate_two_ouput_cells_per_input_cell(Grid_config *input_grid, Grid_config *output_grid, Answers *answers);
void check_data_on_device(Grid_config *input_grid, Grid_config *output_grid,
                          int *approx_nxcells_per_ij1, int *ij2_start, int *ij2_end,
                          Grid_cells_struct_config *input_grid_cells,
                          Grid_cells_struct_config *output_grid_cells);
int run_tests(const Grid_config *input_grid, const Grid_config *output_grid, Grid_cells_struct_config *input_grid_cells,
              Grid_cells_struct_config *output_grid_cells, int **approx_nxcells_per_ij1, int **ij2_start, int **ij2_end,
              int *upbound_nxcells);
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

  printf("CHECKING CASE 1:  Input grid = Output_grid\n");
  generate_input_is_same_as_output_grids(&input_grid, &output_grid, &answers);
  run_tests(&input_grid, &output_grid, &input_grid_cells, &output_grid_cells, &approx_nxcells_per_ij1, &ij2_start, &ij2_end,
            &upbound_nxcells);
  check_answers(&answers, approx_nxcells_per_ij1, ij2_start, ij2_end, upbound_nxcells);
  cleanup_test(&answers, &input_grid, &output_grid, &input_grid_cells, &output_grid_cells, approx_nxcells_per_ij1,
               ij2_start, ij2_end);

  generate_two_ouput_cells_per_input_cell(&input_grid, &output_grid, &answers);
  run_tests(&input_grid, &output_grid, &input_grid_cells, &output_grid_cells, &approx_nxcells_per_ij1, &ij2_start, &ij2_end,
            &upbound_nxcells);

  return 0;

}

void generate_input_is_same_as_output_grids(Grid_config *input_grid, Grid_config *output_grid, Answers *answers)
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

}

/*
void generate_two_ouput_cells_per_input_cell(Grid_config *input_grid, Grid_config *output_grid, Answers *answers)
{

  int nlat_input_cells = 2;   //[0, 30, 60]
  int nlon_input_cells = 180; //[0, 1, 2...180]

  int nlat_output_cells = 2*nlat_input_cells ; //[0, 15, 30, 45, 60]
  int nlon_output_cells = 2*nlon_input_cells ; //[0, 0.5, 1..]

  int ncells_input  = nlon_input_cells * nlat_input_cells;
  int ncells_output = nlon_output_cells * nlat_output_cells;
  int ngridpts_input  = (nlon_input_cells+1)*(nlat_input_cells+1);
  int ngridpts_output = (nlon_output_cells+1)*(nlat_output_cells+1);

  double dlat_input = 30.0;
  double dlon_input = 1.0;
  double dlat_output = 15.0;
  double dlon_output = 0.5;

  int i=0;

  input_grid->nxc  = nlon_input_cells  ; input_grid->nyc  = nlat_input_cells;
  output_grid->nxc = nlon_output_cells ; output_grid->nyc = nlat_output_cells;

  input_grid->latc = (double *)malloc(ngridpts_input*sizeof(double));
  input_grid->lonc = (double *)malloc(ngridpts_input*sizeof(double));
  output_grid->latc = (double *)malloc(ngridpts_output*sizeof(double));
  output_grid->lonc = (double *)malloc(ngridpts_output*sizeof(double));

  //input grid
  for( int ilat=0 ; ilat<nlat_input_cells+1 ; ilat++){
    double latitude = (dlat_input*ilat)*D2R ;
    for(int ilon=0; ilon<nlon_input_cells+1 ; ilon++) {
      input_grid->latc[ilat] = latitude;
      input_grid->lonc[ilat] = (ilon * dlon_input)*D2R ;
    }
  }

  for(int ilat=0 ; ilat<nlat_output_cells+1 ; ilat++) {
    double latitude = (dlat_output*ilat)*D2R;
    for(int ilon=0 ; ilon<nlon_output_cells+1; ilon++) {
      output_grid->latc[ilat] = latitude;
      output_grid->lonc[ilat] = (ilon *dlon_output) * D2R;
    }
  }

  //answers
  answers->ncells_input = ncells_input;
  answers->upbound_nxcells = ncells_output;
  answers->approx_nxcells_per_ij1 = (int *)malloc(ncells_input*sizeof(int));
  answers->ij2_start = (int *)malloc(ncells_input*sizeof(int));
  answers->ij2_end = (int *)malloc(ncells_input*sizeof(int));
  for( int ilat=0 ; ilat<nlat_input_cells; ilat++) {
    for( int ilon=0 ; ilon<nlon_input_cells ; ilon++) {
      answers->approx_nxcells_per_ij1[i] = 4;
      answers->ij2_start[i] = (ilat*2)*nlon_output_cells + ilon*2;
      answers->ij2_end[i] = (2*ilat+1)*nlon_output_cells + ilon*2+1;
      i++;
    }
  }

}
*/

void check_data_on_device( Grid_config *input_grid, Grid_config *output_grid,
                           int *approx_nxcells_per_ij1, int *ij2_start, int *ij2_end,
                           Grid_cells_struct_config *input_grid_cells,
                           Grid_cells_struct_config *output_grid_cells)
{

  int nlon_input_cells = input_grid->nxc;
  int nlat_input_cells = input_grid->nyc;
  int nlon_output_cells = output_grid->nxc;
  int nlat_output_cells = output_grid->nyc;

  int input_ngridpts  = (nlon_input_cells+1)*(nlat_input_cells+1);
  int output_ngridpts = (nlon_output_cells+1)*(nlat_output_cells+1);
  int input_ncells  = nlon_input_cells * nlat_input_cells;
  int output_ncells = nlon_output_cells * nlat_output_cells;

  // grid points should be on device.
  if(!acc_is_present(input_grid->lonc, input_ngridpts*sizeof(double))) {
    printf("INPUT GRID LON COORDINATES NOT ON DEVICE!"); exit(1);
  }
  if(!acc_is_present(input_grid->latc, input_ngridpts*sizeof(double))) {
    printf("INPUT GRID LAT COORDINATES NOT ON DEVICE!"); exit(1);
  }
  if(!acc_is_present(output_grid->lonc, output_ngridpts*sizeof(double))) {
    printf("OUTPUT GRID LON COORDINATES NOT ON DEVICE!"); exit(1);
  }
  if(!acc_is_present(output_grid->latc, output_ngridpts*sizeof(double))) {
    printf("OUTPUT GRID LAT COORDINATES NOT ON DEVICE!"); exit(1);
  }

  //skip input cells should be on device
  if(!acc_is_present(input_grid_cells->skip_cells, input_ncells*sizeof(double))) {
    printf("INPUT SKIP CELLS not on device!"); exit(1);
  }

  //arrays used by upbound_nxgrid should be on device
  if(!acc_is_present(approx_nxcells_per_ij1, input_ncells*sizeof(int))) {
    printf("APPROX_NXCELLS_PER_IJ1 NOT ON DEVICE!"); exit(1);
  }
  if(!acc_is_present(ij2_start, input_ncells*sizeof(int))) {
    printf("IJ2_START NOT ON DEVICE!"); exit(1);
  }
  if(!acc_is_present(ij2_end, input_ncells*sizeof(int))) {
    printf("IJ2_END NOT ON DEVICE!"); exit(1);
  }

  //output_grid_cells information should be on device
  if(!acc_is_present(output_grid_cells->lon_min, output_ncells*sizeof(double))) {
    printf("OUTPUT_GRID_CELLS LON_MIN NOT ON DEVICE!"); exit(1);
  }
  if(!acc_is_present(output_grid_cells->lon_max, output_ncells*sizeof(double))) {
    printf("OUTPUT_GRID_CELLS LON_MAX NOT ON DEVICE!"); exit(1);
  }
  if(!acc_is_present(output_grid_cells->lat_min, output_ncells*sizeof(double))) {
    printf("OUTPUT_GRID_CELLS LAT_MIN NOT ON DEVICE!"); exit(1);
  }
  if(!acc_is_present(output_grid_cells->lat_max, output_ncells*sizeof(double))) {
    printf("OUTPUT_GRID_CELLS LAT_MAX NOT ON DEVICE!"); exit(1);
  }
  if(!acc_is_present(output_grid_cells->lon_cent, output_ncells*sizeof(double))) {
    printf("OUTPUT_GRID_CELLS LON_CENT NOT ON DEVICE!"); exit(1);
  }
  if(!acc_is_present(output_grid_cells->area, output_ncells*sizeof(double))) {
    printf("OUTPUT_GRID_CELLS AREA NOT ON DEVICE!"); exit(1);
  }
  if(!acc_is_present(output_grid_cells->nvertices, output_ncells*sizeof(int))) {
    printf("OUTPUT_GRID_CELLS NVERTICES NOT ON DEVICE!"); exit(1);
  }

  for(int icell=0 ; icell<output_ncells; icell++) {
    if(!acc_is_present(output_grid_cells->lon_vertices[icell], MAX_V*sizeof(double))) {
    printf("OUTPUT_GRID_CELLS LON VERTICES NOT ON DEVICE!"); exit(1);
    }
    if(!acc_is_present(output_grid_cells->lat_vertices[icell], MAX_V*sizeof(double))) {
      printf("OUTPUT_GRID_CELLS LAT VERTICES NOT ON DEVICE!"); exit(1);
    }
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

int run_tests(const Grid_config *input_grid, const Grid_config *output_grid, Grid_cells_struct_config *input_grid_cells,
              Grid_cells_struct_config *output_grid_cells, int **approx_nxcells_per_ij1, int **ij2_start, int **ij2_end,
              int *upbound_nxcells)
{

  int nlon_input_cells  = input_grid->nxc;
  int nlat_input_cells  = input_grid->nyc;
  int nlon_output_cells = output_grid->nxc;
  int nlat_output_cells = output_grid->nyc;
  int ncells_input  = nlon_input_cells * nlat_input_cells;
  int ncells_output = nlon_output_cells * nlat_output_cells;
  int ngridpts_input = (nlon_input_cells+1)*(nlat_input_cells+1);
  int ngridpts_output = (nlon_output_cells+1)*(nlat_output_cells+1);

  int *p_approx_nxcells_per_ij1, *p_ij2_start, *p_ij2_end;

  //copy grid to device
  copy_grid_to_device_acc(ngridpts_input, input_grid->latc, input_grid->lonc);
  copy_grid_to_device_acc(ngridpts_output, output_grid->latc, output_grid->lonc);

  //get mask to skip input cells in creating xgrid
  get_skip_cells_acc(ncells_input, &(input_grid_cells->skip_cells));

  //get output grid cell info
  get_grid_cells_struct_acc(nlon_output_cells, nlat_output_cells, output_grid->lonc, output_grid->latc,
                            output_grid_cells);

  //malloc and create arrays
  create_upbound_nxcells_arrays_on_device_acc(ncells_input, &p_approx_nxcells_per_ij1, &p_ij2_start, &p_ij2_end);

  check_data_on_device(input_grid, output_grid, p_approx_nxcells_per_ij1, p_ij2_start, p_ij2_end,
                       input_grid_cells, output_grid_cells);

  *upbound_nxcells = get_upbound_nxcells_2dx2d_acc(input_grid->nxc, input_grid->nyc, output_grid->nxc, output_grid->nyc,
                                                   input_grid->lonc, input_grid->latc, output_grid->lonc, output_grid->latc,
                                                   input_grid_cells->skip_cells, output_grid_cells,
                                                   p_approx_nxcells_per_ij1, p_ij2_start, p_ij2_end);

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
