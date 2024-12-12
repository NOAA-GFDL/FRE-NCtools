/***********************************************************************
 *                   GNU Lesser General Public License
 *
 * This file is part of the GFDL FRE NetCDF tools package (FRE-NCTools).
 *
 * FRE-NCtools is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at
 * your option) any loner version.
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

// This test tests function read_remap_file_gpu to read in a made-up
// remap file.  It ensures the correct initialization of the interp_gpu struct.
// The remap file is generated with the python script
// test_make_remap_file_conserve.py that uses the xarray module.  This
// test also tests the function copy_interp_to_device_gpu which copies interp_gpu
// to device.

#include <openacc.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <unistd.h>
#include "conserve_interp_gpu.h"
#include "interp_utils_gpu.h"
#include "globals_gpu.h"

#define INPUT_GRID_NTILES 6
#define OUTPUT_GRID_NTILES 2

double const tolerance = 1.e-7;

int const input_grid_nlon = 2;
int const output_grid_nlon = 3;

typedef struct {
  char answer_file[15];
  size_t total_nxcells;
  int nxcells[INPUT_GRID_NTILES];
  int *input_parent_cell_index[INPUT_GRID_NTILES];
  int *output_parent_cell_index[INPUT_GRID_NTILES];
  double *xcell_area[INPUT_GRID_NTILES];
  double *dcentroid_lon[INPUT_GRID_NTILES];
  double *dcentroid_lat[INPUT_GRID_NTILES];
} Answers;

typedef enum {
  myCONSERVE_ORDER1,
  myCONSERVE_ORDER2
} myInterp_Method;

typedef enum {
  ON_DEVICE,
  ON_HOST
} Where_Am_I;

char answer_files1[OUTPUT_GRID_NTILES][30] = { "answers_conserve1.tile1.txt", "answers_conserve1.tile2.txt" };
char answer_files2[OUTPUT_GRID_NTILES][30] = { "answers_conserve2.tile1.txt", "answers_conserve2.tile2.txt" };

char remap_files1[OUTPUT_GRID_NTILES][30] = { "remap_conserve1.tile1.nc", "remap_conserve1.tile2.nc" };
char remap_files2[OUTPUT_GRID_NTILES][30] = { "remap_conserve2.tile1.nc", "remap_conserve2.tile2.nc" };

void read_all_answers(Answers *answers, int myinterp_method);
void read_ianswers( FILE *myfile, int *nxcells, int **ianswer);
void read_ranswers( FILE *myfile, int *nxcells, double **ianswer);
void check_answers_on_device(Answers *answers, Interp_config_gpu *interp_gpu, int myinterp_method);
void check_answers_on_host(Answers *answers, Interp_config_gpu *interp_gpu, int myinterp_method);
void reset_interp_gpu_on_host( Interp_config_gpu *interp_gpu, Answers *answers, int myinterp_method );
void check_ianswers(int n, int *answers, int *checkme, int host_or_device);
void check_ranswers(int n, double *answers, double *checkme, int host_or_device);
void error(char *error_message);

// start program
int main(int argc, char *argv[]) {

  Interp_config_gpu interp_gpu[OUTPUT_GRID_NTILES];
  Grid_config       input_grid[INPUT_GRID_NTILES], output_grid[OUTPUT_GRID_NTILES];
  Answers           answers[OUTPUT_GRID_NTILES];

  unsigned int opcode;
  myInterp_Method myinterp_method;

  if( atoi(argv[1]) == 1) {
    opcode = CONSERVE_ORDER1;
    myinterp_method = myCONSERVE_ORDER1;
  }
  else if( atoi(argv[1]) == 2 ) {
    opcode = CONSERVE_ORDER2;
    myinterp_method = myCONSERVE_ORDER2;
  }

  // assign number of cells in lon direction to Grid_config
  for(int otile=0 ; otile<OUTPUT_GRID_NTILES ; otile++) output_grid[otile].nxc = output_grid_nlon;
  for(int itile=0 ; itile<INPUT_GRID_NTILES ; itile++) input_grid[itile].nxc = input_grid_nlon;


  // assign remap files
  for(int n=0 ; n<OUTPUT_GRID_NTILES ; n++) {
    interp_gpu[n].input_tile = ( Interp_per_input_tile *)malloc(INPUT_GRID_NTILES*sizeof(Interp_per_input_tile));
    if(myinterp_method == myCONSERVE_ORDER1) strcpy(interp_gpu[n].remap_file, remap_files1[n]);
    if(myinterp_method == myCONSERVE_ORDER2) strcpy(interp_gpu[n].remap_file, remap_files2[n]);
    if(access(interp_gpu[n].remap_file, F_OK)==0) interp_gpu[n].file_exist=1;
  }

  // assign answer files
  for(int n=0 ; n<OUTPUT_GRID_NTILES ; n++) {
    if(myinterp_method == myCONSERVE_ORDER1) strcpy(answers[n].answer_file, answer_files1[n]);
    if(myinterp_method == myCONSERVE_ORDER2) strcpy(answers[n].answer_file, answer_files2[n]);
  }

  // read in remap and transfer data to device
  read_remap_file_gpu(INPUT_GRID_NTILES, OUTPUT_GRID_NTILES, &output_grid, &input_grid, interp_gpu, opcode);
  copy_interp_to_device_gpu(INPUT_GRID_NTILES, OUTPUT_GRID_NTILES, interp_gpu, opcode) ;

  // read in answers from txt file
  read_all_answers(answers, myinterp_method);

  // check answers on device
  check_answers_on_device(answers, interp_gpu, myinterp_method);

  // compare answers on host
  check_answers_on_host( answers, interp_gpu, myinterp_method);

  //test success
  return 0;

}
//------------------------------------------
//------------------------------------------
void read_ianswers( FILE *myfile, int *nxcells, int **ianswer)
{

  for( int m=0 ; m<INPUT_GRID_NTILES ; m++ ) {
    int inxcells = nxcells[m] ;
    for( int j=0 ; j<inxcells ; j++ ) fscanf( myfile, "%d", ianswer[m]+j );
  }

}
//------------------------------------------
//------------------------------------------
void read_ranswers( FILE *myfile, int *nxcells, double **ranswer)
{

  for( int m=0 ; m<INPUT_GRID_NTILES ; m++ ) {
    int inxcells = nxcells[m] ;
    for( int j=0 ; j<inxcells ; j++ ) fscanf( myfile, "%lf", ranswer[m]+j );
  }

}
//------------------------------------------
//------------------------------------------
void read_all_answers(Answers *answers, int myinterp_method)
{

  for(int n=0 ; n<OUTPUT_GRID_NTILES ; n++) {

    FILE * myfile = fopen(answers[n].answer_file, "r");
    int **lon_index, **lat_index;

    printf("READING IN ANSWERS FOR n=%d\n", n);

    // get answers for interp_gpu[n].nxcells of type size_t
    fscanf( myfile, "%zu", &(answers[n].total_nxcells) );

    // get answers for interp_gpu[n].input_tile[m].nxcells
    for( int m=0 ; m<INPUT_GRID_NTILES ; m++ ) fscanf( myfile, "%d", answers[n].nxcells+m );

    // temporary arrays to read in lon and lat indices
    lon_index = (int **)malloc(INPUT_GRID_NTILES*sizeof(int *));
    lat_index = (int **)malloc(INPUT_GRID_NTILES*sizeof(int *));
    for(int m=0 ; m<INPUT_GRID_NTILES ; m++) {
      lon_index[m] = (int *)malloc(answers[n].nxcells[m]*sizeof(int));
      lat_index[m] = (int *)malloc(answers[n].nxcells[m]*sizeof(int));
    }

    // allocate arrays for answerss
    for( int m=0 ; m<INPUT_GRID_NTILES ; m++ ) {
      int inxcells = answers[n].nxcells[m] ;
      answers[n].input_parent_cell_index[m]  = (int *)malloc( inxcells * sizeof(int) ) ;
      answers[n].output_parent_cell_index[m] = (int *)malloc( inxcells * sizeof(int) ) ;
      answers[n].xcell_area[m]  = (double *)malloc( inxcells * sizeof(double) );
      if( myinterp_method == myCONSERVE_ORDER2 ) {
        answers[n].dcentroid_lon[m] = (double *)malloc( inxcells * sizeof(double) );
        answers[n].dcentroid_lat[m] = (double *)malloc( inxcells * sizeof(double) );
      }
    }

    //read in input parent grid indices
    read_ianswers( myfile, answers[n].nxcells, lon_index) ;
    read_ianswers( myfile, answers[n].nxcells, lat_index);
    //compute ij index
    for( int m=0 ; m<INPUT_GRID_NTILES ; m++) {
      for( int i=0 ; i<answers[n].nxcells[m] ; i++ ) {
        answers[n].input_parent_cell_index[m][i] = lat_index[m][i]*input_grid_nlon + lon_index[m][i];
      }
    }

    //read in output parent grid indices
    read_ianswers( myfile, answers[n].nxcells, lon_index) ;
    read_ianswers( myfile, answers[n].nxcells, lat_index);
    //sort ij index
    for( int m=0 ; m<INPUT_GRID_NTILES ; m++) {
      for( int i=0 ; i<answers[n].nxcells[m] ; i++ ) {
        answers[n].output_parent_cell_index[m][i] = lat_index[m][i]*output_grid_nlon + lon_index[m][i];
      }
    }

    // read in xcell_area answers
    read_ranswers( myfile, answers[n].nxcells, answers[n].xcell_area);

    //read dcentroid_lon dcentroid_lat if conserve_order2
    if( myinterp_method == myCONSERVE_ORDER2 ) {
      read_ranswers( myfile, answers[n].nxcells, answers[n].dcentroid_lon);
      read_ranswers( myfile, answers[n].nxcells, answers[n].dcentroid_lat);
    }

    for(int m=0 ; m<INPUT_GRID_NTILES; m++) {
      free(lon_index[m]);
      free(lat_index[m]);
    }
    free(lon_index); lon_index = NULL;
    free(lat_index); lat_index = NULL;

    fclose( myfile );

  }

}
//------------------------------------------
//------------------------------------------
void check_answers_on_device( Answers *answers, Interp_config_gpu *interp_gpu, int myinterp_method )
{

  int *p_input_parent_cell_index, *p_output_parent_cell_index;
  double *p_xcell_area, *p_dcentroid_lon, *p_dcentroid_lat;

  printf("CHECKING ANSWERS ON DEVICE\n");

  if( !acc_is_present(interp_gpu, OUTPUT_GRID_NTILES*sizeof(Interp_config_gpu)) )
    error("ERROR interp_gpu is not present") ;

  for(int n=0 ; n<OUTPUT_GRID_NTILES ; n++) {

    size_t nxcells = answers[n].total_nxcells;

    printf("checking for n=%d\n", n);

#pragma acc parallel loop seq copyin(nxcells) present(interp_gpu[n].nxcells)
    for(int i=0 ; i<1 ; i++ ){
      if( nxcells != interp_gpu[n].nxcells ) {
        printf("ERROR checking total number of nxcells for n=0", n);
        nxcells = 0;
      }
    }
    if( nxcells == 0 ) error("ERROR with the value of nxcells on device");

    for( int m=0 ; m<INPUT_GRID_NTILES ; m++ ) {

      int inxcells = answers[n].nxcells[m];

      //copy in answers
      p_input_parent_cell_index  = answers[n].input_parent_cell_index[m];
      p_output_parent_cell_index = answers[n].output_parent_cell_index[m];
      acc_copyin(&inxcells, sizeof(int));
      acc_copyin(p_input_parent_cell_index, inxcells*sizeof(int));
      acc_copyin(p_output_parent_cell_index, inxcells*sizeof(int));
      p_xcell_area  = answers[n].xcell_area[m];  acc_copyin(p_xcell_area, inxcells*sizeof(double));
      if(myinterp_method == myCONSERVE_ORDER2) {
        p_dcentroid_lon = answers[n].dcentroid_lon[m]; acc_copyin(p_dcentroid_lon, inxcells*sizeof(double));
        p_dcentroid_lat = answers[n].dcentroid_lat[m]; acc_copyin(p_dcentroid_lat, inxcells*sizeof(double));
      }

      // check answers
      printf("m=%d ", m);
      printf("nxcells "); check_ianswers(1, &inxcells, &(interp_gpu[n].input_tile[m].nxcells),  ON_DEVICE);

      printf("input_parent_cell_index ");
      check_ianswers(inxcells, p_input_parent_cell_index, interp_gpu[n].input_tile[m].input_parent_cell_index,  ON_DEVICE);

      printf("output_parent_cell_index ");
      check_ianswers(inxcells, p_output_parent_cell_index, interp_gpu[n].input_tile[m].output_parent_cell_index, ON_DEVICE);

      printf("xcell_area ");
      check_ranswers(inxcells, p_xcell_area,  interp_gpu[n].input_tile[m].xcell_area,  ON_DEVICE);

      if( myinterp_method == myCONSERVE_ORDER2 ){
        printf("dcentroid_lon ");
        check_ranswers(inxcells, p_dcentroid_lon, interp_gpu[n].input_tile[m].dcentroid_lon, ON_DEVICE);
        printf("dcentroid_lat ");
        check_ranswers(inxcells, p_dcentroid_lat, interp_gpu[n].input_tile[m].dcentroid_lat, ON_DEVICE);
      }
      printf("\n");

      acc_delete(&inxcells, sizeof(int));
      acc_delete(p_input_parent_cell_index, inxcells*sizeof(int));
      acc_delete(p_output_parent_cell_index, inxcells*sizeof(int));
      acc_delete(p_xcell_area, inxcells*sizeof(double));
      if(myinterp_method == myCONSERVE_ORDER2) {
        acc_delete(p_dcentroid_lon, inxcells*sizeof(double));
        acc_delete(p_dcentroid_lat, inxcells*sizeof(double));
      }
    }
  }

}
//------------------------------------------
//------------------------------------------
void check_answers_on_host( Answers *answers, Interp_config_gpu *interp_gpu, int myinterp_method )
{

  printf("CHECKING COPIED OUT DATA ON HOST\n");

  //cautious step.  to ensure data is actually copied out
  reset_interp_gpu_on_host(interp_gpu, answers, myinterp_method);

  // copyout answers
  for(int n=0 ; n<OUTPUT_GRID_NTILES; n++) {
    acc_update_host( &interp_gpu[n].nxcells, sizeof(size_t) );
    for(int m=0 ; m<INPUT_GRID_NTILES ; m++) {
      int inxcells = answers[n].nxcells[m];
      acc_update_host( &interp_gpu[n].input_tile[m].nxcells, sizeof(int));
      acc_copyout( interp_gpu[n].input_tile[m].input_parent_cell_index,  inxcells*sizeof(int));
      acc_copyout( interp_gpu[n].input_tile[m].output_parent_cell_index, inxcells*sizeof(int));
      acc_copyout( interp_gpu[n].input_tile[m].xcell_area,  inxcells*sizeof(double));
      if(myinterp_method == myCONSERVE_ORDER2) {
        acc_copyout( interp_gpu[n].input_tile[m].dcentroid_lon, inxcells*sizeof(double));
        acc_copyout( interp_gpu[n].input_tile[m].dcentroid_lat, inxcells*sizeof(double));
      }
    }
  }
  acc_delete(interp_gpu, OUTPUT_GRID_NTILES*sizeof(Interp_config_gpu));

  for( int n=0 ; n<OUTPUT_GRID_NTILES; n++ ) {

    printf("checking nxcells for n=%d\n", n);
    if( answers[n].total_nxcells != interp_gpu[n].nxcells ) {
      printf("ERROR nxcells for n=%d: %zu vs %zu\n", answers[n].total_nxcells, interp_gpu[n].nxcells);
      error("goodbye!");
    }

    for( int m=0 ; m<INPUT_GRID_NTILES ; m++ ) {
      printf("m=%d ", m);
      int inxcells = interp_gpu[n].input_tile[m].nxcells;

      printf("nxcells "); check_ianswers(1, answers[n].nxcells+m, &(interp_gpu[n].input_tile[m].nxcells), ON_HOST);

      printf("input_parent_lon_index ");
      check_ianswers(inxcells, answers[n].input_parent_cell_index[m],  interp_gpu[n].input_tile[m].input_parent_cell_index,  ON_HOST);

      printf("output_parent_lon_index ");
      check_ianswers(inxcells, answers[n].output_parent_cell_index[m], interp_gpu[n].input_tile[m].output_parent_cell_index, ON_HOST);

      printf("xcell_area ");
      check_ranswers(inxcells, answers[n].xcell_area[m],  interp_gpu[n].input_tile[m].xcell_area,  ON_HOST);

      if(myinterp_method == myCONSERVE_ORDER2) {
        printf("dcentroid_lon ");
        check_ranswers(inxcells, answers[n].dcentroid_lon[m],  interp_gpu[n].input_tile[m].dcentroid_lon, ON_HOST);
        printf("dcentroid_lat");
        check_ranswers(inxcells, answers[n].dcentroid_lat[m],  interp_gpu[n].input_tile[m].dcentroid_lat, ON_HOST);
      }
      printf("\n");
    }
  }

}
//------------------------------------------
//------------------------------------------
void reset_interp_gpu_on_host( Interp_config_gpu *interp_gpu, Answers *answers, int myinterp_method )
{

  for(int n=0 ; n<OUTPUT_GRID_NTILES ; n++) {
    interp_gpu[n].nxcells=0;
    for(int m=0 ; m<INPUT_GRID_NTILES ; m++) {
      int inxcells = answers[n].nxcells[m];
      interp_gpu[n].input_tile[m].nxcells = -99;
      for( int i=0 ; i<inxcells ; i++) {
        interp_gpu[n].input_tile[m].input_parent_cell_index[i]  = -99 ;
        interp_gpu[n].input_tile[m].output_parent_cell_index[i] = -99 ;
        interp_gpu[n].input_tile[m].xcell_area[i]  = -99.9 ;
        if(myinterp_method == myCONSERVE_ORDER2) {
          interp_gpu[n].input_tile[m].dcentroid_lon[i] = -99.9;
          interp_gpu[n].input_tile[m].dcentroid_lat[i] = -99.9;
        }
      }
    }
  }

}
//------------------------------------------
//------------------------------------------
void check_ianswers( int n, int *answers, int *checkme, int host_or_device )
{

  if( host_or_device == ON_HOST ) {
    for(int i=0 ; i<n ; i++ ){
      if( answers[i] != checkme[i] ) {
        printf("ERROR element %d:  %d vs %d\n", i, answers[i], checkme[i]);
        exit(1);
      }
    }
  }
  else {

    int theres_an_error=0;

    if( !acc_is_present(answers, n*sizeof(int)) ) error("answers are not on device");
    if( !acc_is_present(checkme, n*sizeof(int)) ) error("answers are not on device");

#pragma acc data present(answers[:n], checkme[:n]) copy(theres_an_error)
#pragma acc parallel loop independent reduction(+:theres_an_error)
    for(int i=0 ; i<n ; i++ ){
      if( answers[i] != checkme[i] ) {
        printf("ERROR element %d:  %d vs %d\n", i, answers[i], checkme[i]);
        theres_an_error++;
      }
    }

    if(theres_an_error > 0 ) error("ERROR with data on device");

  }

}
//------------------------------------------
//------------------------------------------
void check_ranswers( int n, double *answers, double *checkme, int host_or_device)
{

  if( host_or_device == ON_HOST ) {
    for(int i=0 ; i<n ; i++ ){
      if( abs(answers[i]-checkme[i]) > tolerance ) {
        printf("ERROR element %d:  %lf vs %lf\n", i, answers[i], checkme[i]);
        exit(1);
      }
    }
  }
  else {

    int theres_an_error=-99;

    if( !acc_is_present(answers, n*sizeof(int)) ) error("answers are not on device");
    if( !acc_is_present(checkme, n*sizeof(int)) ) error("answers are not on device");

#pragma acc data present(answers[:n], checkme[:n]) copy(theres_an_error)
#pragma acc parallel loop independent reduction(+:theres_an_error)
    for(int i=0 ; i<n ; i++ ){
      if( abs(answers[i]-checkme[i]) > tolerance ) {
        printf("ERROR element %d:  %lf vs %lf\n", i, answers[i], checkme[i]);
        theres_an_error++;
      }
    }
    if(theres_an_error != -99) error("ERROR with data on device");
  }

}
//------------------------------------------
//------------------------------------------
void error( char *error_message )
{
  printf("%s\n", error_message);
  exit(1);
}
