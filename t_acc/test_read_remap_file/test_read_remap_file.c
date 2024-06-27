#include <openacc.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <unistd.h>
#include "conserve_interp_acc.h"
#include "interp_utils_acc.h"
#include "globals.h"

#define NTILES_IN 6
#define NTILES_OUT 2

double const tolerance = 1.e-7;

typedef struct {
  char answer_file[15];
  size_t total_nxcells;
  int nxcells[NTILES_IN];
  int *input_parent_lon_indices[NTILES_IN];
  int *input_parent_lat_indices[NTILES_IN];
  int *output_parent_lon_indices[NTILES_IN];
  int *output_parent_lat_indices[NTILES_IN];
  double *xcell_area[NTILES_IN];
  double *dcentroid_lon[NTILES_IN];
  double *dcentroid_lat[NTILES_IN];
} Answers;

typedef enum {
  myCONSERVE_ORDER1,
  myCONSERVE_ORDER2
} myInterp_Method;

typedef enum {
  ON_DEVICE,
  ON_HOST
} Where_Am_I;

char answer_files1[NTILES_OUT][30] = { "answers_conserve1.tile1.txt", "answers_conserve1.tile2.txt" };
char answer_files2[NTILES_OUT][30] = { "answers_conserve2.tile1.txt", "answers_conserve2.tile2.txt" };

char remap_files1[NTILES_OUT][30] = { "remap_conserve1.tile1.nc", "remap_conserve1.tile2.nc" };
char remap_files2[NTILES_OUT][30] = { "remap_conserve2.tile1.nc", "remap_conserve2.tile2.nc" };

void read_all_answers(Answers *answers, int myinterp_method);
void read_ianswers( FILE *myfile, int *nxcells, int **ianswer);
void read_ranswers( FILE *myfile, int *nxcells, double **ianswer);
void check_answers_on_device(Answers *answers, Xgrid_config *xgrid, int myinterp_method);
void check_answers_on_host(Answers *answers, Xgrid_config *xgrid, int myinterp_method);
void reset_xgrid_on_host( Xgrid_config *xgrid, Answers *answers, int myinterp_method );
void check_ianswers(int n, int *answers, int *checkme, int host_or_device);
void check_ranswers(int n, double *answers, double *checkme, int host_or_device);
void error(char *error_message);

// start program
int main(int argc, char *argv[]) {

  Xgrid_config xgrid[NTILES_OUT];
  Answers answers[NTILES_OUT];

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

  // assign remap files
  for(int n=0 ; n<NTILES_OUT ; n++) {
    xgrid[n].per_intile = ( Xinfo_per_input_tile *)malloc(NTILES_IN*sizeof(Xinfo_per_input_tile));
    if(myinterp_method == myCONSERVE_ORDER1) strcpy(xgrid[n].remap_file, remap_files1[n]);
    if(myinterp_method == myCONSERVE_ORDER2) strcpy(xgrid[n].remap_file, remap_files2[n]);
    if(access(xgrid[n].remap_file, F_OK)==0) xgrid[n].file_exist=1;
  }

  // read in remap and transfer data to device
  read_remap_file_acc(NTILES_IN, NTILES_OUT, xgrid, opcode) ;
  copy_xgrid_to_device_acc(NTILES_IN, NTILES_OUT, xgrid, opcode) ;

  // get all answers
  for(int n=0 ; n<NTILES_OUT ; n++) {
    if(myinterp_method == myCONSERVE_ORDER1) strcpy(answers[n].answer_file, answer_files1[n]);
    if(myinterp_method == myCONSERVE_ORDER2) strcpy(answers[n].answer_file, answer_files2[n]);
  }
  read_all_answers(answers, myinterp_method);

  // check answers on device
  check_answers_on_device(answers, xgrid, myinterp_method);

  // compare answers on host
  check_answers_on_host( answers, xgrid, myinterp_method);

  return 0; //test success

}
//------------------------------------------
//------------------------------------------
void read_ianswers( FILE *myfile, int *nxcells, int **ianswer)
{

  for( int m=0 ; m<NTILES_IN ; m++ ) {
    int inxcells = nxcells[m] ;
    for( int j=0 ; j<inxcells ; j++ ) fscanf( myfile, "%d", ianswer[m]+j );
  }

}
//------------------------------------------
//------------------------------------------
void read_ranswers( FILE *myfile, int *nxcells, double **ranswer)
{

  for( int m=0 ; m<NTILES_IN ; m++ ) {
    int inxcells = nxcells[m] ;
    for( int j=0 ; j<inxcells ; j++ ) fscanf( myfile, "%lf", ranswer[m]+j );
  }

}
//------------------------------------------
//------------------------------------------
void read_all_answers(Answers *answers, int myinterp_method)
{

  for(int n=0 ; n<NTILES_OUT ; n++) {

    FILE * myfile = fopen(answers[n].answer_file, "r");

    printf("READING IN ANSWERS FOR n=%d\n", n);

    // get answers for xgrid[n].nxcells of type size_t
    fscanf( myfile, "%zu", &(answers[n].total_nxcells) );

    // get answers for xgrid[n].per_intile[m].nxcells
    for( int m=0 ; m<NTILES_IN ; m++ ) fscanf( myfile, "%d", answers[n].nxcells+m );

    // allocate arrays for answerss
    for( int m=0 ; m<NTILES_IN ; m++ ) {
      int inxcells = answers[n].nxcells[m] ;
      answers[n].input_parent_lon_indices[m]  = (int *)malloc( inxcells * sizeof(int) ) ;
      answers[n].input_parent_lat_indices[m]  = (int *)malloc( inxcells * sizeof(int) ) ;
      answers[n].output_parent_lon_indices[m] = (int *)malloc( inxcells * sizeof(int) ) ;
      answers[n].output_parent_lat_indices[m] = (int *)malloc( inxcells * sizeof(int) ) ;
      answers[n].xcell_area[m]  = (double *)malloc( inxcells * sizeof(double) );
      if( myinterp_method == myCONSERVE_ORDER2 ) {
        answers[n].dcentroid_lon[m] = (double *)malloc( inxcells * sizeof(double) );
        answers[n].dcentroid_lat[m] = (double *)malloc( inxcells * sizeof(double) );
      }
    }

    //read in answers
    read_ianswers( myfile, answers[n].nxcells, answers[n].input_parent_lon_indices) ;
    read_ianswers( myfile, answers[n].nxcells, answers[n].input_parent_lat_indices);
    read_ianswers( myfile, answers[n].nxcells, answers[n].output_parent_lon_indices);
    read_ianswers( myfile, answers[n].nxcells, answers[n].output_parent_lat_indices);
    read_ranswers( myfile, answers[n].nxcells, answers[n].xcell_area);

    //read dcentroid_lon dcentroid_lat if conserve_order2
    if( myinterp_method == myCONSERVE_ORDER2 ) {
      read_ranswers( myfile, answers[n].nxcells, answers[n].dcentroid_lon);
      read_ranswers( myfile, answers[n].nxcells, answers[n].dcentroid_lat);
    }

    fclose( myfile );

  }

}
//------------------------------------------
//------------------------------------------
void check_answers_on_device( Answers *answers, Xgrid_config *xgrid, int myinterp_method )
{

  int *p_input_parent_lon_indices, *p_input_parent_lat_indices;
  int *p_output_parent_lon_indices, *p_output_parent_lat_indices;
  double *p_xcell_area, *p_dcentroid_lon, *p_dcentroid_lat;

  printf("CHECKING ANSWERS ON DEVICE\n");

  if( !acc_is_present(xgrid, NTILES_OUT*sizeof(Xgrid_config)) ) error("ERROR interp is not present") ;

  for(int n=0 ; n<NTILES_OUT ; n++) {

    size_t nxcells = answers[n].total_nxcells;

    printf("checking for n=%d\n", n);

#pragma acc parallel loop seq copyin(nxcells) present(xgrid[n].nxcells)
    for(int i=0 ; i<1 ; i++ ){
      if( nxcells != xgrid[n].nxcells ) {
        printf("ERROR checking total number of nxcells for n=0", n);
        nxcells = 0;
      }
    }
    if( nxcells == 0 ) error("ERROR with the value of nxcells on device");

    for( int m=0 ; m<NTILES_IN ; m++ ) {

      int inxcells = answers[n].nxcells[m];

      //copy in answers
      p_input_parent_lon_indices  = answers[n].input_parent_lon_indices[m];
      p_input_parent_lat_indices  = answers[n].input_parent_lat_indices[m];
      p_output_parent_lon_indices = answers[n].output_parent_lon_indices[m];
      p_output_parent_lat_indices = answers[n].output_parent_lat_indices[m];
      acc_copyin(&inxcells, sizeof(int));
      acc_copyin(p_input_parent_lon_indices, inxcells*sizeof(int));
      acc_copyin(p_input_parent_lat_indices, inxcells*sizeof(int));
      acc_copyin(p_output_parent_lon_indices, inxcells*sizeof(int));
      acc_copyin(p_output_parent_lat_indices, inxcells*sizeof(int));
      p_xcell_area  = answers[n].xcell_area[m];  acc_copyin(p_xcell_area, inxcells*sizeof(double));
      if(myinterp_method == myCONSERVE_ORDER2) {
        p_dcentroid_lon = answers[n].dcentroid_lon[m]; acc_copyin(p_dcentroid_lon, inxcells*sizeof(double));
        p_dcentroid_lat = answers[n].dcentroid_lat[m]; acc_copyin(p_dcentroid_lat, inxcells*sizeof(double));
      }

      printf("m=%d ", m);
      printf("nxcells "); check_ianswers(1, &inxcells, &(xgrid[n].per_intile[m].nxcells),  ON_DEVICE);
      printf("input_parent_lon_indices ");
      check_ianswers(inxcells, p_input_parent_lon_indices, xgrid[n].per_intile[m].input_parent_lon_indices,  ON_DEVICE);
      printf("input_parent_lat_indices ");
      check_ianswers(inxcells, p_input_parent_lat_indices, xgrid[n].per_intile[m].input_parent_lat_indices,  ON_DEVICE);
      printf("output_parent_lon_indices ");
      check_ianswers(inxcells, p_output_parent_lon_indices, xgrid[n].per_intile[m].output_parent_lon_indices, ON_DEVICE);
      printf("output_parent_lat_indices ");
      check_ianswers(inxcells, p_output_parent_lat_indices, xgrid[n].per_intile[m].output_parent_lat_indices, ON_DEVICE);
      printf("xcell_area ");
      check_ranswers(inxcells, p_xcell_area,  xgrid[n].per_intile[m].xcell_area,  ON_DEVICE);
      if( myinterp_method == myCONSERVE_ORDER2 ){
        printf("dcentroid_lon ");
        check_ranswers(inxcells, p_dcentroid_lon, xgrid[n].per_intile[m].dcentroid_lon, ON_DEVICE);
        printf("dcentroid_lat ");
        check_ranswers(inxcells, p_dcentroid_lat, xgrid[n].per_intile[m].dcentroid_lat, ON_DEVICE);
      }
      printf("\n");

      acc_delete(&inxcells, sizeof(int));
      acc_delete(p_input_parent_lon_indices, inxcells*sizeof(int));
      acc_delete(p_input_parent_lat_indices, inxcells*sizeof(int));
      acc_delete(p_output_parent_lon_indices, inxcells*sizeof(int));
      acc_delete(p_output_parent_lat_indices, inxcells*sizeof(int));
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
void check_answers_on_host( Answers *answers, Xgrid_config *xgrid, int myinterp_method )
{

  printf("CHECKING COPIED OUT DATA ON HOST\n");

  //cautious step.  to ensure data is actually copied out
  reset_xgrid_on_host(xgrid, answers, myinterp_method);

  for(int n=0 ; n<NTILES_OUT; n++) {
    acc_update_host( &xgrid[n].nxcells, sizeof(size_t) );
    for(int m=0 ; m<NTILES_IN ; m++) {
      int inxcells = answers[n].nxcells[m];
      acc_update_host( &xgrid[n].per_intile[m].nxcells, sizeof(int));
      acc_copyout( xgrid[n].per_intile[m].input_parent_lon_indices,  inxcells*sizeof(int));
      acc_copyout( xgrid[n].per_intile[m].input_parent_lat_indices,  inxcells*sizeof(int));
      acc_copyout( xgrid[n].per_intile[m].output_parent_lon_indices, inxcells*sizeof(int));
      acc_copyout( xgrid[n].per_intile[m].output_parent_lat_indices, inxcells*sizeof(int));
      acc_copyout( xgrid[n].per_intile[m].xcell_area,  inxcells*sizeof(double));
      if(myinterp_method == myCONSERVE_ORDER2) {
        acc_copyout( xgrid[n].per_intile[m].dcentroid_lon, inxcells*sizeof(double));
        acc_copyout( xgrid[n].per_intile[m].dcentroid_lat, inxcells*sizeof(double));
      }
    }
  }
  acc_delete(xgrid, NTILES_OUT*sizeof(Xgrid_config));

  for( int n=0 ; n<NTILES_OUT; n++ ) {

    printf("checking nxcells for n=%d\n", n);
    if( answers[n].total_nxcells != xgrid[n].nxcells ) {
      printf("ERROR nxcells for n=%d: %zu vs %zu\n", answers[n].total_nxcells, xgrid[n].nxcells);
      error("goodbye!");
    }

    for( int m=0 ; m<NTILES_IN ; m++ ) {
      printf("m=%d ", m);
      int inxcells = xgrid[n].per_intile[m].nxcells;
      printf("nxcells "); check_ianswers(1, answers[n].nxcells+m, &(xgrid[n].per_intile[m].nxcells), ON_HOST);
      printf("input_parent_lon_indices ");
      check_ianswers(inxcells, answers[n].input_parent_lon_indices[m],  xgrid[n].per_intile[m].input_parent_lon_indices,  ON_HOST);
      printf("input_parent_lat_indices ");
      check_ianswers(inxcells, answers[n].input_parent_lat_indices[m],  xgrid[n].per_intile[m].input_parent_lat_indices,  ON_HOST);
      printf("output_parent_lon_indices ");
      check_ianswers(inxcells, answers[n].output_parent_lon_indices[m], xgrid[n].per_intile[m].output_parent_lon_indices, ON_HOST);
      printf("output_parent_lat_indices ");
      check_ianswers(inxcells, answers[n].output_parent_lat_indices[m], xgrid[n].per_intile[m].output_parent_lat_indices, ON_HOST);
      printf("xcell_area ");
      check_ranswers(inxcells, answers[n].xcell_area[m],  xgrid[n].per_intile[m].xcell_area,  ON_HOST);
      if(myinterp_method == myCONSERVE_ORDER2) {
        printf("dcentroid_lon ");
        check_ranswers(inxcells, answers[n].dcentroid_lon[m],  xgrid[n].per_intile[m].dcentroid_lon, ON_HOST);
        printf("dcentroid_lat");
        check_ranswers(inxcells, answers[n].dcentroid_lat[m],  xgrid[n].per_intile[m].dcentroid_lat, ON_HOST);
      }
      printf("\n");
    }
  }

}
//------------------------------------------
//------------------------------------------
void reset_xgrid_on_host( Xgrid_config *xgrid, Answers *answers, int myinterp_method )
{

  for(int n=0 ; n<NTILES_OUT ; n++) {
    xgrid[n].nxcells=0;
    for(int m=0 ; m<NTILES_IN ; m++) {
      int inxcells = answers[n].nxcells[m];
      xgrid[n].per_intile[m].nxcells = -99;
      for( int i=0 ; i<inxcells ; i++) {
        xgrid[n].per_intile[m].input_parent_lon_indices[i]  = -99 ;
        xgrid[n].per_intile[m].input_parent_lat_indices[i]  = -99 ;
        xgrid[n].per_intile[m].output_parent_lon_indices[i] = -99 ;
        xgrid[n].per_intile[m].output_parent_lat_indices[i] = -99 ;
        xgrid[n].per_intile[m].xcell_area[i]  = -99.9 ;
        if(myinterp_method == myCONSERVE_ORDER2) {
          xgrid[n].per_intile[m].dcentroid_lon[i] = -99.9;
          xgrid[n].per_intile[m].dcentroid_lat[i] = -99.9;
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
