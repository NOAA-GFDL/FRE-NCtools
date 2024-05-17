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
  size_t total_nxgrid;
  int itile_nxgrid[NTILES_IN];
  int *i_in[NTILES_IN];
  int *j_in[NTILES_IN];
  int *i_out[NTILES_IN];
  int *j_out[NTILES_IN];
  double *area[NTILES_IN];
  double *di_in[NTILES_IN];
  double *dj_in[NTILES_IN];
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
void read_ianswers( FILE *myfile, int *itile_nxgrid, int **ianswer);
void read_ranswers( FILE *myfile, int *itile_nxgrid, double **ianswer);
void check_answers_on_device(Answers *answers, Interp_config *interp, int myinterp_method);
void check_answers_on_host(Answers *answers, Interp_config *interp, int myinterp_method);
void reset_interp_on_host( Interp_config *interp, Answers *answers, int myinterp_method );
void check_ianswers(int n, int *answers, int *checkme, int host_or_device);
void check_ranswers(int n, double *answers, double *checkme, int host_or_device);
void error(char *error_message);

// start program
int main(int argc, char *argv[]) {

  Interp_config interp[NTILES_OUT];
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
    interp[n].interp_mini = ( Interp_config_mini *)malloc(NTILES_IN*sizeof(Interp_config_mini));
    if(myinterp_method == myCONSERVE_ORDER1) strcpy(interp[n].remap_file, remap_files1[n]);
    if(myinterp_method == myCONSERVE_ORDER2) strcpy(interp[n].remap_file, remap_files2[n]);
    if(access(interp[n].remap_file, F_OK)==0) interp[n].file_exist=1;
  }

  // read in remap and transfer data to device
  read_remap_file_acc(NTILES_IN, NTILES_OUT, interp, opcode) ;
  copy_interp_to_device_acc(NTILES_IN, NTILES_OUT, interp, opcode) ;

  // get all answers
  for(int n=0 ; n<NTILES_OUT ; n++) {
    if(myinterp_method == myCONSERVE_ORDER1) strcpy(answers[n].answer_file, answer_files1[n]);
    if(myinterp_method == myCONSERVE_ORDER2) strcpy(answers[n].answer_file, answer_files2[n]);
  }
  read_all_answers(answers, myinterp_method);


  // check answers on device
  check_answers_on_device(answers, interp, myinterp_method);

  // compare answers on host
  check_answers_on_host( answers, interp, myinterp_method);

  return 0; //test success

}
//------------------------------------------
//------------------------------------------
void read_ianswers( FILE *myfile, int *itile_nxgrid, int **ianswer)
{

  for( int m=0 ; m<NTILES_IN ; m++ ) {
    int inxgrid = itile_nxgrid[m] ;
    for( int j=0 ; j<inxgrid ; j++ ) fscanf( myfile, "%d", ianswer[m]+j );
  }

}
//------------------------------------------
//------------------------------------------
void read_ranswers( FILE *myfile, int *itile_nxgrid, double **ranswer)
{

  for( int m=0 ; m<NTILES_IN ; m++ ) {
    int inxgrid = itile_nxgrid[m] ;
    for( int j=0 ; j<inxgrid ; j++ ) fscanf( myfile, "%lf", ranswer[m]+j );
  }

}
//------------------------------------------
//------------------------------------------
void read_all_answers(Answers *answers, int myinterp_method)
{

  for(int n=0 ; n<NTILES_OUT ; n++) {

    FILE * myfile = fopen(answers[n].answer_file, "r");

    printf("READING IN ANSWERS FOR n=%d\n", n);

    // get answers for interp[n].nxgrid of type size_t
    fscanf( myfile, "%zu", &(answers[n].total_nxgrid) );

    // get answers for interp[n].interp_mini[m].nxgrid
    for( int m=0 ; m<NTILES_IN ; m++ ) fscanf( myfile, "%d", answers[n].itile_nxgrid+m );

    // allocate arrays for answerss
    for( int m=0 ; m<NTILES_IN ; m++ ) {
      int inxgrid = answers[n].itile_nxgrid[m] ;
      answers[n].i_in[m]  = (int *)malloc( inxgrid * sizeof(int) ) ;
      answers[n].j_in[m]  = (int *)malloc( inxgrid * sizeof(int) ) ;
      answers[n].i_out[m] = (int *)malloc( inxgrid * sizeof(int) ) ;
      answers[n].j_out[m] = (int *)malloc( inxgrid * sizeof(int) ) ;
      answers[n].area[m]  = (double *)malloc( inxgrid * sizeof(double) );
      if( myinterp_method == myCONSERVE_ORDER2 ) {
        answers[n].di_in[m] = (double *)malloc( inxgrid * sizeof(double) );
        answers[n].dj_in[m] = (double *)malloc( inxgrid * sizeof(double) );
      }
    }

    //read in answers
    read_ianswers( myfile, answers[n].itile_nxgrid, answers[n].i_in) ;
    read_ianswers( myfile, answers[n].itile_nxgrid, answers[n].j_in);
    read_ianswers( myfile, answers[n].itile_nxgrid, answers[n].i_out);
    read_ianswers( myfile, answers[n].itile_nxgrid, answers[n].j_out);
    read_ranswers( myfile, answers[n].itile_nxgrid, answers[n].area);

    //read di_in dj_in if conserve_order2
    if( myinterp_method == myCONSERVE_ORDER2 ) {
      read_ranswers( myfile, answers[n].itile_nxgrid, answers[n].di_in);
      read_ranswers( myfile, answers[n].itile_nxgrid, answers[n].dj_in);
    }

    fclose( myfile );

  }

}
//------------------------------------------
//------------------------------------------
void check_answers_on_device( Answers *answers, Interp_config *interp, int myinterp_method )
{

  int *p_i_in, *p_j_in, *p_i_out, *p_j_out;
  double *p_area, *p_di_in, *p_dj_in;

  printf("CHECKING ANSWERS ON DEVICE\n");

  if( !acc_is_present(interp, NTILES_OUT*sizeof(Interp_config)) ) error("ERROR interp is not present") ;

  for(int n=0 ; n<NTILES_OUT ; n++) {

    size_t nxgrid = answers[n].total_nxgrid;

    printf("checking for n=%d\n", n);

#pragma acc parallel loop seq copyin(nxgrid) present(interp[n].nxgrid)
    for(int i=0 ; i<1 ; i++ ){
      if( nxgrid != interp[n].nxgrid ) {
        printf("ERROR checking total number of nxgrid for n=0", n);
        nxgrid = 0;
      }
    }
    if( nxgrid == 0 ) error("ERROR with the value of nxgrid on device");

    for( int m=0 ; m<NTILES_IN ; m++ ) {

      int inxgrid = answers[n].itile_nxgrid[m];

      //copy in answers
      acc_copyin(&inxgrid, sizeof(int));
      p_i_in  = answers[n].i_in[m];  acc_copyin(p_i_in, inxgrid*sizeof(int));
      p_j_in  = answers[n].j_in[m];  acc_copyin(p_j_in, inxgrid*sizeof(int));
      p_i_out = answers[n].i_out[m]; acc_copyin(p_i_out, inxgrid*sizeof(int));
      p_j_out = answers[n].j_out[m]; acc_copyin(p_j_out, inxgrid*sizeof(int));
      p_area  = answers[n].area[m];  acc_copyin(p_area, inxgrid*sizeof(double));
      if(myinterp_method == myCONSERVE_ORDER2) {
        p_di_in = answers[n].di_in[m]; acc_copyin(p_di_in, inxgrid*sizeof(double));
        p_dj_in = answers[n].dj_in[m]; acc_copyin(p_dj_in, inxgrid*sizeof(double));
      }

      printf("m=%d ", m);
      printf("nxgrid "); check_ianswers(1, &inxgrid, &(interp[n].interp_mini[m].nxgrid),  ON_DEVICE);
      printf("i_in ");   check_ianswers(inxgrid, p_i_in,  interp[n].interp_mini[m].i_in,  ON_DEVICE);
      printf("j_in ");   check_ianswers(inxgrid, p_j_in,  interp[n].interp_mini[m].j_in,  ON_DEVICE);
      printf("i_out ");  check_ianswers(inxgrid, p_i_out, interp[n].interp_mini[m].i_out, ON_DEVICE);
      printf("j_out ");  check_ianswers(inxgrid, p_j_out, interp[n].interp_mini[m].j_out, ON_DEVICE);
      printf("area ");   check_ranswers(inxgrid, p_area,  interp[n].interp_mini[m].area,  ON_DEVICE);
      if( myinterp_method == myCONSERVE_ORDER2 ){
        printf("di_in "); check_ranswers(inxgrid, p_di_in, interp[n].interp_mini[m].di_in,  ON_DEVICE);
        printf("dj_in "); check_ranswers(inxgrid, p_dj_in, interp[n].interp_mini[m].dj_in,  ON_DEVICE);
      }
      printf("\n");

      acc_delete(&inxgrid, sizeof(int));
      acc_delete(p_i_in, inxgrid*sizeof(int));
      acc_delete(p_j_in, inxgrid*sizeof(int));
      acc_delete(p_i_out, inxgrid*sizeof(int));
      acc_delete(p_j_out, inxgrid*sizeof(int));
      acc_delete(p_area, inxgrid*sizeof(double));
      if(myinterp_method == myCONSERVE_ORDER2) {
        acc_delete(p_di_in, inxgrid*sizeof(double));
        acc_delete(p_dj_in, inxgrid*sizeof(double));
      }
    }
  }

}
//------------------------------------------
//------------------------------------------
void check_answers_on_host( Answers *answers, Interp_config *interp, int myinterp_method )
{

  printf("CHECKING COPIED OUT DATA ON HOST\n");

  //cautious step.  to ensure data is actually copied out
  reset_interp_on_host(interp, answers, myinterp_method);

  for(int n=0 ; n<NTILES_OUT; n++) {
    acc_update_host( &interp[n].nxgrid, sizeof(size_t) );
    for(int m=0 ; m<NTILES_IN ; m++) {
      int inxgrid = answers[n].itile_nxgrid[m];
      acc_update_host( &interp[n].interp_mini[m].nxgrid, sizeof(int));
      acc_copyout( interp[n].interp_mini[m].i_in,  inxgrid*sizeof(int));
      acc_copyout( interp[n].interp_mini[m].j_in,  inxgrid*sizeof(int));
      acc_copyout( interp[n].interp_mini[m].i_out, inxgrid*sizeof(int));
      acc_copyout( interp[n].interp_mini[m].j_out, inxgrid*sizeof(int));
      acc_copyout( interp[n].interp_mini[m].area,  inxgrid*sizeof(double));
      if(myinterp_method == myCONSERVE_ORDER2) {
        acc_copyout( interp[n].interp_mini[m].di_in, inxgrid*sizeof(double));
        acc_copyout( interp[n].interp_mini[m].dj_in, inxgrid*sizeof(double));
      }
    }
  }
  acc_delete(interp, NTILES_OUT*sizeof(Interp_config));

  for( int n=0 ; n<NTILES_OUT; n++ ) {

    printf("checking nxgrid for n=%d\n", n);
    if( answers[n].total_nxgrid != interp[n].nxgrid ) {
      printf("ERROR nxgrid for n=%d: %zu vs %zu\n", answers[n].total_nxgrid, interp[n].nxgrid);
      error("goodbye!");
    }

    for( int m=0 ; m<NTILES_IN ; m++ ) {
      printf("m=%d ", m);
      int inxgrid = interp[n].interp_mini[m].nxgrid;
      printf("nxgrid "); check_ianswers(1, answers[n].itile_nxgrid+m, &(interp[n].interp_mini[m].nxgrid), ON_HOST);
      printf("i_in ");   check_ianswers(inxgrid, answers[n].i_in[m],  interp[n].interp_mini[m].i_in,  ON_HOST);
      printf("j_in ");   check_ianswers(inxgrid, answers[n].j_in[m],  interp[n].interp_mini[m].j_in,  ON_HOST);
      printf("i_out ");  check_ianswers(inxgrid, answers[n].i_out[m], interp[n].interp_mini[m].i_out, ON_HOST);
      printf("j_out ");  check_ianswers(inxgrid, answers[n].j_out[m], interp[n].interp_mini[m].j_out, ON_HOST);
      printf("area ");   check_ranswers(inxgrid, answers[n].area[m],  interp[n].interp_mini[m].area,  ON_HOST);
      if(myinterp_method == myCONSERVE_ORDER2) {
        printf("di_in "); check_ranswers(inxgrid, answers[n].di_in[m],  interp[n].interp_mini[m].di_in, ON_HOST);
        printf("dj_in");  check_ranswers(inxgrid, answers[n].dj_in[m],  interp[n].interp_mini[m].dj_in, ON_HOST);
      }
      printf("\n");
    }
  }

}
//------------------------------------------
//------------------------------------------
void reset_interp_on_host( Interp_config *interp, Answers *answers, int myinterp_method )
{

  for(int n=0 ; n<NTILES_OUT ; n++) {
    interp[n].nxgrid=0;
    for(int m=0 ; m<NTILES_IN ; m++) {
      int inxgrid = answers[n].itile_nxgrid[m];
      interp[n].interp_mini[m].nxgrid = -99;
      for( int i=0 ; i<inxgrid ; i++) {
        interp[n].interp_mini[m].i_in[i]  = -99 ;
        interp[n].interp_mini[m].j_in[i]  = -99 ;
        interp[n].interp_mini[m].i_out[i] = -99 ;
        interp[n].interp_mini[m].j_out[i] = -99 ;
        interp[n].interp_mini[m].area[i]  = -99.9 ;
        if(myinterp_method == myCONSERVE_ORDER2) {
          interp[n].interp_mini[m].di_in[i] = -99.9;
          interp[n].interp_mini[m].dj_in[i] = -99.9;
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
