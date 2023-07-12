/***********************************************************************
 *                   GNU Lesser General Public License
 *
 * This file is part of the GFDL FRE NetCDF tools package (FRE-NCTools).
 *
 * FRE-NCTools is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at
 * your option) any later version.
 *
 * FRE-NCTools is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with FRE-NCTools.  If not, see
 * <http://www.gnu.org/licenses/>.
 **********************************************************************/
/* The following is a test program to test subroutines in create_xgrid.c */

#include <math.h>
#include <string.h>
#include "create_xgrid.h"
#include "mosaic_util.h"
#include "constant.h"

#define D2R (M_PI/180)
#define R2D (180/M_PI)
#define MAXPOINT 1000

int main(int argc, char* argv[])
{

  double lon1_in[MAXPOINT], lat1_in[MAXPOINT];
  double lon2_in[MAXPOINT], lat2_in[MAXPOINT];
  double x1_in[MAXPOINT], y1_in[MAXPOINT], z1_in[MAXPOINT];
  double x2_in[MAXPOINT], y2_in[MAXPOINT], z2_in[MAXPOINT];
  double lon_out[20], lat_out[20];
  double x_out[20], y_out[20], z_out[20];
  int    n1_in, n2_in, n_out, i, j;
  int    nlon1=0, nlat1=0, nlon2=0, nlat2=0;
  int    n;
  int    ntest = 26;


  for(n=11; n<=ntest; n++) {

    switch (n) {
    case 1:
      /****************************************************************

       test clip_2dx2d_great_cirle case 1:
       box 1: (20,10), (20,12), (22,12), (22,10)
       box 2: (21,11), (21,14), (24,14), (24,11)
       out  : (21, 12.0018), (22, 12), (22, 11.0033), (21, 11)

      ****************************************************************/
      n1_in = 4; n2_in = 4;
      /* first a simple lat-lon grid box to clip another lat-lon grid box */
      lon1_in[0] = 20; lat1_in[0] = 10;
      lon1_in[1] = 20; lat1_in[1] = 12;
      lon1_in[2] = 22; lat1_in[2] = 12;
      lon1_in[3] = 22; lat1_in[3] = 10;
      lon2_in[0] = 21; lat2_in[0] = 11;
      lon2_in[1] = 21; lat2_in[1] = 14;
      lon2_in[2] = 24; lat2_in[2] = 14;
      lon2_in[3] = 24; lat2_in[3] = 11;
      break;

    case 2:
      /****************************************************************

        test clip_2dx2d_great_cirle case 2: two identical box
        box 1: (20,10), (20,12), (22,12), (22,10)
        box 2: (20,10), (20,12), (22,12), (22,10)
        out  : (20,10), (20,12), (22,12), (22,10)

      ****************************************************************/
      n1_in = 4; n2_in = 4;
      lon1_in[0] = 20; lat1_in[0] =-90;
      lon1_in[1] = 20; lat1_in[1] =-89;
      lon1_in[2] = 22; lat1_in[2] =-88;
      lon1_in[3] = 24; lat1_in[3] =-89;

      for(i=0; i<n2_in; i++) {
	lon2_in[i] = lon1_in[i];
	lat2_in[i] = lat1_in[i];
      }
      break;

    case 3:
      /****************************************************************

       test clip_2dx2d_great_cirle case 3: one cubic sphere grid close to the pole with lat-lon grid.
       box 1: (251.7, 88.98), (148.3, 88.98), (57.81, 88.72), (342.2, 88.72)
       box 2: (150, 88), (150, 90), (152.5, 90), (152.5, 88)
       out  : (152.5, 89.0642), (150, 89.0165), (0, 90)

      ****************************************************************/
      n1_in = 4; n2_in = 4;
      /* first a simple lat-lon grid box to clip another lat-lon grid box */
      lon1_in[0] = 251.7; lat1_in[0] = 88.98;
      lon1_in[1] = 148.3; lat1_in[1] = 88.98;
      lon1_in[2] = 57.81; lat1_in[2] = 88.72;
      lon1_in[3] = 342.2; lat1_in[3] = 88.72;

      lon2_in[0] = 150; lat2_in[0] = 88;
      lon2_in[1] = 150; lat2_in[1] = 90;
      lon2_in[2] = 152.5; lat2_in[2] = 90;
      lon2_in[3] = 152.5; lat2_in[3] = 88;
      /*
      for(i=0; i<4; i++) {
	lon2_in[i] = lon1_in[i];
	lat2_in[i] = lat1_in[i];
      }
      */
      break;

    case 4:
      /****************************************************************

       test clip_2dx2d_great_cirle case 4: One box contains the pole
       box 1: (-160, 88.5354), (152.011, 87.8123) , (102.985, 88.4008), (20, 89.8047)
       box 2: (145,88), (145,90), (150,90), (150,88)
       out  : (145.916, 88.0011), (145, 88.0249), (0, 90), (150, 88)

      ****************************************************************/
      n1_in = 4; n2_in = 4;
      /* first a simple lat-lon grid box to clip another lat-lon grid box */

      lon1_in[0] = -160;  lat1_in[0] = 88.5354;
      lon1_in[1] = 152.011; lat1_in[1] = 87.8123;
      lon1_in[2] = 102.985; lat1_in[2] = 88.4008;
      lon1_in[3] = 20; lat1_in[3] = 89.8047;

      lon2_in[0] = 145; lat2_in[0] = 88;
      lon2_in[1] = 145; lat2_in[1] = 90;
      lon2_in[2] = 150; lat2_in[2] = 90;
      lon2_in[3] = 150; lat2_in[3] = 88;
      break;

    case 5:
      /****************************************************************

       test clip_2dx2d_great_cirle case 5: One tripolar grid around the pole with lat-lon grid.
       box 1: (-202.6, 87.95), (-280, 89.56), (-100, 90), (-190, 88)
       box 2: (21,11), (21,14), (24,14), (24,11)
       out  : (150, 88.7006), (145,  88.9507), (0, 90)

      ****************************************************************/
      n1_in = 4; n2_in = 4;
      /* first a simple lat-lon grid box to clip another lat-lon grid box */

      lon1_in[0] = -202.6;  lat1_in[0] = 87.95;
      lon1_in[1] = -280.;   lat1_in[1] = 89.56;
      lon1_in[2] = -100.0; lat1_in[2] = 90;
      lon1_in[3] = -190.; lat1_in[3] = 88;

      lon2_in[0] = 145; lat2_in[0] = 88;
      lon2_in[1] = 145; lat2_in[1] = 90;
      lon2_in[2] = 150; lat2_in[2] = 90;
      lon2_in[3] = 150; lat2_in[3] = 88;
      break;

    case 6:
      /****************************************************************

       test clip_2dx2d_great_cirle case 6: One cubic sphere grid arounc the pole with one tripolar grid box
                                       around the pole.
       box 1: (-160, 88.5354), (152.011, 87.8123) , (102.985, 88.4008), (20, 89.8047)
       box 2: (-202.6, 87.95), (-280, 89.56), (-100, 90), (-190, 88)
       out  : (170, 88.309), (157.082, 88.0005), (83.714, 89.559), (80, 89.6094), (0, 90), (200, 88.5354)


      ****************************************************************/
      n1_in = 4; n2_in = 4;
      /* first a simple lat-lon grid box to clip another lat-lon grid box */

      lon1_in[0] = -160;  lat1_in[0] = 88.5354;
      lon1_in[1] = 152.011; lat1_in[1] = 87.8123;
      lon1_in[2] = 102.985; lat1_in[2] = 88.4008;
      lon1_in[3] = 20; lat1_in[3] = 89.8047;

      lon2_in[0] = -202.6;  lat2_in[0] = 87.95;
      lon2_in[1] = -280.;   lat2_in[1] = 89.56;
      lon2_in[2] = -100.0;  lat2_in[2] = 90;
      lon2_in[3] = -190.;   lat2_in[3] = 88;
      break;

    case 7:
      /****************************************************************

       test clip_2dx2d_great_cirle case 7: One small grid box inside a big grid box.
       box 1: (20,10), (20,12), (22,12), (22,10)
       box 2: (18,8), (18,14), (24,14), (24,8)
       out  : (20,10), (20,12), (22,12), (22,10)

      ****************************************************************/
      n1_in = 4; n2_in = 4;
      /* first a simple lat-lon grid box to clip another lat-lon grid box */
      lon1_in[0] = 20; lat1_in[0] = 10;
      lon1_in[1] = 20; lat1_in[1] = 12;
      lon1_in[2] = 22; lat1_in[2] = 12;
      lon1_in[3] = 22; lat1_in[3] = 10;
      lon2_in[0] = 18; lat2_in[0] = 8;
      lon2_in[1] = 18; lat2_in[1] = 14;
      lon2_in[2] = 24; lat2_in[2] = 14;
      lon2_in[3] = 24; lat2_in[3] = 8;
      break;

    case 8:
      /****************************************************************

       test clip_2dx2d_great_cirle case 8: Cubic sphere grid at tile = 1, point (i=25,j=1)
          with N45 at (i=141,j=23)
       box 1:
       box 2:
       out  : None

      ****************************************************************/
      n1_in = 4; n2_in = 4;
      /* first a simple lat-lo
	 n grid box to clip another lat-lon grid box */
      lon1_in[0] = 350.0; lat1_in[0] = -45;
      lon1_in[1] = 350.0; lat1_in[1] = -43.43;
      lon1_in[2] = 352.1; lat1_in[2] = -43.41;
      lon1_in[3] = 352.1; lat1_in[3] = -44.98;
      lon2_in[0] = 350.0;   lat2_in[0] = -46;
      lon2_in[1] = 350.0;   lat2_in[1] = -44;
      lon2_in[2] = 352.5; lat2_in[2] = -44;
      lon2_in[3] = 352.5; lat2_in[3] = -46;
      break;

    case 9:
      /****************************************************************

       test clip_2dx2d_great_cirle case 9: Cubic sphere grid at tile = 1, point (i=1,j=1)
          with N45 at (i=51,j=61)
       box 1:
       box 2:
       out  : None

      ****************************************************************/
      n1_in = 4; n2_in = 4;

      lon1_in[0] = 305.0; lat1_in[0] = -35.26;
      lon1_in[1] = 305.0; lat1_in[1] = -33.80;
      lon1_in[2] = 306.6; lat1_in[2] = -34.51;
      lon1_in[3] = 306.6; lat1_in[3] = -35.99;
      lon2_in[0] = 125;   lat2_in[0] = 32;
      lon2_in[1] = 125;   lat2_in[1] = 34;
      lon2_in[2] = 127.5; lat2_in[2] = 34;
      lon2_in[3] = 127.5; lat2_in[3] = 32;
      break;

    case 10:
      /****************************************************************

       test clip_2dx2d_great_cirle case 10: Cubic sphere grid at tile = 3, point (i=24,j=1)
          with N45 at (i=51,j=46)
       box 1:
       box 2:
       out  : None

      ****************************************************************/
      n1_in = 4; n2_in = 4;

      lon1_in[0] = 125.0; lat1_in[0] = 1.46935;
      lon1_in[1] = 126.573; lat1_in[1] = 1.5091;
      lon1_in[2] = 126.573; lat1_in[2] = 0;
      lon1_in[3] = 125.0; lat1_in[3] = 0;
      lon2_in[0] = 125;   lat2_in[0] = 0;
      lon2_in[1] = 125;   lat2_in[1] = 2;
      lon2_in[2] = 127.5; lat2_in[2] = 2;
      lon2_in[3] = 127.5; lat2_in[3] = 0;
      break;

    case 11:
      /****************************************************************

       test clip_2dx2d_great_cirle case 10: Cubic sphere grid at tile = 3, point (i=24,j=1)
          with N45 at (i=51,j=46)
       box 1:
       box 2:
       out  :

      ****************************************************************/
      nlon1 = 1;
      nlat1 = 1;
      nlon2 = 1;
      nlat2 = 1;
      n1_in = (nlon1+1)*(nlat1+1);
      n2_in = (nlon2+1)*(nlat2+1);

      lon1_in[0] = 350.0; lat1_in[0] = 90.00;
      lon1_in[1] = 170.0; lat1_in[1] = 87.92;
      lon1_in[2] = 260.0; lat1_in[2] = 87.92;
      lon1_in[3] = 215.0;  lat1_in[3] = 87.06;

/*       lon1_in[0] = 35.0; lat1_in[0] = 87.06; */
/*       lon1_in[1] = 80.0; lat1_in[1] = 87.92; */
/*       lon1_in[2] = 125.0; lat1_in[2] = 87.06; */
/*       lon1_in[3] = 350.0; lat1_in[3] = 87.92; */
/*       lon1_in[4] = 350.0; lat1_in[4] = 90.00; */
/*       lon1_in[5] = 170.0; lat1_in[5] = 87.92; */
/*       lon1_in[6] = 305.0; lat1_in[6] = 87.06; */
/*       lon1_in[7] = 260.0; lat1_in[7] = 87.92; */
/*       lon1_in[8] = 215.0;  lat1_in[8] = 87.06; */

      lon2_in[0] = 167.5; lat2_in[0] = 88;
      lon2_in[1] = 170;   lat2_in[1] = 88;
      lon2_in[2] = 167.5; lat2_in[2] = 90;
      lon2_in[3] = 170;   lat2_in[3] = 90;

/*       nlon1 = 3; */
/*       nlat1 = 2; */
/*       nlon2 = 1; */
/*       nlat2 = 1; */
/*       n1_in = (nlon1+1)*(nlat1+1); */
/*       n2_in = (nlon2+1)*(nlat2+1); */

/*       lon1_in[0] = 35.00;     lat1_in[0] = -59.90; */
/*       lon1_in[1] = 37.64;     lat1_in[1] = -58.69; */
/*       lon1_in[2] = 40.07;     lat1_in[2] = -57.44; */
/*       lon1_in[3] = 42.32;     lat1_in[3] = -56.15; */
/*       lon1_in[4] = 32.36;     lat1_in[4] = -58.69; */
/*       lon1_in[5] = 35.00;     lat1_in[5] = -57.56; */
/*       lon1_in[6] = 37.45;     lat1_in[6] = -56.39; */
/*       lon1_in[7] = 39.74;     lat1_in[7] = -55.18; */
/*       lon1_in[8] = 29.93;     lat1_in[8] = -57.44; */
/*       lon1_in[9] = 32.55;     lat1_in[9] = -56.39; */
/*       lon1_in[10] = 35.00;     lat1_in[10] = -55.29; */
/*       lon1_in[11] = 37.30;     lat1_in[11] = -54.16; */
/*       lon2_in[0] = 35;   lat2_in[0] = -58; */
/*       lon2_in[1] = 37.5; lat2_in[1] = -58; */
/*       lon2_in[2] = 35;   lat2_in[2] = -56; */
/*       lon2_in[3] = 37.5; lat2_in[3] = -56; */

/*       nlon1 = 1; */
/*       nlat1 = 1; */
/*       nlon2 = 1; */
/*       nlat2 = 1; */
/*       n1_in = (nlon1+1)*(nlat1+1); */
/*       n2_in = (nlon2+1)*(nlat2+1); */

/*       lon1_in[0] = 305;     lat1_in[0] = -35.26; */
/*       lon1_in[1] = 306;     lat1_in[1] = -35.99; */
/*       lon1_in[2] = 305;     lat1_in[2] = -33.80; */
/*       lon1_in[3] = 306;     lat1_in[3] = -34.51; */
/*       lon2_in[0] = 305;   lat2_in[0] = -34; */
/*       lon2_in[1] = 307.5; lat2_in[1] = -34; */
/*       lon2_in[2] = 305;   lat2_in[2] = -32; */
/*       lon2_in[3] = 307.5; lat2_in[3] = -32; */

       nlon1 = 2;
       nlat1 = 2;
       nlon2 = 1;
       nlat2 = 1;
      n1_in = (nlon1+1)*(nlat1+1);
      n2_in = (nlon2+1)*(nlat2+1);

      lon1_in[0] = 111.3; lat1_in[0] = 1.591;
      lon1_in[1] = 109.7; lat1_in[1] = 2.926;
      lon1_in[2] = 108.2; lat1_in[2] = 4.256;
      lon1_in[3] = 110.0; lat1_in[3] = 0.000;
      lon1_in[4] = 108.4; lat1_in[4] = 1.335;
      lon1_in[5] = 106.8; lat1_in[5] = 2.668;
      lon1_in[6] = 108.7; lat1_in[6] = -1.591;
      lon1_in[7] = 107.1; lat1_in[7] = -0.256;
      lon1_in[8] = 105.5;  lat1_in[8] = 1.078;

      lon2_in[0] = 107.5; lat2_in[0] = 0;
      lon2_in[1] = 110;   lat2_in[1] = 0;
      lon2_in[2] = 107.5; lat2_in[2] = 2;
      lon2_in[3] = 110;   lat2_in[3] = 2;

      break;

    case 12:
      /****************************************************************

       test : create_xgrid_great_circle
       box 1: (20,10), (20,12), (22,12), (22,10)
       box 2: (21,11), (21,14), (24,14), (24,11)
       out  : (21, 12.0018), (22, 12), (22, 11.0033), (21, 11)

      ****************************************************************/
      nlon1 = 2;
      nlat1 = 2;
      nlon2 = 3;
      nlat2 = 3;
      n1_in = (nlon1+1)*(nlat1+1);
      n2_in = (nlon2+1)*(nlat2+1);

      /* first a simple lat-lon grid box to clip another lat-lon grid box */
      for(j=0; j<=nlat1; j++) for(i=0; i<=nlon1; i++){
	lon1_in[j*(nlon1+1)+i] = 20.0 + (i-1)*2.0;
	lat1_in[j*(nlon1+1)+i] = 10.0 + (j-1)*2.0;
      }
       for(j=0; j<=nlat2; j++) for(i=0; i<=nlon2; i++){
	lon2_in[j*(nlon2+1)+i] = 19.0 + (i-1)*2.0;
	lat2_in[j*(nlon2+1)+i] = 9.0 + (j-1)*2.0;
      }

      break;

    case 13:

      nlon1 = 1;
      nlat1 = 1;
      nlon2 = 1;
      nlat2 = 1;
      n1_in = (nlon1+1)*(nlat1+1);
      n2_in = (nlon2+1)*(nlat2+1);

/*       lon1_in[0] = ; lat1_in[0] = ; */
/*       lon1_in[1] = ; lat1_in[1] = ; */
/*       lon1_in[2] = ; lat1_in[2] = ; */
/*       lon1_in[3] = ; lat1_in[3] = ; */
/*       lon2_in[0] = ; lat2_in[0] = ; */
/*       lon2_in[1] = ; lat2_in[1] = ; */
/*       lon2_in[2] = ; lat2_in[2] = ; */
/*       lon2_in[3] = ; lat2_in[3] = ;     */

/*       lon1_in[0] = 1.35536; lat1_in[0] = 1.16251; */
/*       lon1_in[1] = 1.36805; lat1_in[1] = 1.15369; */
/*       lon1_in[2] = 1.37843; lat1_in[2] = 1.16729; */
/*       lon1_in[3] = 1.39048; lat1_in[3] = 1.15826; */
/*       lon2_in[0] = 1.34611; lat2_in[0] = 1.16372; */
/*       lon2_in[1] = 1.35616; lat2_in[1] = 1.15802;    */
/*       lon2_in[2] = 1.35143; lat2_in[2] = 1.16509; */
/*       lon2_in[3] = 1.36042; lat2_in[3] = 1.15913; */

/*       lon1_in[0] = 12.508065121288551; lat1_in[0] = -87.445883646793547; */
/*       lon1_in[1] = 325.425637772; lat1_in[1] = -86.481216821859505; */
/*       lon1_in[2] = 97.5; lat1_in[2] = -89.802136057677174; */
/*       lon1_in[3] = 277.5; lat1_in[3] = -87.615232005344637; */

/*       for(j=0; j<=nlat2; j++) for(i=0; i<=nlon2; i++) { */
/* 	lon2_in[j*(nlon2+1)+i] = -280.0 + i*1.0; */
/* 	lat2_in[j*(nlon2+1)+i] = -90.0 + j*8.0; */
/*       } */
      lon1_in[0] = 120.369397984526174; lat1_in[0] = 16.751543427495864;
      lon1_in[1] = 119.999999999999986; lat1_in[1] = 16.751871929590038;
      lon1_in[2] = 120.369397846883501; lat1_in[2] = 16.397797979598028;
      lon1_in[3] = 119.999999999999986; lat1_in[3] = 16.398120477217255;
      lon2_in[0] = 120.369415056522087; lat2_in[0] = 16.752176828509153;
      lon2_in[1] = 119.999999999999986; lat2_in[1] = 16.752505523196167;
      lon2_in[2] = 120.369415056522087; lat2_in[2] = 16.397797949548146;
      lon2_in[3] = 119.999999999999986; lat2_in[3] = 16.398120477217255;

      break;

    case 14:
      /****************************************************************
       test clip_2dx2d_great_cirle case 14: Cubic sphere grid at tile = 3, point (i=24,j=1)
         identical grid boxes
      ****************************************************************/
      /*
      nlon1 = 1;
      nlat1 = 1;
      nlon2 = 1;
      nlat2 = 1;
      n1_in = (nlon1+1)*(nlat1+1);
      n2_in = (nlon2+1)*(nlat2+1);

      lon1_in[0] = 350.0; lat1_in[0] = 90.00;
      lon1_in[1] = 170.0; lat1_in[1] = 87.92;
      lon1_in[2] = 260.0; lat1_in[2] = 87.92;
      lon1_in[3] = 215.0; lat1_in[3] = 87.06;

      lon2_in[0] = 350.0; lat2_in[0] = 90.00;
      lon2_in[1] = 170.0; lat2_in[1] = 87.92;
      lon2_in[2] = 260.0; lat2_in[2] = 87.92;
      lon2_in[3] = 215.0; lat2_in[3] = 87.06;
      */
      n1_in = 4; n2_in = 4;

      //double lon1_14[] = {82.400,82.400,262.400,262.400,326.498,379.641};
      //double lat1_14[] = {89.835,90.000, 90.000, 89.847, 89.648, 89.642};
      double lon1_14[] = {350.,170.,260.,215.};
      double lat1_14[] = {90.,87.92,87.92,87.06};
      //double lon1_14[] = {82.400,262.400,326.498};
      //double lat1_14[] = {89.835, 90.000, 89.648};
      memcpy(lon1_in,lon1_14,sizeof(lon1_in));
      memcpy(lat1_in,lat1_14,sizeof(lat1_in));
      memcpy(lon2_in,lon1_14,sizeof(lon2_in));
      memcpy(lat2_in,lat1_14,sizeof(lat2_in));
      break;

    case 15:
      n1_in = 6; n2_in = 5;
      double lon1_15[] = {145.159, 198.302, 262.400, 262.400,  82.400,  82.400};
      double lat1_15[] = { 89.642,  89.648,  89.847,  90.000,  90.000,  89.835};

      double lon2_15[] = {150.000, 177.824, 240.000, 240.000, 150.000};
      double lat2_15[] = { 89.789,  89.761,  89.889,  90.000,  90.000};
      memcpy(lon1_in,lon1_15,sizeof(lon1_in));
      memcpy(lat1_in,lat1_15,sizeof(lat1_in));
      memcpy(lon2_in,lon2_15,sizeof(lon2_in));
      memcpy(lat2_in,lat2_15,sizeof(lat2_in));
      //Must give the second box
      break;

    case 16:
      /*Must give [[-57.748, -30, -30, -97.6, -97.6],
                   [89.876, 89.891, 90, 90, 89.9183]]*/
      n1_in = 6; n2_in = 5;
      double lon1_16[] = {82.400,  82.400, 262.400, 262.400, 326.498, 379.641};
      double lat1_16[] = {89.835,  90.000,  90.000,  89.847,  89.648,  89.642};
      double lon2_16[] = {302.252, 330.000, 330.000, 240.000, 240.000};
      double lat2_16[] = {89.876,  89.891,  90.000,  90.000,  89.942};
      memcpy(lon1_in,lon1_16,sizeof(lon1_in));
      memcpy(lat1_in,lat1_16,sizeof(lat1_in));
      memcpy(lon2_in,lon2_16,sizeof(lon2_in));
      memcpy(lat2_in,lat2_16,sizeof(lat2_in));
      break;

    case 17:
      /*Must give the second square
	 -30, -2.252, 60, 60, -30,
	 89.891, 89.876, 89.942, 90, 90,
      */
      n1_in = 6; n2_in = 5;
      lon1_in[0]=82.400;  lon1_in[1]=82.400; lon1_in[2]=262.400; lon1_in[3]=262.400; lon1_in[4]=326.498; lon1_in[5]=379.641;
      lat1_in[0]=89.835;  lat1_in[1]=90.000; lat1_in[2]=90.000;  lat1_in[3]=89.847;  lat1_in[4]=89.648;  lat1_in[5]=89.642;

      double lon2_17[] = {-30.000,  -2.252,  60.000,  60.000, -30.000};
      double lat2_17[] = { 89.891,  89.876,  89.942,  90.000,  90.000};
      memcpy(lon2_in,lon2_17,sizeof(lon2_in));
      memcpy(lat2_in,lat2_17,sizeof(lat2_in));
      break;

    case 18:
      n1_in = 6; n2_in = 5;
      double lon1_18[] = {82.400,  82.400, 262.400, 262.400, 326.498, 379.641};
      double lat1_18[] = {89.835,  90.000,  90.000,  89.847,  89.648,  89.642};

      double lon2_18[] = {150.000, 177.824, 240.000, 240.000, 150.000};
      double lat2_18[] = {89.789,  89.761,  89.889,  90.000,  90.000};
      memcpy(lon1_in,lon1_18,sizeof(lon1_in));
      memcpy(lat1_in,lat1_18,sizeof(lat1_in));
      memcpy(lon2_in,lon2_18,sizeof(lon2_in));
      memcpy(lat2_in,lat2_18,sizeof(lat2_in));
      break;
      /*Must give nothing*/
    case 19:
      /****************************************************************
        test clip_2dx2d 2: two boxes that include the North Pole
                           one has vertices on the tripolar fold
                           the other is totally outside the first
                           This actually happens for some stretched grid
                           configurations  mosaic_c256r25tlat32.0_om4p25
        The test gives wrong answers!
      ****************************************************************/
      n1_in = 6; n2_in = 5;
      double lon1_19[] = {145.159, 198.302, 262.400, 262.400,  82.400,  82.400};
      double lat1_19[] = {89.642,  89.648,  89.847,  90.000,  90.000,  89.835};

      double lon2_19[] = {-30.000,  -2.176,  60.000,  60.000, -30.000};
      double lat2_19[] = {89.789,  89.761,  89.889,  90.000,  90.000};
      memcpy(lon1_in,lon1_19,sizeof(lon1_in));
      memcpy(lat1_in,lat1_19,sizeof(lat1_in));
      memcpy(lon2_in,lon2_19,sizeof(lon2_in));
      memcpy(lat2_in,lat2_19,sizeof(lat2_in));
      break;

    case 20:
      /*Must give
 n_out= 5
 122.176, 150, 150, 82.4, 82.4,
 89.761, 89.789, 90, 90, 89.8429,
       */      n1_in = 6; n2_in = 5;
      double lon1_20[] = {145.159, 198.302, 262.400, 262.400,  82.400,  82.400};
      double lat1_20[] = {89.642,  89.648,  89.847,  90.000,  90.000,  89.835};

      double lon2_20[] = {122.176, 150.000, 150.000,  60.000,  60.000};
      double lat2_20[] = { 89.761,  89.789,  90.000,  90.000,  89.889};
      memcpy(lon1_in,lon1_20,sizeof(lon1_in));
      memcpy(lat1_in,lat1_20,sizeof(lat1_in));
      memcpy(lon2_in,lon2_20,sizeof(lon2_in));
      memcpy(lat2_in,lat2_20,sizeof(lat2_in));
      break;

    case 21:
      /*Must give
 n_out= 5
 60.000,  82.400,  82.400,  60.000],
 89.889,  89.843,  90.000,  90.000]
       */
      n1_in = 6; n2_in = 5;
      double lon1_21[] = {82.400,  82.400, 262.400, 262.400, 326.498, 379.641};
      double lat1_21[] = {89.835,  90.000,  90.000,  89.847,  89.648,  89.642};

      double lon2_21[] = {122.176, 150.000, 150.000,  60.000,  60.000};
      double lat2_21[] = { 89.761,  89.789,  90.000,  90.000,  89.889};
      memcpy(lon1_in,lon1_21,sizeof(lon1_in));
      memcpy(lat1_in,lat1_21,sizeof(lat1_in));
      memcpy(lon2_in,lon2_21,sizeof(lon2_in));
      memcpy(lat2_in,lat2_21,sizeof(lat2_in));
      break;

    case 26:
      /*Side crosses SP (Right cell).
	Must give same box
      */
      n1_in = 4; n2_in = 4;
      double lon1_22[] = {209.68793552504,158.60256162113,82.40000000000,262.40000000000};
      double lat1_22[] = {-89.11514201451,-89.26896927380,-89.82370183256, -89.46584623220};

      double lon2_22[] = {209.68793552504,158.60256162113,82.40000000000,262.40000000000};
      double lat2_22[] = {-89.11514201451,-89.26896927380,-89.82370183256, -89.46584623220};
      memcpy(lon1_in,lon1_22,sizeof(lon1_in));
      memcpy(lat1_in,lat1_22,sizeof(lat1_in));
      memcpy(lon2_in,lon2_22,sizeof(lon2_in));
      memcpy(lat2_in,lat2_22,sizeof(lat2_in));
      break;

    case 23:
      /*Side does not cross SP (Right cell).
	Must give same box
      */

      n1_in = 4; n2_in = 4;
      double lon1_23[] = {158.60256162113,121.19651597620,82.40000000000,82.40000000000};
      double lat1_23[] = {-89.26896927380,-88.85737639760,-89.10746816044,-89.82370183256};

      double lon2_23[] = {158.60256162113,121.19651597620,82.40000000000,82.40000000000};
      double lat2_23[] = {-89.26896927380,-88.85737639760,-89.10746816044,-89.82370183256};
      memcpy(lon1_in,lon1_23,sizeof(lon1_in));
      memcpy(lat1_in,lat1_23,sizeof(lat1_in));
      memcpy(lon2_in,lon2_23,sizeof(lon2_in));
      memcpy(lat2_in,lat2_23,sizeof(lat2_in));
      break;

    case 24:
      /*Side crosses SP (Left cell). Added twin poles.
	Must give the same box
      */
      n1_in = 6; n2_in = 6;
      double lon1_24[] = {262.40000000000,262.40000000000,82.4,82.4,6.19743837887,-44.88793552504};
      double lat1_24[] = {-89.46584623220,-90.0,         -90.0,-89.82370183256, -89.26896927380, -89.11514201451};

      double lon2_24[] = {262.40000000000,262.40000000000,82.4,82.4,6.19743837887,-44.88793552504};
      double lat2_24[] = {-89.46584623220,-90.0,         -90.0,-89.82370183256, -89.26896927380, -89.11514201451};
      memcpy(lon1_in,lon1_24,sizeof(lon1_in));
      memcpy(lat1_in,lat1_24,sizeof(lat1_in));
      memcpy(lon2_in,lon2_24,sizeof(lon2_in));
      memcpy(lat2_in,lat2_24,sizeof(lat2_in));
      break;
    case 25:
      /*Side crosses SP (Left cell).
	Must givethe same box
      */
      n1_in = 4; n2_in = 4;
      double lon1_25[] = {262.40000000000,82.4,6.19743837887,-44.88793552504};
      double lat1_25[] = {-89.46584623220, -89.82370183256, -89.26896927380, -89.11514201451};

      double lon2_25[] = {262.40000000000,82.4,6.19743837887,-44.88793552504};
      double lat2_25[] = {-89.46584623220, -89.82370183256, -89.26896927380, -89.11514201451};
      memcpy(lon1_in,lon1_25,sizeof(lon1_in));
      memcpy(lat1_in,lat1_25,sizeof(lat1_in));
      memcpy(lon2_in,lon2_25,sizeof(lon2_in));
      memcpy(lat2_in,lat2_25,sizeof(lat2_in));
      break;
    case 22:
      /*Side does not cross SP (Left cell).
	Must give same box
      */
      n1_in = 4; n2_in = 4;
      double lon1_26[] = {82.4,82.4,43.60348402380,6.19743837887};
      double lat1_26[] = {-89.82370183256, -89.10746816044, -88.85737639760, -89.26896927380};

      double lon2_26[] = {82.4,82.4,43.60348402380,6.19743837887};
      double lat2_26[] = {-89.82370183256, -89.10746816044, -88.85737639760, -89.26896927380};
      memcpy(lon1_in,lon1_26,sizeof(lon1_in));
      memcpy(lat1_in,lat1_26,sizeof(lat1_in));
      memcpy(lon2_in,lon2_26,sizeof(lon2_in));
      memcpy(lat2_in,lat2_26,sizeof(lat2_in));
      break;
    default:
      error_handler("test_create_xgrid: incorrect case number");
    }

    /* convert to radian */

    for(i=0; i<n1_in; i++) {
      lon1_in[i] *= D2R; lat1_in[i] *=D2R;
    }
    for(i=0; i<n2_in; i++) {
      lon2_in[i] *= D2R; lat2_in[i] *=D2R;
    }


    printf("\n*********************************************************\n");
    printf("               Case %d                                    \n", n);


    if( n > 10 && n <= 14) {
      int nxgrid;
      int *i1, *j1, *i2, *j2;
      double *xarea, *xclon, *xclat, *mask1;

      mask1 = (double *)malloc(nlon1*nlat1*sizeof(double));
      i1    = (int    *)malloc(MAXXGRID*sizeof(int));
      j1    = (int    *)malloc(MAXXGRID*sizeof(int));
      i2    = (int    *)malloc(MAXXGRID*sizeof(int));
      j2    = (int    *)malloc(MAXXGRID*sizeof(int));
      xarea = (double *)malloc(MAXXGRID*sizeof(double));
      xclon = (double *)malloc(MAXXGRID*sizeof(double));
      xclat = (double *)malloc(MAXXGRID*sizeof(double));

      for(i=0; i<nlon1*nlat1; i++) mask1[i] = 1.0;

      nxgrid = create_xgrid_great_circle(&nlon1, &nlat1, &nlon2, &nlat2, lon1_in, lat1_in,
					 lon2_in, lat2_in, mask1, i1, j1, i2, j2,
					 xarea, xclon, xclat);
      printf("     First input grid box longitude, latitude   \n");
      for(i=0; i<n1_in; i++) printf(" %g,  %g \n", lon1_in[i]*R2D, lat1_in[i]*R2D);

      printf("     Second input grid box longitude, latitude \n");
      for(i=0; i<n2_in; i++) printf(" %g,  %g \n", lon2_in[i]*R2D, lat2_in[i]*R2D);

      printf("  Number of exchange grid is %d\n", nxgrid);
      for(i=0; i<nxgrid; i++) {
	printf("(i1,j1)=(%d,%d), (i2,j2)=(%d, %d), xgrid_area=%g, xgrid_clon=%g, xgrid_clat=%g\n",
	       i1[i], j1[i], i2[i], j2[i], xarea[i], xclon[i], xclat[i]);
      }

      /* comparing the area sum of exchange grid and grid1 area */
      {
	double *x1, *y1, *z1, *area1;
	double area_sum;
	int    i;
	area_sum = 0.0;

	for(i=0; i<nxgrid; i++) {
	  area_sum+= xarea[i];
	}

	area1 = (double *)malloc((nlon1)*(nlat1)*sizeof(double));
	get_grid_great_circle_area_(&nlon1, &nlat1, lon1_in, lat1_in, area1);

	printf("xgrid area sum is %g, grid 1 area is %g\n", area_sum, area1[0]);
      }

      printf("\n");
      free(i1);
      free(i2);
      free(j1);
      free(j2);
      free(xarea);
      free(xclon);
      free(xclat);
      free(mask1);
    }
    else if(n>14) {
     // latlon2xyz(n1_in, lon1_in, lat1_in, x1_in, y1_in, z1_in);
     // latlon2xyz(n2_in, lon2_in, lat2_in, x2_in, y2_in, z2_in);

      n_out = clip_2dx2d(lon1_in, lat1_in, n1_in, lon2_in, lat2_in, n2_in, lon_out, lat_out);

      n1_in = fix_lon(lon1_in, lat1_in, n1_in, M_PI);
      n2_in = fix_lon(lon2_in, lat2_in, n2_in, M_PI);
      n_out = fix_lon(lon_out, lat_out, n_out, M_PI);

      double area1 = poly_area (lon1_in, lat1_in, n1_in );
      double area2 = poly_area (lon2_in, lat2_in, n2_in );
      double area_out = poly_area (lon_out, lat_out, n_out );

      printf("     First input grid box longitude, latitude, area= %g \n",area1);
      for(i=0; i<n1_in; i++) printf(" %g,", lon1_in[i]*R2D);
      printf("\n");
      for(i=0; i<n1_in; i++) printf(" %g,", lat1_in[i]*R2D);
      printf("\n");

      printf("     Second input grid box longitude, latitude,area= %g \n",area2);
      for(i=0; i<n2_in; i++) printf(" %g,", lon2_in[i]*R2D);
      printf("\n");
      for(i=0; i<n2_in; i++) printf(" %g,", lat2_in[i]*R2D);
      printf("\n");


      printf("     output clip grid box longitude, latitude, area= %g \n ",area_out);
      printf("n_out= %d \n",n_out);
      for(i=0; i<n_out; i++) printf(" %g,", lon_out[i]*R2D);
      printf("\n");
      for(i=0; i<n_out; i++) printf(" %g,", lat_out[i]*R2D);
      printf("\n");
      if(area1>1.0e14 || area2>1.0e14 || area_out>1.0e14) printf("Error in calculating area !\n");
      if(n==16 || n==20) printf("Must result n_out=5!\n");
      if(n==21) printf("Must result n_out=4!\n");
      if(n==15 || n==17) printf("Must result the second box!\n");
      if(n==18 || n==19) printf("Must result n_out=0!\n");
      if(n==22 || n==23) printf("Same box! area22=area23\n");
      if(n==24 || n==25 || n==26) printf("Same box! area24=area25=area26\n");
    }
    else {
      latlon2xyz(n1_in, lon1_in, lat1_in, x1_in, y1_in, z1_in);
      latlon2xyz(n2_in, lon2_in, lat2_in, x2_in, y2_in, z2_in);

      n_out = clip_2dx2d_great_circle(x1_in, y1_in, z1_in, 4, x2_in, y2_in, z2_in, n2_in,
				      x_out, y_out,  z_out);
      xyz2latlon(n_out, x_out, y_out, z_out, lon_out, lat_out);

      printf("\n*********************************************************\n");
      printf("\n     First input grid box longitude, latitude   \n \n");
      for(i=0; i<n1_in; i++) printf(" %g,", lon1_in[i]*R2D);
      printf("\n");
      for(i=0; i<n1_in; i++) printf(" %g,", lat1_in[i]*R2D);
      printf("\n");

      printf("\n     Second input grid box longitude, latitude \n \n");
      for(i=0; i<n2_in; i++) printf(" %g,", lon2_in[i]*R2D);
      printf("\n");
      for(i=0; i<n2_in; i++) printf(" %g,", lat2_in[i]*R2D);
      printf("\n");

      printf("\n     output clip grid box longitude, latitude for case %d \n \n",n);
      printf("n_out= %d \n",n_out);
      for(i=0; i<n_out; i++) printf(" %g,", lon_out[i]*R2D);
      printf("\n");
      for(i=0; i<n_out; i++) printf(" %g,", lat_out[i]*R2D);
      printf("\n");
    }
  }
}
