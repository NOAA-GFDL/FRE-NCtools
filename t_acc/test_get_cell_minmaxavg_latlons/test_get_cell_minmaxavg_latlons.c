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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "create_xgrid_utils_acc.h"
#include "parameters.h"

#define NLON 36 // 36 cells in lon direction (36+1 grid points in the lon direction for each lat point)
#define NLAT 4  // 4 cells in lat direction ( 4+1 grid points in the lat direction for each lon point)

double dlon = 360.0/NLON;
double dlat = 30.0;

double const tolerance = 1.e-7;

typedef struct {
  double lon_min[NLON*NLAT];
  double lon_max[NLON*NLAT];
  double lat_min[NLON*NLAT];
  double lat_max[NLON*NLAT];
  double lon_center[NLON*NLAT];
  double *lon_vertices[NLON*NLAT];
  double *lat_vertices[NLON*NLAT];
} Answers;

void get_answers(const double *lon, const double *lat, Answers *answers);
void check_answers( const Minmaxavg_list *minmaxavg_list, const Answers *answers);
void check_ranswer( const int n, const double *answer, const double *checkme);


int main(){

  double lon[(NLON+1)*(NLON+1)], lat[(NLON+1)*(NLON+1)];
  Answers answers;
  Minmaxavg_list minmaxavg_list;

  // no poles, do not need to worry about fix_lon
  for(int ilat=0 ; ilat<NLAT+1 ; ilat++) {
    for(int ilon=0 ; ilon<NLON+1 ; ilon++) {
      int ipt=ilat*(NLON+1)+ilon;
      lat[ipt] = (-30.0 + dlat*ilat)*D2R;
      lon[ipt] = (0.0 + dlon*ilon)*D2R;
    }
  }

  get_cell_minmaxavg_latlons( NLON, NLAT, lon, lat, &minmaxavg_list);
  get_answers(lon, lat, &answers);
  check_answers( &minmaxavg_list, &answers);

  return 0;
}


void get_answers(const double *lon, const double *lat, Answers *answers ){

  for(int icell=0 ; icell<NLON*NLAT ; icell++){
    answers->lon_vertices[icell] = (double *)calloc(MAX_V, sizeof(double));
    answers->lat_vertices[icell] = (double *)calloc(MAX_V, sizeof(double));
  }

  //get min max avg replaceme lon values of a cell
  for(int ilat=0 ; ilat<NLAT ; ilat++) {
    for(int ilon=0 ; ilon<NLON ; ilon++) {
      int icell = ilat*NLON + ilon;
      answers->lat_min[icell] = (-30.0 + dlat*ilat)* D2R;
      answers->lat_max[icell] = (-30.0 + dlat*(ilat+1))* D2R;
      answers->lon_min[icell] = (0 + dlon*ilon)* D2R;
      answers->lon_max[icell] = (0 + dlon*(ilon+1))* D2R;
      answers->lon_center[icell] = (dlon*ilon+0.5*dlon)* D2R;
    }
  }

  //get vertices for each cell
  for(int ilat=0 ; ilat<NLAT ; ilat++) {
    for( int ilon=0 ; ilon<NLON ; ilon++) {
      int icell=ilat*(NLON) + ilon;
      //   3-----2
      //   |     |
      //   |     |
      //   0-----1
      answers->lat_vertices[icell][0]=answers->lat_min[icell];
      answers->lat_vertices[icell][1]=answers->lat_min[icell];
      answers->lat_vertices[icell][2]=answers->lat_max[icell];
      answers->lat_vertices[icell][3]=answers->lat_max[icell];

      answers->lon_vertices[icell][0]=answers->lon_min[icell];
      answers->lon_vertices[icell][1]=answers->lon_max[icell];
      answers->lon_vertices[icell][2]=answers->lon_max[icell];
      answers->lon_vertices[icell][3]=answers->lon_min[icell];
    }
  }

}

void check_answers( const Minmaxavg_list *minmaxavg_list, const Answers *answers)
{

  printf("Checking for nubmer of cell vertices\n");
  for(int i=0 ; i<NLON*NLAT ; i++) {
    if( minmaxavg_list->n_vertices[i] != 4 ) {
      printf("ERROR element %d:  %d vs %d\n", i, minmaxavg_list->n_vertices[i], 4);
      exit(1);
    }
  }

  printf("Checking for min longitudinal point for each cell\n");
  check_ranswer(NLON*NLAT, answers->lon_min, minmaxavg_list->lon_min);

  printf("Checking for max longitudinal point for each cell\n");
  check_ranswer(NLON*NLAT, answers->lon_max, minmaxavg_list->lon_max);

  printf("Checking for min latitudinal point for each cell\n");
  check_ranswer(NLON*NLAT, answers->lat_min, minmaxavg_list->lat_min);

  printf("Checking for max latitudinal point for each cell\n");
  check_ranswer(NLON*NLAT, answers->lat_max, minmaxavg_list->lat_max);

  printf("Checking for center longitudinal point for each cell\n");
  check_ranswer(NLON*NLAT, answers->lon_center, minmaxavg_list->lon_center);

  printf("Checking for longitudinal vertices for each cell\n");
  for(int icell=0 ; icell<NLAT*NLON ; icell++) {
    check_ranswer(MAX_V, answers->lon_vertices[icell], minmaxavg_list->vertices[icell].lon);
  }

  printf("Checking for latitudinal vertices for each cell\n");
  for(int icell=0 ; icell<NLAT*NLON ; icell++) {
    check_ranswer(MAX_V, answers->lat_vertices[icell], minmaxavg_list->vertices[icell].lat);
  }


}

void check_ranswer( const int n, const double *answers, const double *checkme) {
  for(int i=0 ; i<n ; i++) {
    if( abs(answers[i]-checkme[i])>tolerance ) {
      printf("   ERROR element %d:  %lf vs %lf, difference of %e\n",
             i, answers[i], checkme[i], abs(answers[i]-checkme[i]));
      exit(1);
    }
  }
}
