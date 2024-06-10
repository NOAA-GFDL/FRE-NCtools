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
#include "interp_utils_acc.h"
#include "parameters.h"
#include "globals_acc.h"

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
  double lon_cent[NLON*NLAT];
  double *lon_vertices[NLON*NLAT];
  double *lat_vertices[NLON*NLAT];
} Answers;

void copy_cell_to_host(Grid_cells_struct_config *grid_cells);
void reset_cell_on_host(Grid_cells_struct_config *grid_cells);
void get_answers(const double *lon, const double *lat, Answers *answers);
void check_answers( const Grid_cells_struct_config *grid_cells, const Answers *answers);
void check_ranswer( const int n, const double *answer, const double *checkme);

//TODO:  test for compute poly_area

int main(){

  Answers answers;
  Grid_cells_struct_config grid_cells;
  Grid_config grid[1];

  grid[0].nxc=NLON;
  grid[0].nyc=NLAT;
  grid[0].lonc = (double *)calloc( (NLON+1)*(NLAT+1), sizeof(double) );
  grid[0].latc = (double *)calloc( (NLON+1)*(NLAT+1), sizeof(double) );

  // no poles, do not need to worry about fix_lon
  for(int ilat=0 ; ilat<NLAT+1 ; ilat++) {
    double latitude = (-30.0 + dlat*ilat)*D2R;
    for(int ilon=0 ; ilon<NLON+1 ; ilon++) {
      int ipt=ilat*(NLON+1)+ilon;
      grid[0].latc[ipt] = latitude;
      grid[0].lonc[ipt] = (0.0 + dlon*ilon)*D2R;
    }
  }

  copy_grid_to_device_acc((NLON+1)*(NLAT+1), grid[0].latc, grid[0].lonc);
  get_grid_cells_struct_acc( NLON, NLAT, grid[0].lonc, grid[0].latc, &grid_cells);

  get_answers(grid[0].lonc, grid[0].latc, &answers);

  reset_cell_on_host(&grid_cells); //to verify data is copied out from device
  copy_cell_to_host(&grid_cells);

  check_answers(&grid_cells, &answers);

  printf("TODO:  add check for areas");

  return 0;
}

void get_answers(const double *lon, const double *lat, Answers *answers ){

  for(int icell=0 ; icell<NLON*NLAT ; icell++){
    answers->lon_vertices[icell] = (double *)calloc(MAX_V, sizeof(double));
    answers->lat_vertices[icell] = (double *)calloc(MAX_V, sizeof(double));
  }

  //get min max avg replaceme lon values of a cell
  for(int ilat=0 ; ilat<NLAT ; ilat++) {
    double latitude_min = (-30.0 + dlat*ilat)* D2R;
    double latitude_max = (-30.0 + dlat*(ilat+1))* D2R;
    for(int ilon=0 ; ilon<NLON ; ilon++) {
      int icell = ilat*NLON + ilon;
      answers->lat_min[icell] = latitude_min;
      answers->lat_max[icell] = latitude_max;
      answers->lon_min[icell] = (0 + dlon*ilon)* D2R;
      answers->lon_max[icell] = (0 + dlon*(ilon+1))* D2R;
      answers->lon_cent[icell] = (dlon*ilon+0.5*dlon)* D2R;
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

void reset_cell_on_host(Grid_cells_struct_config *grid_cells)
{

  for( int icell=0 ; icell<NLON*NLAT; icell++) {
    grid_cells->lon_min[icell]=-99.;
    grid_cells->lon_max[icell]=-99.;
    grid_cells->lat_min[icell]=-99.;
    grid_cells->lat_max[icell]=-99.;
    grid_cells->lon_cent[icell]=-99.;
    grid_cells->nvertices[icell]=-99;
    for(int ipt=0 ; ipt<MAX_V ; ipt++){
      grid_cells->lon_vertices[icell][ipt]=-99.;
      grid_cells->lat_vertices[icell][ipt]=-99.;
    }
  }

}

void copy_cell_to_host(Grid_cells_struct_config *grid_cells)
{

  int ncell = NLON*NLAT;
#pragma acc exit data copyout( grid_cells->lon_min[:ncell], grid_cells->lon_max[:ncell], \
                               grid_cells->lat_min[:ncell], grid_cells->lat_max[:ncell], \
                               grid_cells->lon_cent[:ncell], grid_cells->nvertices[:ncell], \
                               grid_cells->lon_vertices[:ncell][:MAX_V],\
                               grid_cells->lat_vertices[:ncell][:MAX_V])
}

void check_answers( const Grid_cells_struct_config *grid_cells, const Answers *answers)
{

  printf("Checking for nubmer of cell vertices\n");
  for(int i=0 ; i<NLON*NLAT ; i++) {
    if( grid_cells->nvertices[i] != 4 ) {
      printf("ERROR element %d:  %d vs %d\n", i, grid_cells->nvertices[i], 4);
      exit(1);
    }
  }

  printf("Checking for min longitudinal point for each cell\n");
  check_ranswer(NLON*NLAT, answers->lon_min, grid_cells->lon_min);

  printf("Checking for max longitudinal point for each cell\n");
  check_ranswer(NLON*NLAT, answers->lon_max, grid_cells->lon_max);

  printf("Checking for min latitudinal point for each cell\n");
  check_ranswer(NLON*NLAT, answers->lat_min, grid_cells->lat_min);

  printf("Checking for max latitudinal point for each cell\n");
  check_ranswer(NLON*NLAT, answers->lat_max, grid_cells->lat_max);

  printf("Checking for avg longitudinal point for each cell\n");
  check_ranswer(NLON*NLAT, answers->lon_cent, grid_cells->lon_cent);

  printf("Checking for longitudinal vertices for each cell\n");
  for(int icell=0 ; icell<NLAT*NLON ; icell++) {
    check_ranswer(MAX_V, answers->lon_vertices[icell], grid_cells->lon_vertices[icell]);
  }

  printf("Checking for latitudinal vertices for each cell\n");
  for(int icell=0 ; icell<NLAT*NLON ; icell++) {
    check_ranswer(MAX_V, answers->lat_vertices[icell], grid_cells->lat_vertices[icell]);
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