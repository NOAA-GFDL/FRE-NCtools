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
#include "mosaic_util.h"
#include "create_xgrid_util.h"
#include "constant.h"

#define AREA_RATIO_THRESH (1.e-6)
#define MASK_THRESH       (0.5)
#define EPSLN8            (1.e-8)
#define EPSLN30           (1.0e-30)
#define EPSLN10           (1.0e-10)
double grid_box_radius(const double *x, const double *y, const double *z, int n);
double dist_between_boxes(const double *x1, const double *y1, const double *z1, int n1,
			  const double *x2, const double *y2, const double *z2, int n2);
int inside_edge(double x0, double y0, double x1, double y1, double x, double y);
int line_intersect_2D_3D(double *a1, double *a2, double *q1, double *q2, double *q3,
		         double *intersect, double *u_a, double *u_q, int *inbound);


/*******************************************************************************
  int get_maxxgrid
  return constants MAXXGRID.
*******************************************************************************/
int get_maxxgrid(void)
{
  return MAXXGRID;
}

int get_maxxgrid_(void)
{
  return get_maxxgrid();
}

/*******************************************************************************
void get_grid_area(const int *nlon, const int *nlat, const double *lon, const double *lat, const double *area)
  return the grid area.
*******************************************************************************/
#ifndef __AIX
void get_grid_area_(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area)
{
  get_grid_area(nlon, nlat, lon, lat, area);
}
#endif

void get_grid_area(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area)
{
  int nx, ny, nxp, i, j, n_in;
  double x_in[20], y_in[20];

  nx = *nlon;
  ny = *nlat;
  nxp = nx + 1;

  for(j=0; j<ny; j++) for(i=0; i < nx; i++) {
    x_in[0] = lon[j*nxp+i];
    x_in[1] = lon[j*nxp+i+1];
    x_in[2] = lon[(j+1)*nxp+i+1];
    x_in[3] = lon[(j+1)*nxp+i];
    y_in[0] = lat[j*nxp+i];
    y_in[1] = lat[j*nxp+i+1];
    y_in[2] = lat[(j+1)*nxp+i+1];
    y_in[3] = lat[(j+1)*nxp+i];
    n_in = fix_lon(x_in, y_in, 4, M_PI);
    area[j*nx+i] = poly_area(x_in, y_in, n_in);
  }

};  /* get_grid_area */

#ifndef __AIX
void get_grid_great_circle_area_(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area)
{
  get_grid_great_circle_area(nlon, nlat, lon, lat, area);

}
#endif

void get_grid_great_circle_area(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area)
{
  int nx, ny, nxp, nyp, i, j;
  int n0, n1, n2, n3;
  struct Node *grid=NULL;
  double *x=NULL, *y=NULL, *z=NULL;


  nx = *nlon;
  ny = *nlat;
  nxp = nx + 1;
  nyp = ny + 1;

  x = (double *)malloc(nxp*nyp*sizeof(double));
  y = (double *)malloc(nxp*nyp*sizeof(double));
  z = (double *)malloc(nxp*nyp*sizeof(double));

  latlon2xyz(nxp*nyp, lon, lat, x, y, z);

  for(j=0; j<ny; j++) for(i=0; i < nx; i++) {
    /* clockwise */
    n0 = j*nxp+i;
    n1 = (j+1)*nxp+i;
    n2 = (j+1)*nxp+i+1;
    n3 = j*nxp+i+1;
    rewindList();
    grid = getNext();
    addEnd(grid, x[n0], y[n0], z[n0], 0, 0, 0, -1);
    addEnd(grid, x[n1], y[n1], z[n1], 0, 0, 0, -1);
    addEnd(grid, x[n2], y[n2], z[n2], 0, 0, 0, -1);
    addEnd(grid, x[n3], y[n3], z[n3], 0, 0, 0, -1);
    area[j*nx+i] = gridArea(grid);
  }

  free(x);
  free(y);
  free(z);

};  /* get_grid_great_circle_area */


void get_grid_area_dimensionless(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area)
{
  int nx, ny, nxp, i, j, n_in;
  double x_in[20], y_in[20];

  nx = *nlon;
  ny = *nlat;
  nxp = nx + 1;

  for(j=0; j<ny; j++) for(i=0; i < nx; i++) {
    x_in[0] = lon[j*nxp+i];
    x_in[1] = lon[j*nxp+i+1];
    x_in[2] = lon[(j+1)*nxp+i+1];
    x_in[3] = lon[(j+1)*nxp+i];
    y_in[0] = lat[j*nxp+i];
    y_in[1] = lat[j*nxp+i+1];
    y_in[2] = lat[(j+1)*nxp+i+1];
    y_in[3] = lat[(j+1)*nxp+i];
    n_in = fix_lon(x_in, y_in, 4, M_PI);
    area[j*nx+i] = poly_area_dimensionless(x_in, y_in, n_in);
  }

};  /* get_grid_area */



void get_grid_area_no_adjust(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area)
{
  int nx, ny, nxp, i, j, n_in;
  double x_in[20], y_in[20];

  nx = *nlon;
  ny = *nlat;
  nxp = nx + 1;

  for(j=0; j<ny; j++) for(i=0; i < nx; i++) {
    x_in[0] = lon[j*nxp+i];
    x_in[1] = lon[j*nxp+i+1];
    x_in[2] = lon[(j+1)*nxp+i+1];
    x_in[3] = lon[(j+1)*nxp+i];
    y_in[0] = lat[j*nxp+i];
    y_in[1] = lat[j*nxp+i+1];
    y_in[2] = lat[(j+1)*nxp+i+1];
    y_in[3] = lat[(j+1)*nxp+i];
    n_in = 4;
    area[j*nx+i] = poly_area_no_adjust(x_in, y_in, n_in);
  }

};  /* get_grid_area_no_adjust */

/*******************************************************************************
   Sutherland-Hodgeman algorithm sequentially clips parts outside 4 boundaries
*******************************************************************************/

int clip(const double lon_in[], const double lat_in[], int n_in, double ll_lon, double ll_lat,
   double ur_lon, double ur_lat, double lon_out[], double lat_out[])
{
  double x_tmp[MV], y_tmp[MV], x_last, y_last;
  int i_in, i_out, n_out, inside_last, inside;

  /* clip polygon with LEFT boundary - clip V_IN to V_TMP */
  x_last = lon_in[n_in-1];
  y_last = lat_in[n_in-1];
  inside_last = (x_last >= ll_lon);
  for (i_in=0,i_out=0;i_in<n_in;i_in++) {

    /* if crossing LEFT boundary - output intersection */
    if ((inside=(lon_in[i_in] >= ll_lon))!=inside_last) {
      x_tmp[i_out] = ll_lon;
      y_tmp[i_out++] = y_last + (ll_lon - x_last) * (lat_in[i_in] - y_last) / (lon_in[i_in] - x_last);
    }

    /* if "to" point is right of LEFT boundary, output it */
    if (inside) {
      x_tmp[i_out]   = lon_in[i_in];
      y_tmp[i_out++] = lat_in[i_in];
    }
    x_last = lon_in[i_in];
    y_last = lat_in[i_in];
    inside_last = inside;
  }
  if (!(n_out=i_out)) return(0);

  /* clip polygon with RIGHT boundary - clip V_TMP to V_OUT */
  x_last = x_tmp[n_out-1];
  y_last = y_tmp[n_out-1];
  inside_last = (x_last <= ur_lon);
  for (i_in=0,i_out=0;i_in<n_out;i_in++) {

    /* if crossing RIGHT boundary - output intersection */
    if ((inside=(x_tmp[i_in] <= ur_lon))!=inside_last) {
      lon_out[i_out]   = ur_lon;
      lat_out[i_out++] = y_last + (ur_lon - x_last) * (y_tmp[i_in] - y_last)
                                                 / (x_tmp[i_in] - x_last);
    }

    /* if "to" point is left of RIGHT boundary, output it */
    if (inside) {
      lon_out[i_out]   = x_tmp[i_in];
      lat_out[i_out++] = y_tmp[i_in];
    }

    x_last = x_tmp[i_in];
    y_last = y_tmp[i_in];
    inside_last = inside;
  }
  if (!(n_out=i_out)) return(0);

  /* clip polygon with BOTTOM boundary - clip V_OUT to V_TMP */
  x_last = lon_out[n_out-1];
  y_last = lat_out[n_out-1];
  inside_last = (y_last >= ll_lat);
  for (i_in=0,i_out=0;i_in<n_out;i_in++) {

    /* if crossing BOTTOM boundary - output intersection */
    if ((inside=(lat_out[i_in] >= ll_lat))!=inside_last) {
      y_tmp[i_out]   = ll_lat;
      x_tmp[i_out++] = x_last + (ll_lat - y_last) * (lon_out[i_in] - x_last) / (lat_out[i_in] - y_last);
    }

    /* if "to" point is above BOTTOM boundary, output it */
    if (inside) {
      x_tmp[i_out]   = lon_out[i_in];
      y_tmp[i_out++] = lat_out[i_in];
    }
    x_last = lon_out[i_in];
    y_last = lat_out[i_in];
    inside_last = inside;
  }
  if (!(n_out=i_out)) return(0);

  /* clip polygon with TOP boundary - clip V_TMP to V_OUT */
  x_last = x_tmp[n_out-1];
  y_last = y_tmp[n_out-1];
  inside_last = (y_last <= ur_lat);
  for (i_in=0,i_out=0;i_in<n_out;i_in++) {

    /* if crossing TOP boundary - output intersection */
    if ((inside=(y_tmp[i_in] <= ur_lat))!=inside_last) {
      lat_out[i_out]   = ur_lat;
      lon_out[i_out++] = x_last + (ur_lat - y_last) * (x_tmp[i_in] - x_last) / (y_tmp[i_in] - y_last);
    }

    /* if "to" point is below TOP boundary, output it */
    if (inside) {
      lon_out[i_out]   = x_tmp[i_in];
      lat_out[i_out++] = y_tmp[i_in];
    }
    x_last = x_tmp[i_in];
    y_last = y_tmp[i_in];
    inside_last = inside;
  }
  return(i_out);
}; /* clip */


/*******************************************************************************
   Revise Sutherland-Hodgeman algorithm to find the vertices of the overlapping
   between any two grid boxes. It return the number of vertices for the exchange grid.
*******************************************************************************/
#ifndef __GNUC__
#pragma acc routine seq
#endif
int clip_2dx2d(const double lon1_in[], const double lat1_in[], int n1_in,
   const double lon2_in[], const double lat2_in[], int n2_in,
   double lon_out[], double lat_out[])
{
  double lon_tmp[MV], lat_tmp[MV];
  double x1_0, y1_0, x1_1, y1_1, x2_0, y2_0, x2_1, y2_1;
  double dx1, dy1, dx2, dy2, determ, ds1, ds2;
  int i_out, n_out, inside_last, inside, i1, i2;
  // used by pimod update
  double lon2_tmp[MV], lat2_tmp[MV];
  int gttwopi=0;

  /* clip polygon with each boundary of the polygon */
  /* We treat lon1_in/lat1_in as clip polygon and lon2_in/lat2_in as subject polygon */
  n_out = n1_in;
#ifndef __GNUC__
#pragma acc loop seq
#endif
  for(i1=0; i1<n1_in; i1++) {
    lon_tmp[i1] = lon1_in[i1];
    lat_tmp[i1] = lat1_in[i1];
    if(lon_tmp[i1]>TPI || lon_tmp[i1]<0.0) gttwopi = 1;
  }
#ifndef __GNUC__
#pragma acc loop seq
#endif
  for(i2=0; i2<n2_in; i2++) {
    lon2_tmp[i2] = lon2_in[i2];
    lat2_tmp[i2] = lat2_in[i2];
  }
  //Some grid boxes near North Pole are clipped wrong (issue #42 )
  //The following heuristic fix seems to work. Why?
  if(gttwopi){pimod(lon_tmp,n1_in);pimod(lon2_tmp,n2_in);}

  x2_0 = lon2_tmp[n2_in-1];
  y2_0 = lat2_tmp[n2_in-1];
  for(i2=0; i2<n2_in; i2++) {
    x2_1 = lon2_tmp[i2];
    y2_1 = lat2_tmp[i2];
    x1_0 = lon_tmp[n_out-1];
    y1_0 = lat_tmp[n_out-1];
    inside_last = inside_edge( x2_0, y2_0, x2_1, y2_1, x1_0, y1_0);
#ifndef __GNUC__
#pragma acc loop seq
#endif
    for(i1=0, i_out=0; i1<n_out; i1++) {
      x1_1 = lon_tmp[i1];
      y1_1 = lat_tmp[i1];
      if((inside = inside_edge(x2_0, y2_0, x2_1, y2_1, x1_1, y1_1)) != inside_last ) {
        /* there is intersection, the line between <x1_0,y1_0> and  <x1_1,y1_1>
           should not parallel to the line between <x2_0,y2_0> and  <x2_1,y2_1>
           may need to consider truncation error */
  dy1 = y1_1-y1_0;
  dy2 = y2_1-y2_0;
  dx1 = x1_1-x1_0;
  dx2 = x2_1-x2_0;
  ds1 = y1_0*x1_1 - y1_1*x1_0;
  ds2 = y2_0*x2_1 - y2_1*x2_0;
  determ = dy2*dx1 - dy1*dx2;
        if(fabs(determ) < EPSLN30) {
          // TODO error handling needs to be openacc friendly
          //error_handler("the line between <x1_0,y1_0> and  <x1_1,y1_1> should not parallel to "
          //              "the line between <x2_0,y2_0> and  <x2_1,y2_1>");
  }
  lon_out[i_out]   = (dx2*ds1 - dx1*ds2)/determ;
  lat_out[i_out++] = (dy2*ds1 - dy1*ds2)/determ;


      }
      if(inside) {
  lon_out[i_out]   = x1_1;
  lat_out[i_out++] = y1_1;
      }
      x1_0 = x1_1;
      y1_0 = y1_1;
      inside_last = inside;
    }
    if(!(n_out=i_out)) return 0;
#ifndef __GNUC__
#pragma acc loop seq
#endif
    for(i1=0; i1<n_out; i1++) {
      lon_tmp[i1] = lon_out[i1];
      lat_tmp[i1] = lat_out[i1];
    }
    /* shift the starting point */
    x2_0 = x2_1;
    y2_0 = y2_1;
  }
  return(n_out);
}; /* clip */

void pimod(double x[],int nn)
{
  for (int i=0;i<nn;i++) {
    if      (x[i] < -M_PI) x[i] += TPI;
    else if (x[i] >  M_PI) x[i] -= TPI;
  }
}

/*******************************************************************************
   Revise Sutherland-Hodgeman algorithm to find the vertices of the overlapping
   between any two grid boxes. It return the number of vertices for the exchange grid.
   Each edge of grid box is a part of great circle. All the points are cartesian
   coordinates. Here we are assuming each polygon is convex.
   RANGE_CHECK_CRITERIA is used to determine if the two grid boxes are possible to be
   overlap. The size should be between 0 and 0.5. The larger the range_check_criteria,
   the more expensive of the computatioin. When the value is close to 0,
   some small exchange grid might be lost. Suggest to use value 0.05 for C48.
*******************************************************************************/

int clip_2dx2d_great_circle(const double x1_in[], const double y1_in[], const double z1_in[], int n1_in,
          const double x2_in[], const double y2_in[], const double z2_in [], int n2_in,
          double x_out[], double y_out[], double z_out[])
{
  struct Node *grid1List=NULL;
  struct Node *grid2List=NULL;
  struct Node *intersectList=NULL;
  struct Node *polyList=NULL;
  struct Node *curList=NULL;
  struct Node *firstIntersect=NULL, *curIntersect=NULL;
  struct Node *temp1=NULL, *temp2=NULL, *temp=NULL;

  int    i1, i2, i1p, i2p, i2p2, npts1, npts2;
  int    nintersect, n_out;
  int    maxiter1, maxiter2, iter1, iter2;
  int    found1, found2, curListNum;
  int    has_inbound, inbound;
  double pt1[MV][3], pt2[MV][3];
  double *p1_0=NULL, *p1_1=NULL;
  double *p2_0=NULL, *p2_1=NULL, *p2_2=NULL;
  double intersect[3];
  double u1, u2;
  double min_x1, max_x1, min_y1, max_y1, min_z1, max_z1;
  double min_x2, max_x2, min_y2, max_y2, min_z2, max_z2;


  /* first check the min and max of (x1_in, y1_in, z1_in) with (x2_in, y2_in, z2_in) */
  min_x1 = minval_double(n1_in, x1_in);
  max_x2 = maxval_double(n2_in, x2_in);
  if(min_x1 >= max_x2+RANGE_CHECK_CRITERIA) return 0;
  max_x1 = maxval_double(n1_in, x1_in);
  min_x2 = minval_double(n2_in, x2_in);
  if(min_x2 >= max_x1+RANGE_CHECK_CRITERIA) return 0;

  min_y1 = minval_double(n1_in, y1_in);
  max_y2 = maxval_double(n2_in, y2_in);
  if(min_y1 >= max_y2+RANGE_CHECK_CRITERIA) return 0;
  max_y1 = maxval_double(n1_in, y1_in);
  min_y2 = minval_double(n2_in, y2_in);
  if(min_y2 >= max_y1+RANGE_CHECK_CRITERIA) return 0;

  min_z1 = minval_double(n1_in, z1_in);
  max_z2 = maxval_double(n2_in, z2_in);
  if(min_z1 >= max_z2+RANGE_CHECK_CRITERIA) return 0;
  max_z1 = maxval_double(n1_in, z1_in);
  min_z2 = minval_double(n2_in, z2_in);
  if(min_z2 >= max_z1+RANGE_CHECK_CRITERIA) return 0;

  rewindList();

  grid1List = getNext();
  grid2List = getNext();
  intersectList = getNext();
  polyList = getNext();

  /* insert points into SubjList and ClipList */
  for(i1=0; i1<n1_in; i1++) addEnd(grid1List, x1_in[i1], y1_in[i1], z1_in[i1], 0, 0, 0, -1);
  for(i2=0; i2<n2_in; i2++) addEnd(grid2List, x2_in[i2], y2_in[i2], z2_in[i2], 0, 0, 0, -1);
  npts1 = length(grid1List);
  npts2 = length(grid2List);

  n_out = 0;
  /* set the inside value */
#ifdef debug_test_create_xgrid
  printf("\nNOTE from clip_2dx2d_great_circle: begin to set inside value grid1List\n");
#endif
  /* first check number of points in grid1 is inside grid2 */

  temp = grid1List;
  while(temp) {
    if(insidePolygon(temp, grid2List))
      temp->isInside = 1;
    else
      temp->isInside = 0;
    temp = getNextNode(temp);
  }

#ifdef debug_test_create_xgrid
  printf("\nNOTE from clip_2dx2d_great_circle: begin to set inside value of grid2List\n");
#endif
  /* check if grid2List is inside grid1List */
  temp = grid2List;

  while(temp) {
    if(insidePolygon(temp, grid1List))
      temp->isInside = 1;
    else
      temp->isInside = 0;
    temp = getNextNode(temp);
  }

  /* make sure the grid box is clockwise */

  /*make sure each polygon is convex, which is equivalent that the great_circle_area is positive */
  if( gridArea(grid1List) <= 0 )
    error_handler("create_xgrid.c(clip_2dx2d_great_circle): grid box 1 is not convex");
  if( gridArea(grid2List) <= 0 )
    error_handler("create_xgrid.c(clip_2dx2d_great_circle): grid box 2 is not convex");

#ifdef debug_test_create_xgrid
  printNode(grid1List, "grid1List");
  printNode(grid2List, "grid2List");
#endif

  /* get the coordinates from grid1List and grid2List.
     Please not npts1 might not equal n1_in, npts2 might not equal n2_in because of pole
  */

  temp = grid1List;
  for(i1=0; i1<npts1; i1++) {
    getCoordinates(temp, pt1[i1]);
    temp = temp->Next;
  }
  temp = grid2List;
  for(i2=0; i2<npts2; i2++) {
    getCoordinates(temp, pt2[i2]);
    temp = temp->Next;
  }

  firstIntersect=getNext();
  curIntersect = getNext();

#ifdef debug_test_create_xgrid
  printf("\n\n************************ Start line_intersect_2D_3D ******************************\n");
#endif
  /* first find all the intersection points */
  nintersect = 0;
  for(i1=0; i1<npts1; i1++) {
    i1p = (i1+1)%npts1;
    p1_0 = pt1[i1];
    p1_1 = pt1[i1p];
    for(i2=0; i2<npts2; i2++) {
      i2p = (i2+1)%npts2;
      i2p2 = (i2+2)%npts2;
      p2_0 = pt2[i2];
      p2_1 = pt2[i2p];
      p2_2 = pt2[i2p2];
#ifdef debug_test_create_xgrid
      printf("\n******************************************************************************\n");
      printf(" i1 = %d, i2 = %d \n", i1, i2);
      printf("********************************************************************************\n");
#endif
      if( line_intersect_2D_3D(p1_0, p1_1, p2_0, p2_1, p2_2, intersect, &u1, &u2, &inbound) ) {

  /* from the value of u1, u2 and inbound, we can partially decide if a point is inside or outside of polygon */

  /* add the intersection into intersetList, The intersection might already be in
     intersectList and will be taken care addIntersect
  */
  if(addIntersect(intersectList, intersect[0], intersect[1], intersect[2], 1, u1, u2, inbound, i1, i1p, i2, i2p)) {
    /* add the intersection into the grid1List */

    if(u1 == 1) {
      insertIntersect(grid1List, intersect[0], intersect[1], intersect[2], 0.0, u2, inbound, p1_1[0], p1_1[1], p1_1[2]);
    }
    else
      insertIntersect(grid1List, intersect[0], intersect[1], intersect[2], u1, u2, inbound, p1_0[0], p1_0[1], p1_0[2]);
    /* when u1 == 0 or 1, need to adjust the vertice to intersect value for roundoff error */
    if(u1==1) {
      p1_1[0] = intersect[0];
      p1_1[1] = intersect[1];
      p1_1[2] = intersect[2];
    }
    else if(u1 == 0) {
      p1_0[0] = intersect[0];
      p1_0[1] = intersect[1];
      p1_0[2] = intersect[2];
    }
    /* add the intersection into the grid2List */
    if(u2==1)
      insertIntersect(grid2List, intersect[0], intersect[1], intersect[2], 0.0, u1, 0, p2_1[0], p2_1[1], p2_1[2]);
    else
      insertIntersect(grid2List, intersect[0], intersect[1], intersect[2], u2, u1, 0, p2_0[0], p2_0[1], p2_0[2]);
    /* when u2 == 0 or 1, need to adjust the vertice to intersect value for roundoff error */
    if(u2==1) {
      p2_1[0] = intersect[0];
      p2_1[1] = intersect[1];
      p2_1[2] = intersect[2];
    }
    else if(u2 == 0) {
      p2_0[0] = intersect[0];
      p2_0[1] = intersect[1];
      p2_0[2] = intersect[2];
    }
  }
      }
    }
  }

  /* set inbound value for the points in intersectList that has inbound == 0,
     this will also set some inbound value of the points in grid1List
  */

  /* get the first point in intersectList has inbound = 2, if not, set inbound value */
  has_inbound = 0;
  /* loop through intersectList to see if there is any has inbound=1 or 2 */
  temp = intersectList;
  nintersect = length(intersectList);
  if(nintersect > 1) {
    getFirstInbound(intersectList, firstIntersect);
    if(firstIntersect->initialized) {
      has_inbound = 1;
    }
  }

  /* when has_inbound == 0, get the grid1List and grid2List */
  if( !has_inbound && nintersect > 1) {
    setInbound(intersectList, grid1List);
    getFirstInbound(intersectList, firstIntersect);
    if(firstIntersect->initialized) has_inbound = 1;
  }

  /* if has_inbound = 1, find the overlapping */
  n_out = 0;

  if(has_inbound) {
    maxiter1 = nintersect;
#ifdef debug_test_create_xgrid
    printf("\nNOTE from clip_2dx2d_great_circle: number of intersect is %d\n", nintersect);
    printf("\n size of grid2List is %d, size of grid1List is %d\n", length(grid2List), length(grid1List));
    printNode(intersectList, "beginning intersection list");
    printNode(grid2List, "beginning clip list");
    printNode(grid1List, "beginning subj list");
    printf("\n************************ End line_intersect_2D_3D **********************************\n\n");
#endif
    temp1 = getNode(grid1List, *firstIntersect);
    if( temp1 == NULL) {
      double lon[10], lat[10];
      int i;
      xyz2latlon(n1_in, x1_in, y1_in, z1_in, lon, lat);
      for(i=0; i< n1_in; i++) printf("lon1 = %g, lat1 = %g\n", lon[i]*R2D, lat[i]*R2D);
      printf("\n");
      xyz2latlon(n2_in, x2_in, y2_in, z2_in, lon, lat);
      for(i=0; i< n2_in; i++) printf("lon2 = %g, lat2 = %g\n", lon[i]*R2D, lat[i]*R2D);
      printf("\n");

      error_handler("firstIntersect is not in the grid1List");
    }
    addNode(polyList, *firstIntersect);
    nintersect--;
#ifdef debug_test_create_xgrid
    printNode(polyList, "polyList at stage 1");
#endif

    /* Loop over the grid1List and grid2List to find again the firstIntersect */
    curList = grid1List;
    curListNum = 0;

    /* Loop through curList to find the next intersection, the loop will end
       when come back to firstIntersect
    */
    copyNode(curIntersect, *firstIntersect);
    iter1 = 0;
    found1 = 0;

    while( iter1 < maxiter1 ) {
#ifdef debug_test_create_xgrid
      printf("\n----------- At iteration = %d\n\n", iter1+1 );
      printNode(curIntersect, "curIntersect at the begining of iter1");
#endif
      /* find the curIntersect in curList and get the next intersection points */
      temp1 =  getNode(curList, *curIntersect);
      temp2 = temp1->Next;
      if( temp2 == NULL ) temp2 = curList;

      maxiter2 = length(curList);
      found2 = 0;
      iter2  = 0;
      /* Loop until find the next intersection */
      while( iter2 < maxiter2 ) {
  int temp2IsIntersect;

  temp2IsIntersect = 0;
  if( isIntersect( *temp2 ) ) { /* copy the point and switch to the grid2List */
    struct Node *temp3;

    /* first check if temp2 is the firstIntersect */
    if( sameNode( *temp2, *firstIntersect) ) {
      found1 = 1;
      break;
    }

    temp3 = temp2->Next;
    if( temp3 == NULL) temp3 = curList;
    if( temp3 == NULL) error_handler("creat_xgrid.c: temp3 can not be NULL");
    found2 = 1;
    /* if next node is inside or an intersection,
       need to keep on curList
    */
    temp2IsIntersect = 1;
    if( isIntersect(*temp3) || (temp3->isInside == 1)  ) found2 = 0;
  }
  if(found2) {
    copyNode(curIntersect, *temp2);
    break;
  }
  else {
    addNode(polyList, *temp2);
#ifdef debug_test_create_xgrid
    printNode(polyList, "polyList at stage 2");
#endif
    if(temp2IsIntersect) {
      nintersect--;
    }
  }
  temp2 = temp2->Next;
  if( temp2 == NULL ) temp2 = curList;
  iter2 ++;
      }
      if(found1) break;

      if( !found2 ) error_handler(" not found the next intersection ");

      /* if find the first intersection, the poly found */
      if( sameNode( *curIntersect, *firstIntersect) ) {
  found1 = 1;
  break;
      }

      /* add curIntersect to polyList and remove it from intersectList and curList */
      addNode(polyList, *curIntersect);
#ifdef debug_test_create_xgrid
      printNode(polyList, "polyList at stage 3");
#endif
      nintersect--;


      /* switch curList */
      if( curListNum == 0) {
  curList = grid2List;
  curListNum = 1;
      }
      else {
  curList = grid1List;
  curListNum = 0;
      }
      iter1++;
    }
    if(!found1) error_handler("not return back to the first intersection");

    /* currently we are only clipping convex polygon to convex polygon */
    if( nintersect > 0) error_handler("After clipping, nintersect should be 0");

    /* copy the polygon to x_out, y_out, z_out */
    temp1 = polyList;
    while (temp1 != NULL) {
      getCoordinate(*temp1, x_out+n_out, y_out+n_out, z_out+n_out);
      temp1 = temp1->Next;
      n_out++;
    }

    /* if(n_out < 3) error_handler(" The clipped region has < 3 vertices"); */
    if( n_out < 3) n_out = 0;
#ifdef debug_test_create_xgrid
    printNode(polyList, "polyList after clipping");
#endif
  }

  /* check if grid1 is inside grid2 */
  if(n_out==0){
    /* first check number of points in grid1 is inside grid2 */
    int n, n1in2;
    /* One possible is that grid1List is inside grid2List */
#ifdef debug_test_create_xgrid
    printf("\nNOTE from clip_2dx2d_great_circle: check if grid1 is inside grid2\n");
#endif
    n1in2 = 0;
    temp = grid1List;
    while(temp) {
      if(temp->intersect != 1) {
#ifdef debug_test_create_xgrid
  printf("grid1->isInside = %d\n", temp->isInside);
#endif
  if( temp->isInside == 1) n1in2++;
      }
      temp = getNextNode(temp);
    }
    if(npts1==n1in2) { /* grid1 is inside grid2 */
      n_out = npts1;
      n = 0;
      temp = grid1List;
      while( temp ) {
  getCoordinate(*temp, &x_out[n], &y_out[n], &z_out[n]);
  n++;
  temp = getNextNode(temp);
      }
    }
    if(n_out>0) return n_out;
  }

  /* check if grid2List is inside grid1List */
  if(n_out ==0){
    int n, n2in1;
#ifdef debug_test_create_xgrid
    printf("\nNOTE from clip_2dx2d_great_circle: check if grid2 is inside grid1\n");
#endif

    temp = grid2List;
    n2in1 = 0;
    while(temp) {
      if(temp->intersect != 1) {
#ifdef debug_test_create_xgrid
  printf("grid2->isInside = %d\n", temp->isInside);
#endif
  if( temp->isInside == 1) n2in1++;
      }
      temp = getNextNode(temp);
    }

    if(npts2==n2in1) { /* grid2 is inside grid1 */
      n_out = npts2;
      n = 0;
      temp = grid2List;
      while( temp ) {
  getCoordinate(*temp, &x_out[n], &y_out[n], &z_out[n]);
  n++;
  temp = getNextNode(temp);
      }

    }
  }


  return n_out;
}


/* Intersects between the line a and the seqment s
   where both line and segment are great circle lines on the sphere represented by
   3D cartesian points.
   [sin sout] are the ends of a line segment
   returns true if the lines could be intersected, false otherwise.
   inbound means the direction of (a1,a2) go inside or outside of (q1,q2,q3)
*/

int line_intersect_2D_3D(double *a1, double *a2, double *q1, double *q2, double *q3,
       double *intersect, double *u_a, double *u_q, int *inbound){

  /* Do this intersection by reprsenting the line a1 to a2 as a plane through the
     two line points and the origin of the sphere (0,0,0). This is the
     definition of a great circle arc.
  */
  double plane[9];
  double plane_p[2];
  double u;
  double p1[3], v1[3], v2[3];
  double c1[3], c2[3], c3[3];
  double coincident, sense, norm;
  int    i;
  int is_inter1, is_inter2;

  *inbound = 0;

  /* first check if any vertices are the same */
  if(samePoint(a1[0], a1[1], a1[2], q1[0], q1[1], q1[2])) {
    *u_a = 0;
    *u_q = 0;
    intersect[0] = a1[0];
    intersect[1] = a1[1];
    intersect[2] = a1[2];
#ifdef debug_test_create_xgrid
    printf("\nNOTE from line_intersect_2D_3D: u_a = %19.15f, u_q=%19.15f, inbound=%d\n", *u_a, *u_q, *inbound);
#endif
    return 1;
   }
   else if (samePoint(a1[0], a1[1], a1[2], q2[0], q2[1], q2[2])) {
    *u_a = 0;
    *u_q = 1;
    intersect[0] = a1[0];
    intersect[1] = a1[1];
    intersect[2] = a1[2];
#ifdef debug_test_create_xgrid
    printf("\nNOTE from line_intersect_2D_3D: u_a = %19.15f, u_q=%19.15f, inbound=%d\n", *u_a, *u_q, *inbound);
#endif
    return 1;
  }
   else if(samePoint(a2[0], a2[1], a2[2], q1[0], q1[1], q1[2])) {
#ifdef debug_test_create_xgrid
    printf("\nNOTE from line_intersect_2D_3D: u_a = %19.15f, u_q=%19.15f, inbound=%d\n", *u_a, *u_q, *inbound);
#endif
    *u_a = 1;
    *u_q = 0;
    intersect[0] = a2[0];
    intersect[1] = a2[1];
    intersect[2] = a2[2];
    return 1;
   }
   else if (samePoint(a2[0], a2[1], a2[2], q2[0], q2[1], q2[2])) {
#ifdef debug_test_create_xgrid
    printf("\nNOTE from line_intersect_2D_3D: u_a = %19.15f, u_q=%19.15f, inbound=%d\n", *u_a, *u_q, *inbound);
#endif
    *u_a = 1;
    *u_q = 1;
    intersect[0] = a2[0];
    intersect[1] = a2[1];
    intersect[2] = a2[2];
    return 1;
  }


  /* Load points defining plane into variable (these are supposed to be in counterclockwise order) */
  plane[0]=q1[0];
  plane[1]=q1[1];
  plane[2]=q1[2];
  plane[3]=q2[0];
  plane[4]=q2[1];
  plane[5]=q2[2];
  plane[6]=0.0;
  plane[7]=0.0;
  plane[8]=0.0;

  /* Intersect the segment with the plane */
  is_inter1 = intersect_tri_with_line(plane, a1, a2, plane_p, u_a);

  if(!is_inter1)
     return 0;

  if(fabs(*u_a) < EPSLN8) *u_a = 0;
  if(fabs(*u_a-1) < EPSLN8) *u_a = 1;


#ifdef debug_test_create_xgrid
    printf("\nNOTE from line_intersect_2D_3D: u_a = %19.15f\n", *u_a);
#endif


  if( (*u_a < 0) || (*u_a > 1) ) return 0;

  /* Load points defining plane into variable (these are supposed to be in counterclockwise order) */
  plane[0]=a1[0];
  plane[1]=a1[1];
  plane[2]=a1[2];
  plane[3]=a2[0];
  plane[4]=a2[1];
  plane[5]=a2[2];
  plane[6]=0.0;
  plane[7]=0.0;
  plane[8]=0.0;

  /* Intersect the segment with the plane */
  is_inter2 = intersect_tri_with_line(plane, q1, q2, plane_p, u_q);

  if(!is_inter2)
     return 0;

  if(fabs(*u_q) < EPSLN8) *u_q = 0;
  if(fabs(*u_q-1) < EPSLN8) *u_q = 1;
#ifdef debug_test_create_xgrid
    printf("\nNOTE from line_intersect_2D_3D: u_q = %19.15f\n", *u_q);
#endif


  if( (*u_q < 0) || (*u_q > 1) ) return 0;

  u =*u_a;

  /* The two planes are coincidental */
  vect_cross(a1, a2, c1);
  vect_cross(q1, q2, c2);
  vect_cross(c1, c2, c3);
  coincident = metric(c3);

  if(fabs(coincident) < EPSLN30) return 0;

  /* Calculate point of intersection */
  intersect[0]=a1[0] + u*(a2[0]-a1[0]);
  intersect[1]=a1[1] + u*(a2[1]-a1[1]);
  intersect[2]=a1[2] + u*(a2[2]-a1[2]);

  norm = metric( intersect );
  for(i = 0; i < 3; i ++) intersect[i] /= norm;

  /* when u_q =0 or u_q =1, the following could not decide the inbound value */
  if(*u_q != 0 && *u_q != 1){

    p1[0] = a2[0]-a1[0];
    p1[1] = a2[1]-a1[1];
    p1[2] = a2[2]-a1[2];
    v1[0] = q2[0]-q1[0];
    v1[1] = q2[1]-q1[1];
    v1[2] = q2[2]-q1[2];
    v2[0] = q3[0]-q2[0];
    v2[1] = q3[1]-q2[1];
    v2[2] = q3[2]-q2[2];

    vect_cross(v1, v2, c1);
    vect_cross(v1, p1, c2);

    sense = dot(c1, c2);
    *inbound = 1;
    if(sense > 0) *inbound = 2; /* v1 going into v2 in CCW sense */
  }
#ifdef debug_test_create_xgrid
    printf("\nNOTE from line_intersect_2D_3D: inbound=%d\n", *inbound);
#endif

  return 1;
}


/*------------------------------------------------------------------------------
  double poly_ctrlat(const double x[], const double y[], int n)
  This routine is used to calculate the latitude of the centroid
  Reference: First- and Second-Order Conservative Remapping Schemes for Grids in
             Spherical Coordinates, P. Jones, Monthly Weather Review, 1998, vol127, p2204
  The following is an implementation of equation (13) in the above paper:
     \int lat.dA = \int_c [-cos(lat)-lat sin(lat)] dlon
  It assumes the sides of the spherical polygons are line segments with tangent
  (lat2-lat1)/(lon2-lon1) between a pair of vertices in approximating the above
  line integral along the sides of the polygon  \int_c.
   ---------------------------------------------------------------------------*/

double poly_ctrlat(const double x[], const double y[], int n)
{
  double ctrlat = 0.0;
  int    i;
  for (i=0;i<n;i++) {
    int ip = (i+1) % n;
    double dx = (x[ip]-x[i]);
    double dy, avg_y, hdy;
    double lat1, lat2;
    lat1 = y[ip];
    lat2 = y[i];
    dy = lat2 - lat1;
    hdy = dy*0.5;
    avg_y = (lat1+lat2)*0.5;
    if      (dx==0.0) continue;
    if(dx > M_PI)  dx = dx - 2.0*M_PI;
    if(dx <= -M_PI) dx = dx + 2.0*M_PI; // flip sign for dx=-pi to fix huge value see comments in function poly_area


    if ( fabs(hdy)< SMALL_VALUE ) /* cheap area calculation along latitude */
      ctrlat -= dx*(2*cos(avg_y) + lat2*sin(avg_y) - cos(lat1) );
    else
      ctrlat -= dx*( (sin(hdy)/hdy)*(2*cos(avg_y) + lat2*sin(avg_y)) - cos(lat1) );
  }
  if(fabs(ctrlat) > HPI) printf("WARNING poly_ctrlat: Large values for ctrlat: %19.15f\n", ctrlat);
  return (ctrlat*RADIUS*RADIUS);
}; /* poly_ctrlat */
/*An alternate implementation of poly_ctrlat for future developments. Under construction.*/
double poly_ctrlat2(const double x[], const double y[], int n)
{
  double ctrlat = 0.0;
  int    i,ip;
  double dx,dy,avg_y, hdy,lat1,lat2,da,dat,dxs= 0.0;
  int hasPole=0, hasBadxm=0, hasBadxp=0;
  for (i=0;i<n;i++) {
    ip = (i+1) % n;
    dx = (x[ip]-x[i]);
    if(fabs(dx+M_PI) < SMALL_VALUE) hasBadxm=1;
    if(fabs(dx-M_PI) < SMALL_VALUE) hasBadxp=1;
    if(y[i]==-HPI || y[i]==HPI) hasPole=1;
  }

  for (i=0;i<n;i++) {
    ip = (i+1) % n;
    dx = (x[ip]-x[i]);
    lat1 = y[ip];
    lat2 = y[i];
    dy = lat2 - lat1;
    hdy = dy*0.5;
    avg_y = (lat1+lat2)*0.5;
    if(dx > M_PI)  dx = dx - 2.0*M_PI;
    if(dx < -M_PI) dx = dx + 2.0*M_PI;

    if ( fabs(hdy)< SMALL_VALUE ) // limit to avoid div by 0
      dat = 1.0;
    else
      dat = sin(hdy)/hdy;

    da = -dx*( dat *(2*cos(avg_y) + lat2*sin(avg_y)) - cos(lat1) );
    ctrlat += da;
    dxs += dx;
    if(hasBadxm || hasBadxp) printf("poly_ctrlat: %19.15f,%19.15f,%19.15f,%19.15f\n", dx,dxs,da,ctrlat);
  }
  if(fabs(ctrlat)>M_PI){
    printf("Error    : Nonzero gridcell dx sum in poly_ctrlat: %19.15f,%19.15f,%19.15f\n",avg_y, dxs,ctrlat);
    ctrlat = fabs(ctrlat) - M_PI*M_PI;
    printf("Corrected: Nonzero gridcell dx sum in poly_ctrlat: %19.15f,%19.15f,%19.15f\n",avg_y, dxs,ctrlat);
  }
  return (ctrlat*RADIUS*RADIUS);
}; /* poly_ctrlat */

/*------------------------------------------------------------------------------
  double poly_ctrlon(const double x[], const double y[], int n, double clon)
  This routine is used to calculate the lontitude of the centroid.
   ---------------------------------------------------------------------------*/
double poly_ctrlon(const double x[], const double y[], int n, double clon)
{
  double ctrlon = 0.0;
  int    i;

  clon = clon;
  for (i=0;i<n;i++) {
    int ip = (i+1) % n;
    double phi1, phi2, dphi, lat1, lat2, dphi1, dphi2;
    double f1, f2, fac, fint;
    phi1   = x[ip];
    phi2   = x[i];
    lat1 = y[ip];
    lat2 = y[i];
    dphi   = phi1 - phi2;

    if      (dphi==0.0) continue;

    f1 = 0.5*(cos(lat1)*sin(lat1)+lat1);
    f2 = 0.5*(cos(lat2)*sin(lat2)+lat2);

     /* this will make sure longitude of centroid is at
        the same interval as the center of any grid */
    if(dphi > M_PI)  dphi = dphi - 2.0*M_PI;
    if(dphi < -M_PI) dphi = dphi + 2.0*M_PI;
    dphi1 = phi1 - clon;
    if( dphi1 > M_PI) dphi1 -= 2.0*M_PI;
    if( dphi1 <-M_PI) dphi1 += 2.0*M_PI;
    dphi2 = phi2 -clon;
    if( dphi2 > M_PI) dphi2 -= 2.0*M_PI;
    if( dphi2 <-M_PI) dphi2 += 2.0*M_PI;

    if(fabs(dphi2 -dphi1) < M_PI) {
      ctrlon -= dphi * (dphi1*f1+dphi2*f2)/2.0;
    }
    else {
      if(dphi1 > 0.0)
  fac = M_PI;
      else
  fac = -M_PI;
      fint = f1 + (f2-f1)*(fac-dphi1)/fabs(dphi);
      ctrlon -= 0.5*dphi1*(dphi1-fac)*f1 - 0.5*dphi2*(dphi2+fac)*f2
  + 0.5*fac*(dphi1+dphi2)*fint;
  }

  }
  return (ctrlon*RADIUS*RADIUS);
};   /* poly_ctrlon */

/* -----------------------------------------------------------------------------
   double box_ctrlat(double ll_lon, double ll_lat, double ur_lon, double ur_lat)
   This routine is used to calculate the latitude of the centroid.
   ---------------------------------------------------------------------------*/
double box_ctrlat(double ll_lon, double ll_lat, double ur_lon, double ur_lat)
{
  double dphi = ur_lon-ll_lon;
  double ctrlat;

  if(dphi > M_PI)  dphi = dphi - 2.0*M_PI;
  if(dphi < -M_PI) dphi = dphi + 2.0*M_PI;
  ctrlat = dphi*(cos(ur_lat) + ur_lat*sin(ur_lat)-(cos(ll_lat) + ll_lat*sin(ll_lat)));
  return (ctrlat*RADIUS*RADIUS);
}; /* box_ctrlat */

/*------------------------------------------------------------------------------
  double box_ctrlon(double ll_lon, double ll_lat, double ur_lon, double ur_lat, double clon)
  This routine is used to calculate the lontitude of the centroid
   ----------------------------------------------------------------------------*/
double box_ctrlon(double ll_lon, double ll_lat, double ur_lon, double ur_lat, double clon)
{
  double phi1, phi2, dphi, lat1, lat2, dphi1, dphi2;
  double f1, f2, fac, fint;
  double ctrlon  = 0.0;
  int i;
  clon = clon;
  for( i =0; i<2; i++) {
    if(i == 0) {
      phi1 = ur_lon;
      phi2 = ll_lon;
      lat1 = lat2 = ll_lat;
    }
    else {
      phi1 = ll_lon;
      phi2 = ur_lon;
      lat1 = lat2 = ur_lat;
    }
    dphi   = phi1 - phi2;
    f1 = 0.5*(cos(lat1)*sin(lat1)+lat1);
    f2 = 0.5*(cos(lat2)*sin(lat2)+lat2);

    if(dphi > M_PI)  dphi = dphi - 2.0*M_PI;
    if(dphi < -M_PI) dphi = dphi + 2.0*M_PI;
    /* make sure the center is in the same grid box. */
    dphi1 = phi1 - clon;
    if( dphi1 > M_PI) dphi1 -= 2.0*M_PI;
    if( dphi1 <-M_PI) dphi1 += 2.0*M_PI;
    dphi2 = phi2 -clon;
    if( dphi2 > M_PI) dphi2 -= 2.0*M_PI;
    if( dphi2 <-M_PI) dphi2 += 2.0*M_PI;

    if(fabs(dphi2 -dphi1) < M_PI) {
      ctrlon -= dphi * (dphi1*f1+dphi2*f2)/2.0;
    }
    else {
      if(dphi1 > 0.0)
  fac = M_PI;
      else
  fac = -M_PI;
      fint = f1 + (f2-f1)*(fac-dphi1)/fabs(dphi);
      ctrlon -= 0.5*dphi1*(dphi1-fac)*f1 - 0.5*dphi2*(dphi2+fac)*f2
  + 0.5*fac*(dphi1+dphi2)*fint;
    }
  }
  return (ctrlon*RADIUS*RADIUS);
} /* box_ctrlon */

/*******************************************************************************
  double grid_box_radius(double *x, double *y, double *z, int n);
  Find the radius of the grid box, the radius is defined the
  maximum distance between any two vertices
*******************************************************************************/
double grid_box_radius(const double *x, const double *y, const double *z, int n)
{
  double radius;
  int i, j;

  radius = 0;
  for(i=0; i<n-1; i++) {
    for(j=i+1; j<n; j++) {
      radius = max(radius, pow(x[i]-x[j],2.)+pow(y[i]-y[j],2.)+pow(z[i]-z[j],2.));
    }
  }

  radius = sqrt(radius);

  return (radius);

}; /* grid_box_radius */

/*******************************************************************************
  double dist_between_boxes(const double *x1, const double *y1, const double *z1, int n1,
          const double *x2, const double *y2, const double *z2, int n2);
  Find the distance between any two grid boxes. The distance is defined by the maximum
  distance between any vertices of these two box
*******************************************************************************/
double dist_between_boxes(const double *x1, const double *y1, const double *z1, int n1,
        const double *x2, const double *y2, const double *z2, int n2)
{
  double dist;
  int i, j;

  dist = 0.0;
  for(i=0; i<n1; i++) {
    for(j=0; j<n2; j++) {
      dist = max(dist, pow(x1[i]-x2[j],2.)+pow(y1[i]-y2[j],2.)+pow(z1[i]-z2[j],2.));
    }
  }

  dist = sqrt(dist);
  return (dist);

}; /* dist_between_boxes */

/*******************************************************************************
 int inside_edge(double x0, double y0, double x1, double y1, double x, double y)
 determine a point(x,y) is inside or outside a given edge with vertex,
 (x0,y0) and (x1,y1). return 1 if inside and 0 if outside. <y1-y0, -(x1-x0)> is
 the outward edge normal from vertex <x0,y0> to <x1,y1>. <x-x0,y-y0> is the vector
 from <x0,y0> to <x,y>.
 if Inner produce <x-x0,y-y0>*<y1-y0, -(x1-x0)> > 0, outside, otherwise inside.
 inner product value = 0 also treate as inside.
*******************************************************************************/
int inside_edge(double x0, double y0, double x1, double y1, double x, double y)
{
   const double SMALL = 1.e-12;
   double product;

   product = ( x-x0 )*(y1-y0) + (x0-x1)*(y-y0);
   return (product<=SMALL) ? 1:0;

 }; /* inside_edge */


/* The following is a test program to test subroutines in create_xgrid.c */
/* To compile: icc -Dtest_create_xgrid -Ddebug_test_create_xgrid -g -O0 ./create_xgrid.c -o test_create_xgrid -L. -lfrencutils */
#ifdef test_create_xgrid

#include "create_xgrid.h"
#include <math.h>
#include <string.h>

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


#endif
