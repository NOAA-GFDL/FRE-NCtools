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
#include <openacc.h>
#include "general_utils_acc.h"
#include "create_xgrid_utils_acc.h"
#include "globals_acc.h"

/*******************************************************************************
void get_grid_area(const int *nlon, const int *nlat, const double *lon, const double *lat, const double *area)
  return the grid area.
*******************************************************************************/
void get_grid_area_acc(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area)
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
      n_in = fix_lon_acc(x_in, y_in, 4, M_PI);
      area[j*nx+i] = poly_area_acc(x_in, y_in, n_in);
    }

};  /* get_grid_area */

void get_grid_great_circle_area_acc(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area)
{
  int nx, ny, nxp, nyp, i, j, n_in;
  int n0, n1, n2, n3;
  double x_in[20], y_in[20], z_in[20];
  struct Node_acc *grid=NULL;
  double *x=NULL, *y=NULL, *z=NULL;


  nx = *nlon;
  ny = *nlat;
  nxp = nx + 1;
  nyp = ny + 1;

  x = (double *)malloc(nxp*nyp*sizeof(double));
  y = (double *)malloc(nxp*nyp*sizeof(double));
  z = (double *)malloc(nxp*nyp*sizeof(double));

  latlon2xyz_acc(nxp*nyp, lon, lat, x, y, z);

  for(j=0; j<ny; j++) for(i=0; i < nx; i++) {
      /* clockwise */
      n0 = j*nxp+i;
      n1 = (j+1)*nxp+i;
      n2 = (j+1)*nxp+i+1;
      n3 = j*nxp+i+1;
      rewindList_acc();
      grid = getNext_acc();
      addEnd_acc(grid, x[n0], y[n0], z[n0], 0, 0, 0, -1);
      addEnd_acc(grid, x[n1], y[n1], z[n1], 0, 0, 0, -1);
      addEnd_acc(grid, x[n2], y[n2], z[n2], 0, 0, 0, -1);
      addEnd_acc(grid, x[n3], y[n3], z[n3], 0, 0, 0, -1);
      area[j*nx+i] = gridArea_acc(grid);
    }

  free(x);
  free(y);
  free(z);

};  /* get_grid_great_circle_area */

/*---------------------------------------------------------------------------
  double poly_ctrlon(const double x[], const double y[], int n, double clon)
  This routine is used to calculate the lontitude of the centroid.
---------------------------------------------------------------------------*/
void poly_ctrlon_acc(const double *x, const double *y, int n, double clon_in, double *ctrlon)
{

  *ctrlon = 0.0;
  for (int i=0;i<n;i++) {
    int ip = (i+1) % n;
    double phi1 = x[ip];
    double phi2 = x[i];
    double dphi = phi1 - phi2;
    double lat1 = y[ip];
    double lat2 = y[i];
    double dphi1, dphi2;
    double f1, f2, fac, fint;

    if (dphi==0.0) continue;

    f1 = 0.5*(cos(lat1)*sin(lat1)+lat1);
    f2 = 0.5*(cos(lat2)*sin(lat2)+lat2);

    //this will make sure longitude of centroid is at
    //the same interval as the center of any grid
    if(dphi > M_PI)  dphi = dphi - 2.0*M_PI;
    if(dphi < -M_PI) dphi = dphi + 2.0*M_PI;
    dphi1 = phi1 - clon_in;
    if( dphi1 > M_PI) dphi1 -= 2.0*M_PI;
    if( dphi1 <-M_PI) dphi1 += 2.0*M_PI;
    dphi2 = phi2 -clon_in;
    if( dphi2 > M_PI) dphi2 -= 2.0*M_PI;
    if( dphi2 <-M_PI) dphi2 += 2.0*M_PI;

    if(fabs(dphi2 -dphi1) < M_PI) *ctrlon -= dphi * (dphi1*f1+dphi2*f2)/2.0;
    else {
      fac = -M_PI; if(dphi1 > 0.0) fac = M_PI;
      fint = f1 + (f2-f1)*(fac-dphi1)/fabs(dphi);
      *ctrlon -= 0.5*dphi1*(dphi1-fac)*f1 - 0.5*dphi2*(dphi2+fac)*f2
        + 0.5*fac*(dphi1+dphi2)*fint;
    }
  }
  *ctrlon = *ctrlon * RADIUS *RADIUS;
  //return (ctrlon*RADIUS*RADIUS);
};   /* poly_ctrlon */

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
void poly_ctrlat_acc(const double *x, const double *y, int n, double *ctrlat)
{

  *ctrlat = 0.0;
  for (int i=0;i<n;i++) {
    int ip = (i+1) % n;
    double dx = (x[ip]-x[i]);
    double lat1 = y[ip];
    double lat2 = y[i];
    double dy = lat2 - lat1;
    double avg_y = (lat1+lat2)*0.5;
    double hdy = dy*0.5;

    if (dx==0.0) continue;
    if(dx > M_PI)   dx = dx - 2.0*M_PI;
    if(dx <= -M_PI) dx = dx + 2.0*M_PI; // flip sign for dx=-pi to fix huge value see comments in function poly_area

    if ( fabs(hdy)< SMALL_VALUE ) /* cheap area calculation along latitude */
      *ctrlat -= dx*(2*cos(avg_y) + lat2*sin(avg_y) - cos(lat1) );
    else
      *ctrlat -= dx*( (sin(hdy)/hdy)*(2*cos(avg_y) + lat2*sin(avg_y)) - cos(lat1) );
  }
  *ctrlat = *ctrlat * RADIUS*RADIUS;
  //if(fabs(ctrlat) > HPI) printf("WARNING poly_ctrlat: Large values for ctrlat: %19.15f\n", ctrlat);
  //return (ctrlat*RADIUS*RADIUS);

}; /* poly_ctrlat */

/*******************************************************************************
   Revise Sutherland-Hodgeman algorithm to find the vertices of the overlapping
   between any two grid boxes. It return the number of vertices for the exchange grid.
*******************************************************************************/
int clip_2dx2d_acc(const double lon1_in[], const double lat1_in[], int n1_in,
                   const double lon2_in[], const double lat2_in[], int n2_in,
                   double lon_out[], double lat_out[])
{
  double lon_tmp[MV], lat_tmp[MV];
  double lon2_tmp[MV], lat2_tmp[MV];
  double x1_0, y1_0, x1_1, y1_1, x2_0, y2_0, x2_1, y2_1;
  double dx1, dy1, dx2, dy2, determ, ds1, ds2;
  int i_out, n_out, inside_last, inside, i1, i2;
  int gttwopi=0;
  /* clip polygon with each boundary of the polygon */
  /* We treat lon1_in/lat1_in as clip polygon and lon2_in/lat2_in as subject polygon */
  n_out = n1_in;
  for(i1=0; i1<n1_in; i1++) {
    lon_tmp[i1] = lon1_in[i1];
    lat_tmp[i1] = lat1_in[i1];
    if(lon_tmp[i1]>TPI || lon_tmp[i1]<0.0) gttwopi = 1;
  }
  for(i2=0; i2<n2_in; i2++) {
    lon2_tmp[i2] = lon2_in[i2];
    lat2_tmp[i2] = lat2_in[i2];
  }
  //Some grid boxes near North Pole are clipped wrong (issue #42 )
  //The following heuristic fix seems to work. Why?
  if(gttwopi){pimod_acc(lon_tmp,n1_in);pimod_acc(lon2_tmp,n2_in);}

  x2_0 = lon2_tmp[n2_in-1];
  y2_0 = lat2_tmp[n2_in-1];
  for(i2=0; i2<n2_in; i2++) {
    x2_1 = lon2_tmp[i2];
    y2_1 = lat2_tmp[i2];
    x1_0 = lon_tmp[n_out-1];
    y1_0 = lat_tmp[n_out-1];
    inside_last = inside_edge_acc( x2_0, y2_0, x2_1, y2_1, x1_0, y1_0);
    for(i1=0, i_out=0; i1<n_out; i1++) {
      x1_1 = lon_tmp[i1];
      y1_1 = lat_tmp[i1];
      if((inside = inside_edge_acc(x2_0, y2_0, x2_1, y2_1, x1_1, y1_1)) != inside_last ) {
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
          printf("ERROR the line between <x1_0,y1_0> and  <x1_1,y1_1> should not parallel to "
                 "the line between <x2_0,y2_0> and  <x2_1,y2_1>\n");
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
int clip_2dx2d_great_circle_acc(const double x1_in[], const double y1_in[], const double z1_in[], int n1_in,
                                const double x2_in[], const double y2_in[], const double z2_in [], int n2_in,
                                double x_out[], double y_out[], double z_out[])
{
  struct Node_acc *subjList=NULL;
  struct Node_acc *clipList=NULL;
  struct Node_acc *grid1List=NULL;
  struct Node_acc *grid2List=NULL;
  struct Node_acc *intersectList=NULL;
  struct Node_acc *polyList=NULL;
  struct Node_acc *curList=NULL;
  struct Node_acc *firstIntersect=NULL, *curIntersect=NULL;
  struct Node_acc *temp1=NULL, *temp2=NULL, *temp=NULL;

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
  static int first_call=1;


  /* first check the min and max of (x1_in, y1_in, z1_in) with (x2_in, y2_in, z2_in) */
  min_x1 = minval_double_acc(n1_in, x1_in);
  max_x2 = maxval_double_acc(n2_in, x2_in);
  if(min_x1 >= max_x2+RANGE_CHECK_CRITERIA) return 0;
  max_x1 = maxval_double_acc(n1_in, x1_in);
  min_x2 = minval_double_acc(n2_in, x2_in);
  if(min_x2 >= max_x1+RANGE_CHECK_CRITERIA) return 0;

  min_y1 = minval_double_acc(n1_in, y1_in);
  max_y2 = maxval_double_acc(n2_in, y2_in);
  if(min_y1 >= max_y2+RANGE_CHECK_CRITERIA) return 0;
  max_y1 = maxval_double_acc(n1_in, y1_in);
  min_y2 = minval_double_acc(n2_in, y2_in);
  if(min_y2 >= max_y1+RANGE_CHECK_CRITERIA) return 0;

  min_z1 = minval_double_acc(n1_in, z1_in);
  max_z2 = maxval_double_acc(n2_in, z2_in);
  if(min_z1 >= max_z2+RANGE_CHECK_CRITERIA) return 0;
  max_z1 = maxval_double_acc(n1_in, z1_in);
  min_z2 = minval_double_acc(n2_in, z2_in);
  if(min_z2 >= max_z1+RANGE_CHECK_CRITERIA) return 0;

  rewindList_acc();

  grid1List = getNext_acc();
  grid2List = getNext_acc();
  intersectList = getNext_acc();
  polyList = getNext_acc();

  /* insert points into SubjList and ClipList */
  for(i1=0; i1<n1_in; i1++) addEnd_acc(grid1List, x1_in[i1], y1_in[i1], z1_in[i1], 0, 0, 0, -1);
  for(i2=0; i2<n2_in; i2++) addEnd_acc(grid2List, x2_in[i2], y2_in[i2], z2_in[i2], 0, 0, 0, -1);
  npts1 = length_acc(grid1List);
  npts2 = length_acc(grid2List);

  n_out = 0;
  /* set the inside value */

  /* first check number of points in grid1 is inside grid2 */

  temp = grid1List;
  while(temp) {
    if(insidePolygon_acc(temp, grid2List))
      temp->isInside = 1;
    else
      temp->isInside = 0;
    temp = getNextNode_acc(temp);
  }

  /* check if grid2List is inside grid1List */
  temp = grid2List;
  while(temp) {
    if(insidePolygon_acc(temp, grid1List))
      temp->isInside = 1;
    else
      temp->isInside = 0;
    temp = getNextNode_acc(temp);
  }

  /* make sure the grid box is clockwise */

  /*make sure each polygon is convex, which is equivalent that the great_circle_area is positive */
  if( gridArea_acc(grid1List) <= 0 )
    printf("ERROR create_xgrid.c(clip_2dx2d_great_circle): grid box 1 is not convex\n");
  if( gridArea_acc(grid2List) <= 0 )
    printf("ERROR create_xgrid.c(clip_2dx2d_great_circle): grid box 2 is not convex\n");

  /* get the coordinates from grid1List and grid2List.
     Please not npts1 might not equal n1_in, npts2 might not equal n2_in because of pole
  */

  temp = grid1List;
  for(i1=0; i1<npts1; i1++) {
    getCoordinates_acc(temp, pt1[i1]);
    temp = temp->Next_acc;
  }
  temp = grid2List;
  for(i2=0; i2<npts2; i2++) {
    getCoordinates_acc(temp, pt2[i2]);
    temp = temp->Next_acc;
  }

  firstIntersect=getNext_acc();
  curIntersect = getNext_acc();

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
      if( line_intersect_2D_3D_acc(p1_0, p1_1, p2_0, p2_1, p2_2, intersect, &u1, &u2, &inbound) ) {
        int n_prev, n_cur;
        int is_in_subj, is_in_clip;

        /* from the value of u1, u2 and inbound, we can partially decide if a point is inside or outside of polygon */

        /* add the intersection into intersetList, The intersection might already be in
           intersectList and will be taken care addIntersect
        */
        if(addIntersect_acc(intersectList, intersect[0], intersect[1], intersect[2], 1, u1, u2, inbound, i1, i1p, i2, i2p)) {
          /* add the intersection into the grid1List */

          if(u1 == 1) {
            insertIntersect_acc(grid1List, intersect[0], intersect[1], intersect[2], 0.0, u2, inbound, p1_1[0], p1_1[1], p1_1[2]);
          }
          else
            insertIntersect_acc(grid1List, intersect[0], intersect[1], intersect[2], u1, u2, inbound, p1_0[0], p1_0[1], p1_0[2]);
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
            insertIntersect_acc(grid2List, intersect[0], intersect[1], intersect[2], 0.0, u1, 0, p2_1[0], p2_1[1], p2_1[2]);
          else
            insertIntersect_acc(grid2List, intersect[0], intersect[1], intersect[2], u2, u1, 0, p2_0[0], p2_0[1], p2_0[2]);
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
  nintersect = length_acc(intersectList);
  if(nintersect > 1) {
    getFirstInbound_acc(intersectList, firstIntersect);
    if(firstIntersect->initialized) {
      has_inbound = 1;
    }
  }

  /* when has_inbound == 0, get the grid1List and grid2List */
  if( !has_inbound && nintersect > 1) {
    setInbound_acc(intersectList, grid1List);
    getFirstInbound_acc(intersectList, firstIntersect);
    if(firstIntersect->initialized) has_inbound = 1;
  }

  /* if has_inbound = 1, find the overlapping */
  n_out = 0;

  if(has_inbound) {
    maxiter1 = nintersect;
    temp1 = getNode_acc(grid1List, *firstIntersect);
    if( temp1 == NULL) {
      double lon[10], lat[10];
      int i;
      xyz2latlon_acc(n1_in, x1_in, y1_in, z1_in, lon, lat);
      for(i=0; i< n1_in; i++) printf("lon1 = %g, lat1 = %g\n", lon[i]*R2D, lat[i]*R2D);
      printf("\n");
      xyz2latlon_acc(n2_in, x2_in, y2_in, z2_in, lon, lat);
      for(i=0; i< n2_in; i++) printf("lon2 = %g, lat2 = %g\n", lon[i]*R2D, lat[i]*R2D);
      printf("\n");

      printf("ERROR firstIntersect is not in the grid1List\n");
    }
    addNode_acc(polyList, *firstIntersect);
    nintersect--;

    /* Loop over the grid1List and grid2List to find again the firstIntersect */
    curList = grid1List;
    curListNum = 0;

    /* Loop through curList to find the next intersection, the loop will end
       when come back to firstIntersect
    */
    copyNode_acc(curIntersect, *firstIntersect);
    iter1 = 0;
    found1 = 0;

    while( iter1 < maxiter1 ) {
      /* find the curIntersect in curList and get the next intersection points */
      temp1 =  getNode_acc(curList, *curIntersect);
      temp2 = temp1->Next_acc;
      if( temp2 == NULL ) temp2 = curList;

      maxiter2 = length_acc(curList);
      found2 = 0;
      iter2  = 0;
      /* Loop until find the next intersection */
      while( iter2 < maxiter2 ) {
        int temp2IsIntersect;

        temp2IsIntersect = 0;
        if( isIntersect_acc( *temp2 ) ) { /* copy the point and switch to the grid2List */
          struct Node_acc *temp3;

          /* first check if temp2 is the firstIntersect */
          if( sameNode_acc( *temp2, *firstIntersect) ) {
            found1 = 1;
            break;
          }

          temp3 = temp2->Next_acc;
          if( temp3 == NULL) temp3 = curList;
          if( temp3 == NULL) printf("ERROR creat_xgrid.c: temp3 can not be NULL\n");
          found2 = 1;
          /* if next node is inside or an intersection,
             need to keep on curList
          */
          temp2IsIntersect = 1;
          if( isIntersect_acc(*temp3) || (temp3->isInside == 1)  ) found2 = 0;
        }
        if(found2) {
          copyNode_acc(curIntersect, *temp2);
          break;
        }
        else {
          addNode_acc(polyList, *temp2);
          if(temp2IsIntersect) {
            nintersect--;
          }
        }
        temp2 = temp2->Next_acc;
        if( temp2 == NULL ) temp2 = curList;
        iter2 ++;
      } // while( iter2<maxiter2 )
      if(found1) break;

      if( !found2 ) printf("ERROR not found the next intersection\n ");

      /* if find the first intersection, the poly found */
      if( sameNode_acc( *curIntersect, *firstIntersect) ) {
        found1 = 1;
        break;
      }

      /* add curIntersect to polyList and remove it from intersectList and curList */
      addNode_acc(polyList, *curIntersect);
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
    } // while(iter1)
    if(!found1) printf("ERROR not return back to the first intersection\n");

    /* currently we are only clipping convex polygon to convex polygon */
    if( nintersect > 0) printf("ERROR After clipping, nintersect should be 0\n");

    /* copy the polygon to x_out, y_out, z_out */
    temp1 = polyList;
    while (temp1 != NULL) {
      getCoordinate_acc(*temp1, x_out+n_out, y_out+n_out, z_out+n_out);
      temp1 = temp1->Next_acc;
      n_out++;
    }

    /* if(n_out < 3) error_handler(" The clipped region has < 3 vertices"); */
    if( n_out < 3) n_out = 0;
  } //if(inbound)

  /* check if grid1 is inside grid2 */
  if(n_out==0){
    /* first check number of points in grid1 is inside grid2 */
    int n, n1in2;
    /* One possible is that grid1List is inside grid2List */
    n1in2 = 0;
    temp = grid1List;
    while(temp) {
      if(temp->intersect != 1) {
        if( temp->isInside == 1) n1in2++;
      }
      temp = getNextNode_acc(temp);
    }
    if(npts1==n1in2) { /* grid1 is inside grid2 */
      n_out = npts1;
      n = 0;
      temp = grid1List;
      while( temp ) {
        getCoordinate_acc(*temp, &x_out[n], &y_out[n], &z_out[n]);
        n++;
        temp = getNextNode_acc(temp);
      }
    }
    if(n_out>0) return n_out;
  }

  /* check if grid2List is inside grid1List */
  if(n_out ==0){
    int n, n2in1;

    temp = grid2List;
    n2in1 = 0;
    while(temp) {
      if(temp->intersect != 1) {
        if( temp->isInside == 1) n2in1++;
      }
      temp = getNextNode_acc(temp);
    }

    if(npts2==n2in1) { /* grid2 is inside grid1 */
      n_out = npts2;
      n = 0;
      temp = grid2List;
      while( temp ) {
        getCoordinate_acc(*temp, &x_out[n], &y_out[n], &z_out[n]);
        n++;
        temp = getNextNode_acc(temp);
      }

    }
  }

  return n_out;
}

void get_grid_cells_struct_acc( const int nlon, const int nlat, const double *lon, const double *lat,
                                Grid_cells_struct_config *grid_cells )
{

  int ncells=nlon*nlat;
  int npts=(nlon+1)*(nlat+1);

  grid_cells->lon_min   = (double *)malloc(ncells*sizeof(double));
  grid_cells->lon_max   = (double *)malloc(ncells*sizeof(double));
  grid_cells->lat_min   = (double *)malloc(ncells*sizeof(double));
  grid_cells->lat_max   = (double *)malloc(ncells*sizeof(double));
  grid_cells->lon_cent  = (double *)malloc(ncells*sizeof(double));
  grid_cells->area      = (double *)malloc(ncells*sizeof(double));
  grid_cells->nvertices = (int *)malloc(ncells*sizeof(int));
  grid_cells->lon_vertices  = (double **)malloc(ncells*sizeof(double));
  grid_cells->lat_vertices  = (double **)malloc(ncells*sizeof(double));
  for(int icell=0 ; icell<ncells ; icell++) {
    grid_cells->lat_vertices[icell] = (double *)malloc(MAX_V*sizeof(double));
    grid_cells->lon_vertices[icell] = (double *)malloc(MAX_V*sizeof(double));
  }

#pragma acc enter data create(grid_cells[:1])
#pragma acc enter data create(grid_cells->lon_min[:ncells], grid_cells->lon_max[:ncells], \
                              grid_cells->lat_min[:ncells], grid_cells->lat_max[:ncells], \
                              grid_cells->lon_cent[:ncells], grid_cells->nvertices[:ncells],\
                              grid_cells->area[:ncells])
#pragma acc enter data create(grid_cells->lon_vertices[:ncells][:MAX_V], \
                              grid_cells->lat_vertices[:ncells][:MAX_V])


#pragma acc data present(grid_cells[:1], lon[:npts], lat[:npts])
#pragma acc parallel loop seq //independent
  for(int icell=0; icell<ncells; icell++){
    int nvertices;
    double lon_vertices[MV], lat_vertices[MV];

    get_cell_vertices_acc( icell, nlon, lon, lat, lon_vertices, lat_vertices );

    grid_cells->lat_min[icell] = minval_double_acc(4, lat_vertices);
    grid_cells->lat_max[icell] = maxval_double_acc(4, lat_vertices);

    nvertices = fix_lon_acc(lon_vertices, lat_vertices, 4, M_PI);
    grid_cells->nvertices[icell] = nvertices;

    if(nvertices>MAX_V) printf("ERROR get_cell_minmaxavg_latlons:  number of cell vertices is greater than MAX_V\n");
    grid_cells->lon_min[icell]  = minval_double_acc(nvertices, lon_vertices);
    grid_cells->lon_max[icell]  = maxval_double_acc(nvertices, lon_vertices);
    grid_cells->lon_cent[icell] = avgval_double_acc(nvertices, lon_vertices);

    grid_cells->area[icell] = poly_area_acc(lon_vertices, lat_vertices, nvertices);

    for(int ivertex=0 ; ivertex<nvertices ; ivertex++) {
      grid_cells->lon_vertices[icell][ivertex] = lon_vertices[ivertex];
      grid_cells->lat_vertices[icell][ivertex] = lat_vertices[ivertex];
    }
  }

}

void get_cell_vertices_acc( const int icell, const int nlon, const double *lon, const double *lat, double *x, double *y )
{

  int i, j;
  int n0, n1, n2, n3;

  i = icell%nlon;
  j = icell/nlon;
  n0 = j*(nlon+1)+i;
  n1 = j*(nlon+1)+i+1;
  n2 = (j+1)*(nlon+1)+i+1;
  n3 = (j+1)*(nlon+1)+i;

  x[0] = lon[n0]; y[0] = lat[n0];
  x[1] = lon[n1]; y[1] = lat[n1];
  x[2] = lon[n2]; y[2] = lat[n2];
  x[3] = lon[n3]; y[3] = lat[n3];

}

void create_upbound_nxcells_arrays_on_device_acc(const int n, int **approx_nxcells_per_ij1,
                                                 int **ij2_start, int **ij2_end)
{

  int *p_approx_nxcells_per_ij1;
  int *p_ij2_start;
  int *p_ij2_end;

  *approx_nxcells_per_ij1 = (int *)malloc(n*sizeof(int));
  *ij2_start = (int *)malloc(n*sizeof(int));
  *ij2_end = (int *)malloc(n*sizeof(int));

  p_approx_nxcells_per_ij1 = *approx_nxcells_per_ij1;
  p_ij2_start = *ij2_start;
  p_ij2_end = *ij2_end;

#pragma acc enter data copyin(p_approx_nxcells_per_ij1[:n],   \
                              p_ij2_start[:n],                \
                              p_ij2_end[:n])

}

void free_upbound_nxcells_array_from_all_acc( const int n, int *approx_nxcells_per_ij1,
                                              int *ij2_start, int *ij2_end)
{
#pragma acc exit data delete(approx_nxcells_per_ij1[:n],  \
                             ij2_start[:n],               \
                             ij2_end[:n])

  free(approx_nxcells_per_ij1);
  free(ij2_start);
  free(ij2_end);
}

void free_output_grid_cell_struct_from_all_acc(const int n, Grid_cells_struct_config *grid_cells)
{

#pragma acc exit data delete(grid_cells->lon_min[:n],   \
                             grid_cells->lon_max[:n],   \
                             grid_cells->lon_min[:n],   \
                             grid_cells->lat_max[:n],   \
                             grid_cells->lat_min[:n],   \
                             grid_cells->area[:n],      \
                             grid_cells->nvertices[:n])
  for(int icell=0 ; icell<n; icell++){
#pragma acc exit data delete(grid_cells->lon_vertices[icell][:MAX_V], \
                             grid_cells->lat_vertices[icell][:MAX_V])
  }
#pragma acc exit data delete(grid_cells->lon_vertices[:n],  \
                             grid_cells->lat_vertices[:n])
#pragma acc exit data delete(grid_cells[:1])

  free(grid_cells->lon_min);
  free(grid_cells->lon_max);
  free(grid_cells->lat_min);
  free(grid_cells->lat_max);
  free(grid_cells->lon_cent);
  free(grid_cells->area);
  free(grid_cells->nvertices);
  for(int icell=0 ; icell<n ; icell++) {
    free(grid_cells->lat_vertices[icell]);
    free(grid_cells->lon_vertices[icell]);
  }
  free(grid_cells->lon_vertices);
  free(grid_cells->lat_vertices);

}

void copy_data_to_xgrid_on_device_acc(const int nxcells, const int input_ncells, const int upbound_nxcells,
                                      int *xcells_per_ij1, double *xcell_dclon, double *xcell_dclat,
                                      int *approx_xcells_per_ij1, int *parent_input_indices, int *parent_output_indices,
                                      double *xcell_areas, Xinfo_per_input_tile *xgrid_for_input_tile)
{
  int copy_xcentroid = (*xcell_dclon != -99.99 ) ? nxcells : 0;
  xgrid_for_input_tile->nxcells = nxcells;

  if(nxcells>0) {

    xgrid_for_input_tile->input_parent_cell_indices = (int *)malloc(nxcells*sizeof(int));
    xgrid_for_input_tile->output_parent_cell_indices = (int *)malloc(nxcells*sizeof(int));
    xgrid_for_input_tile->xcell_area = (double *)malloc(nxcells*sizeof(double));
    if(copy_xcentroid) {
      xgrid_for_input_tile->dcentroid_lon = (double *)malloc(nxcells*sizeof(double));
      xgrid_for_input_tile->dcentroid_lat = (double *)malloc(nxcells*sizeof(double));
    }

#pragma acc enter data copyin(xgrid_for_input_tile)
#pragma acc enter data create(xgrid_for_input_tile->input_parent_cell_indices[:nxcells], \
                              xgrid_for_input_tile->output_parent_cell_indices[:nxcells], \
                              xgrid_for_input_tile->xcell_area[:nxcells])
#pragma acc enter data create(xgrid_for_input_tile->dcentroid_lon[:nxcells], \
                              xgrid_for_input_tile->dcentroid_lat[:nxcells]) if(copy_xcentroid)
#pragma acc data present(xcells_per_ij1[:input_ncells], approx_xcells_per_ij1[:input_ncells], \
                         parent_input_indices[:upbound_nxcells],        \
                         parent_output_indices[:upbound_nxcells],       \
                         xcell_areas[:upbound_nxcells], xgrid_for_input_tile[:1])
#pragma acc data present(xcell_dclon[:nxcells], xcell_dclat[:nxcells]) if(copy_xcentroid)
#pragma acc parallel loop independent
    for(int ij1=0 ; ij1<input_ncells ; ij1++) {
      int xcells_before_ij1 = 0, approx_xcells=0 ;

#pragma acc loop
      for(int i=1 ; i<=ij1 ; i++) {
        xcells_before_ij1 += xcells_per_ij1[i-1];
        approx_xcells += approx_xcells_per_ij1[i-1];
      }

#pragma acc loop independent
      for(int i=0 ; i<xcells_per_ij1[ij1]; i++){
        xgrid_for_input_tile->input_parent_cell_indices[xcells_before_ij1+i]  = parent_input_indices[approx_xcells+i];
        xgrid_for_input_tile->output_parent_cell_indices[xcells_before_ij1+i] = parent_output_indices[approx_xcells+i];
        xgrid_for_input_tile->xcell_area[xcells_before_ij1+i] = xcell_areas[approx_xcells+i];
      }

      if(copy_xcentroid) {
#pragma acc loop independent
        for(int i=0 ; i<xcells_per_ij1[ij1]; i++){
          xgrid_for_input_tile->dcentroid_lon[xcells_before_ij1+i] = xcell_dclon[approx_xcells+i];
          xgrid_for_input_tile->dcentroid_lat[xcells_before_ij1+i] = xcell_dclat[approx_xcells+i];
        }
      }
    }
  }//if nxcells>0

  }
