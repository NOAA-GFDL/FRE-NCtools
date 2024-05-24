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

/**
 * \author Zhi Liang
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "general_utils_acc.h"
#include "globals_acc.h"

const double from_pole_threshold_rad_acc = 0.0174533;  // 1.0 deg

//reproduce siena was removed
int rotate_poly_flag_acc = 0;
double the_rotation_matrix_acc[3][3] = { 0 };

#pragma acc routine seq
void set_rotate_poly_true_acc(void){
  rotate_poly_flag_acc = 1;
  set_the_rotation_matrix_acc();
}

struct Node_acc *nodeList_acc=NULL;
int curListPos_acc=0;

#pragma acc declare copyin(from_pole_threshold_rad_acc)
#pragma acc declare copyin(rotate_poly_flag_acc)
#pragma acc declare copyin(nodeList_acc)
#pragma acc declare copyin(curListPos_acc)
#pragma acc declare copyin(the_rotation_matrix_acc)

/*********************************************************************

   int nearest_index(double value, const double *array, int ia)

   return index of nearest data point within "array" corresponding to "value".
   if "value" is outside the domain of "array" then nearest_index = 0
   or = size(array)-1 depending on whether array(0) or array(ia-1) is
   closest to "value"

   Arguments:
     value:  arbitrary data...same units as elements in "array"
     array:  array of data points  (must be monotonically increasing)
     ia   :  size of array.
 ********************************************************************/
int nearest_index_acc(double value, const double *array, int ia)
{
  int index, i;
  int keep_going;

  for(i=1; i<ia; i++){
    if (array[i] < array[i-1])
      printf("ERROR:  nearest_index: array must be monotonically increasing\n");
  }
  if (value < array[0] ) index = 0;
  else if ( value > array[ia-1]) index = ia-1;
  else {
    i=0;
    keep_going = 1;
    while (i < ia && keep_going) {
      i = i+1;
      if (value <= array[i]) {
        index = i;
        if (array[i]-value > value-array[i-1]) index = i-1;
        keep_going = 0;
      }
    }
  }
  return index;

}

/*******************************************************************************
  double maxval_double(int size, double *data)
  get the maximum value of double array
*******************************************************************************/
double maxval_double_acc(int size, const double *data)
{
  int n;
  double maxval;

  maxval = data[0];
  for(n=1; n<size; n++){
    if( data[n] > maxval ) maxval = data[n];
  }

  return maxval;

} /* maxval_double */


/*******************************************************************************
  double minval_double(int size, double *data)
  get the minimum value of double array
*******************************************************************************/
double minval_double_acc(int size, const double *data)
{
  int n;
  double minval;

  minval = data[0];
  for(n=1; n<size; n++){
    if( data[n] < minval ) minval = data[n];
  }

  return minval;

} /* minval_double */

/*******************************************************************************
  double avgval_double(int size, double *data)
  get the average value of double array
*******************************************************************************/
double avgval_double_acc(int size, const double *data)
{
  int n;
  double avgval;

  avgval = 0;
  for(n=0; n<size; n++) avgval += data[n];
  avgval /= size;

  return avgval;

} /* avgval_double */


/*******************************************************************************
  void latlon2xyz
  Routine to map (lon, lat) to (x,y,z)
******************************************************************************/
void latlon2xyz_acc(int size, const double *lon, const double *lat, double *x, double *y, double *z)
{
  int n;

  for(n=0; n<size; n++) {
    x[n] = cos(lat[n])*cos(lon[n]);
    y[n] = cos(lat[n])*sin(lon[n]);
    z[n] = sin(lat[n]);
  }

} /* latlon2xyz */

/*------------------------------------------------------------
       void xyz2laton(np, p, xs, ys)
   Transfer cartesian coordinates to spherical coordinates
   ----------------------------------------------------------*/
void xyz2latlon_acc( int np, const double *x, const double *y, const double *z, double *lon, double *lat)
{

  double xx, yy, zz;
  double dist;
  int i;

  for(i=0; i<np; i++) {
    xx = x[i];
    yy = y[i];
    zz = z[i];
    dist = sqrt(xx*xx+yy*yy+zz*zz);
    xx /= dist;
    yy /= dist;
    zz /= dist;

    if ( fabs(xx)+fabs(yy)  < EPSLN10 ) lon[i] = 0;
    else lon[i] = atan2(yy, xx);
    lat[i] = asin(zz);

    if ( lon[i] < 0.) lon[i] = 2.*M_PI + lon[i];
  }

} /* xyz2latlon */

/*------------------------------------------------------------------------------
  double poly_area(const x[], const y[], int n)
  obtains area of input polygon by line integrating -sin(lat)d(lon)
  Vertex coordinates must be in Radians.
  Vertices must be listed counter-clockwise around polygon.
  grid is in radians.

  Reference: First- and Second-Order Conservative Remapping Schemes for Grids in
             Spherical Coordinates, P. Jones, Monthly Weather Review, 1998, vol127, p2204
  The following is an implementation of equation (12) in the above paper:
     \int dA = \int_c [-sin(lat)] dlon

  An alternative derivation of line integrating -sin(lat)d(lon) formula:
  Consider a vector function in spherical coordinates (r,lon,lat)
  with only a lon component :
      A=(0, (1-sin(lat))/cos(lat)/r , 0)
  Then
      Curl(A)=(1/r**2 , 0, 0) .
  Now consider any loop C on the suface of the sphere enclosing an area S
  and apply the Stokes theorem:
      \integral_surface_S Curl(A).da = \integral_loop_C A.dl
  where da and dl are the vectorial suface and line elements on sphere:
      da=(da, 0, 0) and dl=(0, R*cos(lat)*d(lon), R*d(lat)).
  We get
      \integral_surface_S Curl(A) = \integral_surface_S (da/r**2) = S/R**2
  and
      \integral_loop_C A.dl = \integral_loop_C (1-sin(lat))*d(lon)

  Hence per Stokes formula:
      S/R**2 = \integral_loop_C d(lon) - \integral_loop_C sin(lat)*d(lon).
             = I1 - I2

  Now the approximation used for the second loop integral
  I2= \integral_loop_C sin(lat)*d(lon)
    = sum_over_loop_segments \integral_loop_C sin(lat)*d(lon)
    = sum_over_loop_segments \integral_loop_C sin(lat)*d(lat) *d(lon)/d(lat)

  If d(lon)/d(lat) is assumed constant over a loop segment, given that the segments
  used here are the sides of a grid cell, then we can take it outside the integral
    sum_over_loop_segments \integral_loop_C sin(lat)*d(lat) *d(lon)/d(lat)
   =sum_over_loop_segments delta(lon)/delta(lat) \int_segment sin(lat)*d(lat)
   =sum_over_grid_cell_sides (lon2-lon1)/(lat2-lat1) * (-(cos(lat2)-cos(lat1)))
   =sum_over_grid_cell_sides (lon2-lon1)/(lat2-lat1) * (2*sin(0.5*(lat2+lat1))*sin(0.5*(lat2-lat1)))

  So, finally:
   S/R**2 = I1 - I2
          = I1 - sum_over_grid_cell_sides (lon2-lon1)* sin(0.5*(lat2+lat1))* sin(0.5*(lat2-lat1))/(0.5*(lat2-lat1))
  We can prove that the first integral I1 above is zero for grid cells that do not include or cross a pole.

  Special cases:
   In general I1=0, but for the grid cells that have the N Pole as a vertex:
      I1=- angle at pole vertex
   Particularly,
      I1=-pi/2 for regular CS grids (or close to -pi/2 for stregtched CS grids)
   and
      I1=-pi for grid cells with a side passing through the N pole
      I1=+pi for grid cells with a side passing through the S pole
   You could easily see that for a spherical grid cell bounded by the North pole, two longitude and a lattitude to the south.
   There are 8 such cells in a CS grids and a few more in the exchange grids.

   Alternatively, we can treat such a pole vertex as the limit of an additional (imagined) grid segment that shrinks
   as its two end points approach the pole. Such segment contributes pi/2 to the line integral I2.
   So, as an implementaion trick, if for each pole vertex we insert an additional twin pole vertex with a longitude the same
   as the next vertex in the cell (which is 90 degrees or close to 90 degrees appart) the contribution to I1 shifts into -I2
   and we can again assume I1=0. This scheme is implemented in subroutine fix_lon() which is always called before poly_area().
   That is why in the implementation below I1 is totally absent.

   In some pathological grid cells that do not have a pole vertex, I1 may not come out zero due to ambiguities in choosing
   d(lon) mod 2pi and we have to be careful in the implementation to deal with those grid cells.
   As we saw the approximation for I2 is based on the assumtion that  d(lon)/d(lat) is
   almost a constant over each segment of the path. This approximation is also equivalent
   to assuming a linear variation of latitude versus the longitude along the loop segments (grid cell sides)
   i.e., lat = lat1 + (lon-lon1)*(lat2-lat1)/(lon2-lon1)
   This assumption breaks down if the segment passes through a pole, e.g. along a longitude great circle where
   lon abruptly changes by pi as it goes through the pole like in the following two grid cells:
       x----x----x
       |    |    |
       |   Pole  |
       |    |    |
       |    |    |
       x----x----x
    x denotes the grid points.
    We can diagnose such situations
      either via I1=sum(dlons) coming out as 2*pi instead of 0, in that case we can correct the area by 2*pi
      or via dx=-pi, in that case we can just flip the sign of dx and continue
    These two turn out to be equivalent and so in the implementation below we add the equal sign inside the if conditional   if(dx <= -M_PI) dx = dx + 2.0*M_PI;
    The reason for the sign flip is delicate. As we saw above the contribution to I1 from a grid cell side
    passing through the S pole is +pi regardless of the direction.

  Side note:
  The above approximation for I2 is not very accurate, patricularly for CS grid cells near the Poles.
  This will manifest itself in the form of a discontinuity or bump in the calculated area of CS grid cells as the
  Poles are approached (see fre-nctools issue #44). This is because the sides of these grid cells are
  great circles but the above approximation ignores that fact and the line integral is calculated
  along an arbitrary line joining the vertices. This situation can indeed be cured by integrating along
  the great circle joining the vertices. However the implementation becomes rather tedious. For future
  reference here is the paramtrized equation of the curves that we should integrate along, that is
  the great circle passing at (lon0,lat0) closest to the N pole (excluding the pole itself):
     tan(lat)=-tan(lat0)*cos(lon-lon0).
  So, I2=\integral dx sin(arctan(tan(lat0)*cos(x-lon0)))  can be shown to give an accurate estimate of grid
  cell areas surrounding the pole without a bump.
   ----------------------------------------------------------------------------*/
double poly_area_main_acc(const double x[], const double y[], int n) {
  double area = 0.0;
  int i;

  for (i = 0; i < n; i++) {
    int ip = (i + 1) % n;
    double dx = (x[ip] - x[i]);
    double lat1, lat2;
    double dy, dat;

    lat1 = y[ip];
    lat2 = y[i];
    if (dx > M_PI)
      dx = dx - 2.0 * M_PI;
    if (dx < -M_PI)
      dx = dx + 2.0 * M_PI;
    /*Sides that go through a pole contribute PI regardless of their direction and extent.*/
    if (fabs(dx + M_PI) < SMALL_VALUE || fabs(dx - M_PI) < SMALL_VALUE) {
      area += M_PI;
      continue;  // next i
    }
    /*
     We have to be careful in implementing sin(0.5*(lat2-lat1))/(0.5*(lat2-lat1))
     in the limit of lat1=lat2 where it becomes 1. Note that in that limit we get the well known formula
     for the area of a section of sphere bounded by two longitudes and two latitude circles.
    */
    if (fabs(lat1 - lat2) < SMALL_VALUE) /* cheap area calculation along latitude */
      area -= dx * sin(0.5 * (lat1 + lat2));
    else {
      dy = 0.5 * (lat1 - lat2);
      dat = sin(dy) / dy;
      area -= dx * sin(0.5 * (lat1 + lat2)) * dat;
    }
  }


  if (fabs(area) > HPI) printf("WARNING poly_area: Large values for area: %19.15f\n", area);
  if (area < 0) return -area * RADIUS * RADIUS;
  else return area * RADIUS * RADIUS;
} /* poly_area_main */

/*
  Calculate the area of a polygon with the original poly_area function,
  exept that for those polygons that are polar, if the global poly_are_flag
  is set, the area is calculated by first rotating a copy of the polygon away
  from the pole and calculating the area of the rotated polygon instead.
  Polar means having any of its verices cloase to a pole, or an edge crossing
  a pole.

  This routine is here merely for improving the accuracy of the are calculation
  in the given conditions.
  TODO: The tiling error reported by make_coupler mosaic may be non-zero when
  using this feature. This may be developped in the future.
*/
double poly_area_acc(const double xo[], const double yo[], int n) {
  double area_pa = 0.0;
  double area_par = 0.0;
  int pole = 0;
  int crosses = 0;

  double xr[8];  // rotated lon
  double yr[8];  // rotated lat

  if (rotate_poly_flag_acc == 0) {
    area_pa = poly_area_main_acc(xo, yo, n);
    return area_pa;
  }
  else {
    // anything near enough to the pole gets rotated
    pole = is_near_pole_acc(yo, n);
    crosses = crosses_pole_acc(xo, n);
    if (crosses == 1 && pole == 0) {
      printf("ERROR poly_area:  crosses == 1 && pole == 0\n");
    }

    if (pole == 1) {
      if (n > 8) {
        printf("ERROR poly_area: n > 8. n=%d,n\n");
      }
      rotate_poly_acc(xo, yo, n, xr, yr);

      int pole2 = is_near_pole_acc(yr, n);
      if (pole2 == 1) {
        printf("ERROR poly_area: pole2 == 1\n");
      }
      area_par = poly_area_main_acc(xr, yr, n);
    } else {
      area_pa = poly_area_main_acc(xo, yo, n);
    }

    if (pole == 1) {
      return area_par;
    } else {
      return area_pa;
    }
  }
}

int delete_vtx_acc(double x[], double y[], int n, int n_del)
{
  for (;n_del<n-1;n_del++) {
    x[n_del] = x[n_del+1];
    y[n_del] = y[n_del+1];
  }

  return (n-1);
} /* delete_vtx */

int insert_vtx_acc(double x[], double y[], int n, int n_ins, double lon_in, double lat_in)
{
  int i;

  for (i=n-1;i>=n_ins;i--) {
    x[i+1] = x[i];
    y[i+1] = y[i];
  }

  x[n_ins] = lon_in;
  y[n_ins] = lat_in;
  return (n+1);
} /* insert_vtx */

void v_print_acc(double x[], double y[], int n)
{
  int i;

  for (i=0;i<n;i++) printf(" %20g   %20g\n", x[i]*R2D, y[i]*R2D);
} /* v_print */

int fix_lon_acc(double x[], double y[], int n, double tlon)
{
  double x_sum, dx;
  int i, nn = n, pole = 0;

  for (i=0;i<nn;i++) if (fabs(y[i])>=HPI-TOLERANCE) pole = 1;
  if (0&&pole) {
    printf("fixing pole cell\n");
    v_print_acc(x, y, nn);
    printf("---------");
  }

  /* all pole points must be paired */
  /* The reason is poly_area() function needs a contribution equal to the angle (in radians)
     between the sides that connect to the pole. */
  for (i=0;i<nn;i++) if (fabs(y[i])>=HPI-TOLERANCE) {
      int im=(i+nn-1)%nn, ip=(i+1)%nn;

      if (y[im]==y[i] && y[ip]==y[i]) {
        nn = delete_vtx_acc(x, y, nn, i);
        i--;
      }
      else if (y[im]!=y[i] && y[ip]!=y[i]) {
        nn = insert_vtx_acc(x, y, nn, i, x[i], y[i]);
        i++;
      }
  }
  /* first of pole pair has longitude of previous vertex */
  /* second of pole pair has longitude of subsequent vertex */
  for (i=0;i<nn;i++) if (fabs(y[i])>=HPI-TOLERANCE) {
    int im=(i+nn-1)%nn, ip=(i+1)%nn;

    if (y[im]!=y[i]) x[i] = x[im];
    if (y[ip]!=y[i]) x[i] = x[ip];
  }

  /*If a polygon side passes through a Pole insert twin vertices at the Pole*/
  /*A fix is also directly applied to poly_area to handle this case.*/
  for (i=0;i<nn;i++) {
    int im=(i+nn-1)%nn;
    // ip=(i+1)%nn;
    double dx = x[i]-x[im];
    if(fabs(dx+M_PI)< SMALL_VALUE || fabs(dx-M_PI)< SMALL_VALUE){
      double x1=x[im];
      double x2=x[i];
      double ypole= HPI;
      if(y[i]<0.0) ypole = -HPI ;
      nn = insert_vtx_acc(x, y, nn, i, x2, ypole);
      nn = insert_vtx_acc(x, y, nn, i, x1, ypole);
      break;
    }
  }
  if (nn) x_sum = x[0]; else return(0);
  for (i=1;i<nn;i++) {
    double dx = x[i]-x[i-1];

    if      (dx < -M_PI) dx = dx + TPI;
    else if (dx >  M_PI) dx = dx - TPI;
    x_sum += (x[i] = x[i-1] + dx);
  }

  dx = (x_sum/nn)-tlon;
  if      (dx < -M_PI) for (i=0;i<nn;i++) x[i] += TPI;
  else if (dx >  M_PI) for (i=0;i<nn;i++) x[i] -= TPI;

  if (0&&pole) {
    printf("area=%g\n", poly_area_acc(x, y,nn));
    v_print_acc(x, y, nn);
    printf("---------");
  }

  return (nn);
} /* fix_lon */


/* Compute the great circle area of a polygon on a sphere */
double great_circle_area_acc(int n, const double *x, const double *y, const double *z) {
  int i;
  double pnt0[3], pnt1[3], pnt2[3];
  double sum, area;

  /* sum angles around polygon */
  sum=0.0;
  for ( i=0; i<n; i++) {
    /* points that make up a side of polygon */
    pnt0[0] = x[i];
    pnt0[1] = y[i];
    pnt0[2] = z[i];
    pnt1[0] = x[(i+1)%n];
    pnt1[1] = y[(i+1)%n];
    pnt1[2] = z[(i+1)%n];
    pnt2[0] = x[(i+2)%n];
    pnt2[1] = y[(i+2)%n];
    pnt2[2] = z[(i+2)%n];

    /* compute angle for pnt1 */
    sum += spherical_angle_acc(pnt1, pnt2, pnt0);

  }
  area = (sum - (n-2.)*M_PI) * RADIUS* RADIUS;
  return area;
}

/*------------------------------------------------------------------------------
  double spherical_angle(const double *p1, const double *p2, const double *p3)
           p3
         /
        /
       p1 ---> angle
         \
          \
           p2
 -----------------------------------------------------------------------------*/
double spherical_angle_acc(const double *v1, const double *v2, const double *v3)
{
  double angle;
  double px, py, pz, qx, qy, qz, ddd;

  /* vector product between v1 and v2 */
  px = v1[1]*v2[2] - v1[2]*v2[1];
  py = v1[2]*v2[0] - v1[0]*v2[2];
  pz = v1[0]*v2[1] - v1[1]*v2[0];
  /* vector product between v1 and v3 */
  qx = v1[1]*v3[2] - v1[2]*v3[1];
  qy = v1[2]*v3[0] - v1[0]*v3[2];
  qz = v1[0]*v3[1] - v1[1]*v3[0];

  ddd = (px*px+py*py+pz*pz)*(qx*qx+qy*qy+qz*qz);
  if ( ddd <= 0.0 ) angle = 0. ;
  else {
    ddd = (px*qx+py*qy+pz*qz) / sqrt(ddd);
    if( fabs(ddd-1) < EPSLN30 ) ddd = 1;
    if( fabs(ddd+1) < EPSLN30 ) ddd = -1;
    if ( ddd>1. || ddd<-1. ) {
      /*FIX (lmh) to correctly handle co-linear points (angle near pi or 0) */
      if (ddd < 0.) angle = M_PI;
      else angle = 0.;
    }
    else
      angle = acos( ddd );
  }

  return angle;
} /* spherical_angle */

/*----------------------------------------------------------------------
    void vect_cross(e, p1, p2)
    Perform cross products of 3D vectors: e = P1 X P2
    -------------------------------------------------------------------*/
void vect_cross_acc(const double *p1, const double *p2, double *e )
{

  e[0] = p1[1]*p2[2] - p1[2]*p2[1];
  e[1] = p1[2]*p2[0] - p1[0]*p2[2];
  e[2] = p1[0]*p2[1] - p1[1]*p2[0];

} /* vect_cross */


/*----------------------------------------------------------------------
    double* vect_cross(p1, p2)
    return cross products of 3D vectors: = P1 X P2
    -------------------------------------------------------------------*/
double dot_acc(const double *p1, const double *p2)
{

  return( p1[0]*p2[0] + p1[1]*p2[1] + p1[2]*p2[2] );

}

double metric_acc(const double *p) {
  return (sqrt(p[0]*p[0] + p[1]*p[1]+p[2]*p[2]) );
}

/* Intersect a line and a plane
   Intersects between the plane ( three points ) (entries in counterclockwise order)
   and the line determined by the endpoints l1 and l2 (t=0.0 at l1 and t=1.0 at l2)
   returns true if the two intersect and the output variables are valid
   outputs p containing the coordinates in the tri and t the coordinate in the line
   of the intersection.
   NOTE: the intersection doesn't have to be inside the tri or line for this to return true
*/
int intersect_tri_with_line_acc(const double *plane, const double *l1, const double *l2, double *p,
                                double *t)
{

  double M[3*3], inv_M[3*3];
  double V[3];
  double X[3];
  int is_invert=0;

  const double *pnt0=plane;
  const double *pnt1=plane+3;
  const double *pnt2=plane+6;

  /* To do intersection just solve the set of linear equations for both
     Setup M
  */

  M[0]=l1[0]-l2[0]; M[1]=pnt1[0]-pnt0[0]; M[2]=pnt2[0]-pnt0[0];
  M[3]=l1[1]-l2[1]; M[4]=pnt1[1]-pnt0[1]; M[5]=pnt2[1]-pnt0[1];
  M[6]=l1[2]-l2[2]; M[7]=pnt1[2]-pnt0[2]; M[8]=pnt2[2]-pnt0[2];
  /* Invert M */
  is_invert = invert_matrix_3x3_acc(M,inv_M);
  if (!is_invert) return 0;

  /* Set variable holding vector */
  V[0]=l1[0]-pnt0[0];
  V[1]=l1[1]-pnt0[1];
  V[2]=l1[2]-pnt0[2];

  /* Calculate solution */
  mult_acc(inv_M, V, X);

  /* Get answer out */
  *t=X[0];
  p[0]=X[1];
  p[1]=X[2];

  return 1;
}

void mult_acc(double m[], double v[], double out_v[])
{

  out_v[0]=m[0]*v[0]+m[1]*v[1]+m[2]*v[2];
  out_v[1]=m[3]*v[0]+m[4]*v[1]+m[5]*v[2];
  out_v[2]=m[6]*v[0]+m[7]*v[1]+m[8]*v[2];

}

/* returns 1 if matrix is inverted, 0 otherwise */
int invert_matrix_3x3_acc(double m[], double m_inv[])
{

  const double det =  m[0] * (m[4]*m[8] - m[5]*m[7])
                     -m[1] * (m[3]*m[8] - m[5]*m[6])
                     +m[2] * (m[3]*m[7] - m[4]*m[6]);

  if (fabs(det) < EPSLN15 ) return 0;

  const double deti = 1.0/det;

  m_inv[0] = (m[4]*m[8] - m[5]*m[7]) * deti;
  m_inv[1] = (m[2]*m[7] - m[1]*m[8]) * deti;
  m_inv[2] = (m[1]*m[5] - m[2]*m[4]) * deti;

  m_inv[3] = (m[5]*m[6] - m[3]*m[8]) * deti;
  m_inv[4] = (m[0]*m[8] - m[2]*m[6]) * deti;
  m_inv[5] = (m[2]*m[3] - m[0]*m[5]) * deti;

  m_inv[6] = (m[3]*m[7] - m[4]*m[6]) * deti;
  m_inv[7] = (m[1]*m[6] - m[0]*m[7]) * deti;
  m_inv[8] = (m[0]*m[4] - m[1]*m[3]) * deti;

  return 1;
}



void rewindList_acc(void)
{
  int n;

  curListPos_acc = 0;
  if(!nodeList_acc) nodeList_acc = (struct Node_acc *)malloc(MAXNODELIST*sizeof(struct Node_acc));
  for(n=0; n<MAXNODELIST; n++) initNode_acc(nodeList_acc+n);

}

struct Node_acc *getNext_acc()
{
  struct Node_acc *temp=NULL;
  int n;

  if(!nodeList_acc) {
    nodeList_acc = (struct Node_acc *)malloc(MAXNODELIST*sizeof(struct Node_acc));
    for(n=0; n<MAXNODELIST; n++) initNode_acc(nodeList_acc+n);
  }

  temp = nodeList_acc+curListPos_acc;
  curListPos_acc++;
  if(curListPos_acc > MAXNODELIST) printf("ERROR getNext_acc: curListPos_acc >= MAXNODELIST\n");

  return (temp);
}


void initNode_acc(struct Node_acc *node)
{
    node->x = 0;
    node->y = 0;
    node->z = 0;
    node->u = 0;
    node->intersect = 0;
    node->inbound = 0;
    node->isInside = 0;
    node->Next_acc = NULL;
    node->initialized=0;

}

void addEnd_acc(struct Node_acc *list, double x, double y, double z, int intersect, double u, int inbound, int inside)
{

  struct Node_acc *temp=NULL;

  if(list == NULL) printf("ERROR addEnd: list is NULL\n");

  if(list->initialized) {

    /* (x,y,z) might already in the list when intersect is true and u=0 or 1 */
      temp = list;
      while (temp) {
        if(samePoint_acc(temp->x, temp->y, temp->z, x, y, z)) return;
        temp=temp->Next_acc;
      }
    temp = list;
    while(temp->Next_acc)
      temp=temp->Next_acc;

    /* Append at the end of the list.  */
    temp->Next_acc = getNext_acc();
    temp = temp->Next_acc;
  }
  else {
    temp = list;
  }

  temp->x = x;
  temp->y = y;
  temp->z = z;
  temp->u = u;
  temp->intersect = intersect;
  temp->inbound = inbound;
  temp->initialized=1;
  temp->isInside = inside;
}

/* return 1 if the point (x,y,z) is added in the list, return 0 if it is already in the list */

int addIntersect_acc(struct Node_acc *list, double x, double y, double z, int intersect, double u1, double u2, int inbound,
       int is1, int ie1, int is2, int ie2)
{

  double u1_cur, u2_cur;
  int    i1_cur, i2_cur;
  struct Node_acc *temp=NULL;

  if(list == NULL) printf("ERROR addEnd: list is NULL\n");

  /* first check to make sure this point is not in the list */
  u1_cur = u1;
  i1_cur = is1;
  u2_cur = u2;
  i2_cur = is2;
  if(u1_cur == 1) {
    u1_cur = 0;
    i1_cur = ie1;
  }
  if(u2_cur == 1) {
    u2_cur = 0;
    i2_cur = ie2;
  }

  if(list->initialized) {
    temp = list;
    while(temp) {
      if( temp->u == u1_cur && temp->subj_index == i1_cur) return 0;
      if( temp->u_clip == u2_cur && temp->clip_index == i2_cur) return 0;
      if( !temp->Next_acc ) break;
      temp=temp->Next_acc;
    }

    /* Append at the end of the list.  */
    temp->Next_acc = getNext_acc();
    temp = temp->Next_acc;
  }
  else {
    temp = list;
  }

  temp->x = x;
  temp->y = y;
  temp->z = z;
  temp->intersect = intersect;
  temp->inbound = inbound;
  temp->initialized=1;
  temp->isInside = 0;
  temp->u = u1_cur;
  temp->subj_index = i1_cur;
  temp->u_clip = u2_cur;
  temp->clip_index = i2_cur;

  return 1;
}


int length_acc(struct Node_acc *list)
{
   struct Node_acc *cur_ptr=NULL;
   int count=0;

   cur_ptr=list;

   while(cur_ptr)
   {
     if(cur_ptr->initialized ==0) break;
      cur_ptr=cur_ptr->Next_acc;
      count++;
   }
   return(count);
}

/* two points are the same if there are close enough */
int samePoint_acc(double x1, double y1, double z1, double x2, double y2, double z2)
{
    if( fabs(x1-x2) > EPSLN10 || fabs(y1-y2) > EPSLN10 || fabs(z1-z2) > EPSLN10 )
      return 0;
    else
      return 1;
}



int sameNode_acc(struct Node_acc node1, struct Node_acc node2)
{
  if( node1.x == node2.x && node1.y == node2.y && node1.z==node2.z )
    return 1;
  else
    return 0;
}


void addNode_acc(struct Node_acc *list, struct Node_acc inNode_acc)
{

  addEnd_acc(list, inNode_acc.x, inNode_acc.y, inNode_acc.z, inNode_acc.intersect, inNode_acc.u, inNode_acc.inbound, inNode_acc.isInside);

}

struct Node_acc *getNode_acc(struct Node_acc *list, struct Node_acc inNode_acc)
{
  struct Node_acc *thisNode_acc=NULL;
  struct Node_acc *temp=NULL;

  temp = list;
  while( temp ) {
    if( sameNode_acc( *temp, inNode_acc ) ) {
      thisNode_acc = temp;
      temp = NULL;
      break;
    }
    temp = temp->Next_acc;
  }

  return thisNode_acc;
}

struct Node_acc *getNextNode_acc(struct Node_acc *list)
{
  return list->Next_acc;
}

void copyNode_acc(struct Node_acc *node_out, struct Node_acc node_in)
{

  node_out->x = node_in.x;
  node_out->y = node_in.y;
  node_out->z = node_in.z;
  node_out->u = node_in.u;
  node_out->intersect = node_in.intersect;
  node_out->inbound   = node_in.inbound;
  node_out->Next_acc = NULL;
  node_out->initialized = node_in.initialized;
  node_out->isInside = node_in.isInside;
}

void printNode_acc(struct Node_acc *list, char *str)
{
  struct Node_acc *temp;

  if(list == NULL) printf("ERROR printNode_acc: list is NULL\n");
  if(str) printf("  %s \n", str);
  temp = list;
  while(temp) {
    if(temp->initialized ==0) break;
    printf(" (x, y, z, interset, inbound, isInside) = (%19.15f,%19.15f,%19.15f,%d,%d,%d)\n",
      temp->x, temp->y, temp->z, temp->intersect, temp->inbound, temp->isInside);
    temp = temp->Next_acc;
  }
  printf("\n");
}

int intersectInList_acc(struct Node_acc *list, double x, double y, double z)
{
  struct Node_acc *temp;
  int found=0;

  temp = list;
  found = 0;
  while ( temp ) {
    if( temp->x == x && temp->y == y && temp->z == z ) {
      found = 1;
      break;
    }
    temp=temp->Next_acc;
  }
  if (!found) printf("ERROR intersectInList: point (x,y,z) is not found in the list\n");
  if( temp->intersect == 2 )
    return 1;
  else
    return 0;

}


/* The following insert a intersection after non-intersect point (x2,y2,z2), if the point
   after (x2,y2,z2) is an intersection, if u is greater than the u value of the intersection,
   insert after, otherwise insert before
*/
void insertIntersect_acc(struct Node_acc *list, double x, double y, double z, double u1, double u2, int inbound,
                     double x2, double y2, double z2)
{
  struct Node_acc *temp1=NULL, *temp2=NULL;
  struct Node_acc *temp;
  double u_cur;
  int found=0;

  temp1 = list;
  found = 0;
  while ( temp1 ) {
    if( temp1->x == x2 && temp1->y == y2 && temp1->z == z2 ) {
      found = 1;
      break;
    }
    temp1=temp1->Next_acc;
  }
  if (!found) printf("ERROR inserAfter: point (x,y,z) is not found in the list\n");

  /* when u = 0 or u = 1, set the grid point to be the intersection point to solve truncation error isuse */
  u_cur = u1;
  if(u1 == 1) {
    u_cur = 0;
    temp1 = temp1->Next_acc;
    if(!temp1) temp1 = list;
  }
  if(u_cur==0) {
    temp1->intersect = 2;
    temp1->isInside = 1;
    temp1->u = u_cur;
    temp1->x = x;
    temp1->y = y;
    temp1->z = z;
    return;
  }

  /* when u2 != 0 and u2 !=1, can decide if one end of the point is outside depending on inbound value */
  if(u2 != 0 && u2 != 1) {
    if(inbound == 1) { /* goes outside, then temp1->Next_acc is an outside point */
      /* find the next non-intersect point */
      temp2 = temp1->Next_acc;
      if(!temp2) temp2 = list;
      while(temp2->intersect) {
         temp2=temp2->Next_acc;
         if(!temp2) temp2 = list;
      }

      temp2->isInside = 0;
    }
    else if(inbound ==2) { /* goes inside, then temp1 is an outside point */
      temp1->isInside = 0;
    }
  }

  temp2 = temp1->Next_acc;
  while ( temp2 ) {
    if( temp2->intersect == 1 ) {
      if( temp2->u > u_cur ) {
   break;
      }
    }
    else
      break;
    temp1 = temp2;
    temp2 = temp2->Next_acc;
  }

  /* assign value */
  temp = getNext_acc();
  temp->x = x;
  temp->y = y;
  temp->z = z;
  temp->u = u_cur;
  temp->intersect = 1;
  temp->inbound = inbound;
  temp->isInside = 1;
  temp->initialized = 1;
  temp1->Next_acc = temp;
  temp->Next_acc = temp2;

}

double gridArea_acc(struct Node_acc *grid) {
  double x[20], y[20], z[20];
  struct Node_acc *temp=NULL;
  double area;
  int n;

  temp = grid;
  n = 0;
  while( temp ) {
    x[n] = temp->x;
    y[n] = temp->y;
    z[n] = temp->z;
    n++;
    temp = temp->Next_acc;
  }

  area = great_circle_area_acc(n, x, y, z);

  return area;

}

int isIntersect_acc(struct Node_acc node) {

  return node.intersect;

}


int getInbound_acc( struct Node_acc node )
{
  return node.inbound;
}

struct Node_acc *getLast_acc(struct Node_acc *list)
{
  struct Node_acc *temp1;

  temp1 = list;
  if( temp1 ) {
    while( temp1->Next_acc ) {
      temp1 = temp1->Next_acc;
    }
  }

  return temp1;
}


int getFirstInbound_acc( struct Node_acc *list, struct Node_acc *nodeOut)
{
  struct Node_acc *temp=NULL;

  temp=list;

  while(temp) {
    if( temp->inbound == 2 ) {
      copyNode_acc(nodeOut, *temp);
      return 1;
    }
    temp=temp->Next_acc;
  }

  return 0;
}

void getCoordinate_acc(struct Node_acc node, double *x, double *y, double *z)
{


  *x = node.x;
  *y = node.y;
  *z = node.z;

}

void getCoordinates_acc(struct Node_acc *node, double *p)
{


  p[0] = node->x;
  p[1] = node->y;
  p[2] = node->z;

}

void setCoordinate_acc(struct Node_acc *node, double x, double y, double z)
{


  node->x = x;
  node->y = y;
  node->z = z;

}

/* set inbound value for the points in interList that has inbound =0,
   this will also set some inbound value of the points in list1
*/

void setInbound_acc(struct Node_acc *interList, struct Node_acc *list)
{

  struct Node_acc *temp1=NULL, *temp=NULL;
  struct Node_acc *temp1_prev=NULL, *temp1_next=NULL;
  int prev_is_inside, next_is_inside;

  /* for each point in interList, search through list to decide the inbound value the interList point */
  /* For each inbound point, the prev node should be outside and the next is inside. */
  if(length_acc(interList) == 0) return;

  temp = interList;

  while(temp) {
    if( !temp->inbound) {
      /* search in grid1 to find the prev and next point of temp, when prev point is outside and next point is inside
    inbound = 2, else inbound = 1*/
      temp1 = list;
      temp1_prev = NULL;
      temp1_next = NULL;
      while(temp1) {
   if(sameNode_acc(*temp1, *temp)) {
     if(!temp1_prev) temp1_prev = getLast_acc(list);
     temp1_next = temp1->Next_acc;
     if(!temp1_next) temp1_next = list;
     break;
   }
   temp1_prev = temp1;
   temp1 = temp1->Next_acc;
      }
      if(!temp1_next) printf("ERROR from create_xgrid.c: temp is not in list1\n");
      if( temp1_prev->isInside == 0 && temp1_next->isInside == 1)
   temp->inbound = 2;   /* go inside */
      else
   temp->inbound = 1;
    }
    temp=temp->Next_acc;
  }
}

int isInside_acc(struct Node_acc *node) {

  if(node->isInside == -1) printf("ERROR from mosaic_util.c: node->isInside is not set\n");
  return(node->isInside);

}

/*  #define debug_test_create_xgrid */

/* check if node is inside polygon list or not */
 int insidePolygon_acc( struct Node_acc *node, struct Node_acc *list)
{
  int i, ip, is_inside;
  double pnt0[3], pnt1[3], pnt2[3];
  double anglesum;
  struct Node_acc *p1=NULL, *p2=NULL;

  anglesum = 0;

  pnt0[0] = node->x;
  pnt0[1] = node->y;
  pnt0[2] = node->z;

  p1 = list;
  p2 = list->Next_acc;
  is_inside = 0;


  while(p1) {
    pnt1[0] = p1->x;
    pnt1[1] = p1->y;
    pnt1[2] = p1->z;
    pnt2[0] = p2->x;
    pnt2[1] = p2->y;
    pnt2[2] = p2->z;
    if(samePoint_acc(pnt0[0], pnt0[1], pnt0[2], pnt1[0], pnt1[1], pnt1[2])) return 1;
    anglesum += spherical_angle_acc(pnt0, pnt2, pnt1);
    p1 = p1->Next_acc;
    p2 = p2->Next_acc;
    if(p2==NULL)p2 = list;
  }

  if( fabs(anglesum - 2*M_PI) < EPSLN8 )
    is_inside = 1;
  else
    is_inside = 0;

  return is_inside;

}

int inside_a_polygon_acc(double *lon1, double *lat1, int *npts, double *lon2, double *lat2)
{

  double x2[20], y2[20], z2[20];
  double x1, y1, z1;
  double min_x2, max_x2, min_y2, max_y2, min_z2, max_z2;
  int isinside, i;

  struct Node_acc *grid1=NULL, *grid2=NULL;

  /* first convert to cartesian grid */
  latlon2xyz_acc(*npts, lon2, lat2, x2, y2, z2);
  latlon2xyz_acc(1, lon1, lat1, &x1, &y1, &z1);

  max_x2 = maxval_double_acc(*npts, x2);
  if(x1 >= max_x2+RANGE_CHECK_CRITERIA) return 0;
  min_x2 = minval_double_acc(*npts, x2);
  if(min_x2 >= x1+RANGE_CHECK_CRITERIA) return 0;

  max_y2 = maxval_double_acc(*npts, y2);
  if(y1 >= max_y2+RANGE_CHECK_CRITERIA) return 0;
  min_y2 = minval_double_acc(*npts, y2);
  if(min_y2 >= y1+RANGE_CHECK_CRITERIA) return 0;

  max_z2 = maxval_double_acc(*npts, z2);
  if(z1 >= max_z2+RANGE_CHECK_CRITERIA) return 0;
  min_z2 = minval_double_acc(*npts, z2);
  if(min_z2 >= z1+RANGE_CHECK_CRITERIA) return 0;


  /* add x2,y2,z2 to a Node_acc */
  rewindList_acc();
  grid1 = getNext_acc();
  grid2 = getNext_acc();

  addEnd_acc(grid1, x1, y1, z1, 0, 0, 0, -1);
  for(i=0; i<*npts; i++) addEnd_acc(grid2, x2[i], y2[i], z2[i], 0, 0, 0, -1);

  isinside = insidePolygon_acc(grid1, grid2);

  return isinside;

}

/*
  is_near_pole() reuturns 1 if a polygon has any point with a latitude
  within a threshold from the CGS poles (i.e. near +- Pi/2).
   Otherwise returns 0.
*/
int is_near_pole_acc(const double y[], int n) {
  int pole = 0;
  for (int i = 0; i < n; i++) {
    if ((fabs(y[i]) + from_pole_threshold_rad_acc) > M_PI_2 ) {
      pole = 1;
      break;
    }
  }
  return pole;
}

/*
  crosses_pole() returns 1 iff one line segment of the polygon has its end at opposit
  sides of a pole. i.e. if the longitudes are seperated by about Pi. Note, for realistic
  data (not huge polygons), if crosses_pole() reutrns 1, so should is_near_pole().
*/
int crosses_pole_acc(const double x[] , int n) {
  int has_cl = 0;
  for (int i = 0; i < n; i++) {
    int im = (i + n - 1) % n;
    //int  ip = (i + 1) % n;
    double dx = x[i] - x[im];
    if (fabs(dx + M_PI) < SMALL_VALUE || fabs(dx - M_PI) < SMALL_VALUE) {
      has_cl = 1;
      break;
    }
  }
  return has_cl;
}

/*
  Set the_rotation_matrix_acc  global variable.
  The rotation is 45 degrees and about the vector with orign at
  earths center and the direction <0,1,1>/SQRT(2).  I.e. a big
  rotation away from the pole if what is being rotaed is near a pole.
  For rotation matricies formulas and examples, see F.S.Hill, Computer
  Graphics Using OpenGL, @nd ed., Chapter 5.3.
*/
void set_the_rotation_matrix_acc() {
  double is2 = 1.0 /M_SQRT2;

  double m00 = 0;
  double m01 = - is2;
  double m02 = is2;
  double m11 = 1.0/2;
  double m12 = 0.5;

  double m[3][3] = { {m00, m01, m02}, {m02, m11, m12},{m01, m12, m11} };

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      the_rotation_matrix_acc[i][j] = m[i][j];
    }
  }
#pragma acc data update device(the_rotation_matrix_acc[:3][:3])
}

/* Rotate point given the passed in rotation matrix  */
void rotate_point_acc(double rv[], double rmat [][3]) {
  double v[3];

  for (int i = 0; i < 3; i++) {
    v[i] = 0.0;
    for (int j = 0; j < 3; j++) {
      v[i] += rmat[i][j] * rv[j];
    }
  }
  for (int i = 0; i < 3; i++) {
    rv[i] = v[i];
  }
}

/*
  Rotate polygon defined by x[], y[] points and store in xr[], yr[]*/
void rotate_poly_acc(const double x[], const double y[], const int n, double xr[], double yr[]) {
  double sv[2]; //a rotated lat/lon
  double rv[3]; //rotated xyz point
  for (int i = 0; i < n; i++) {
    latlon2xyz_acc(1, &x[i], &y[i], &rv[0], &rv[1], &rv[2]);
    rotate_point_acc(rv, the_rotation_matrix_acc);
    xyz2latlon_acc(1, &rv[0], &rv[1], &rv[2], &sv[0], &sv[1]);
    xr[i] = sv[0];
    yr[i] = sv[1];
  }
}

void pimod_acc(double x[],int nn)
{
  for (int i=0;i<nn;i++) {
    if      (x[i] < -M_PI) x[i] += TPI;
    else if (x[i] >  M_PI) x[i] -= TPI;
  }
}

/*******************************************************************************
 int inside_edge(double x0, double y0, double x1, double y1, double x, double y)
 determine a point(x,y) is inside or outside a given edge with vertex,
 (x0,y0) and (x1,y1). return 1 if inside and 0 if outside. <y1-y0, -(x1-x0)> is
 the outward edge normal from vertex <x0,y0> to <x1,y1>. <x-x0,y-y0> is the vector
 from <x0,y0> to <x,y>.
 if Inner produce <x-x0,y-y0>*<y1-y0, -(x1-x0)> > 0, outside, otherwise inside.
 inner product value = 0 also treate as inside.
*******************************************************************************/
int inside_edge_acc(double x0, double y0, double x1, double y1, double x, double y)
{
   const double SMALL = 1.e-12;
   double product;

   product = ( x-x0 )*(y1-y0) + (x0-x1)*(y-y0);
   return (product<=SMALL) ? 1:0;

 }; /* inside_edge */

/* Intersects between the line a and the seqment s
   where both line and segment are great circle lines on the sphere represented by
   3D cartesian points.
   [sin sout] are the ends of a line segment
   returns true if the lines could be intersected, false otherwise.
   inbound means the direction of (a1,a2) go inside or outside of (q1,q2,q3)
*/
int line_intersect_2D_3D_acc(double *a1, double *a2, double *q1, double *q2, double *q3,
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
  if(samePoint_acc(a1[0], a1[1], a1[2], q1[0], q1[1], q1[2])) {
    *u_a = 0;
    *u_q = 0;
    intersect[0] = a1[0];
    intersect[1] = a1[1];
    intersect[2] = a1[2];
    return 1;
   }
   else if (samePoint_acc(a1[0], a1[1], a1[2], q2[0], q2[1], q2[2])) {
    *u_a = 0;
    *u_q = 1;
    intersect[0] = a1[0];
    intersect[1] = a1[1];
    intersect[2] = a1[2];
    return 1;
  }
   else if(samePoint_acc(a2[0], a2[1], a2[2], q1[0], q1[1], q1[2])) {
    *u_a = 1;
    *u_q = 0;
    intersect[0] = a2[0];
    intersect[1] = a2[1];
    intersect[2] = a2[2];
    return 1;
   }
   else if (samePoint_acc(a2[0], a2[1], a2[2], q2[0], q2[1], q2[2])) {
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
  is_inter1 = intersect_tri_with_line_acc(plane, a1, a2, plane_p, u_a);

  if(!is_inter1)
     return 0;

  if(fabs(*u_a) < EPSLN8) *u_a = 0;
  if(fabs(*u_a-1) < EPSLN8) *u_a = 1;

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
  is_inter2 = intersect_tri_with_line_acc(plane, q1, q2, plane_p, u_q);

  if(!is_inter2)
     return 0;

  if(fabs(*u_q) < EPSLN8) *u_q = 0;
  if(fabs(*u_q-1) < EPSLN8) *u_q = 1;

  if( (*u_q < 0) || (*u_q > 1) ) return 0;

  u =*u_a;

  /* The two planes are coincidental */
  vect_cross_acc(a1, a2, c1);
  vect_cross_acc(q1, q2, c2);
  vect_cross_acc(c1, c2, c3);
  coincident = metric_acc(c3);

  if(fabs(coincident) < EPSLN30) return 0;

  /* Calculate point of intersection */
  intersect[0]=a1[0] + u*(a2[0]-a1[0]);
  intersect[1]=a1[1] + u*(a2[1]-a1[1]);
  intersect[2]=a1[2] + u*(a2[2]-a1[2]);

  norm = metric_acc( intersect );
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

    vect_cross_acc(v1, v2, c1);
    vect_cross_acc(v1, p1, c2);

    sense = dot_acc(c1, c2);
    *inbound = 1;
    if(sense > 0) *inbound = 2; /* v1 going into v2 in CCW sense */
  }

  return 1;
}
