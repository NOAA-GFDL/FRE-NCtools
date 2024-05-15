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

#ifndef _GENERAL_UTILS_H
#define _GENERAL_UTILS_H 1
#endif

/* override the `fabs` function based on the type */
#define fabs(X) _Generic((X), \
          double: fabs)(X)

#define min(a,b) (a<b ? a:b)
#define max(a,b) (a>b ? a:b)

struct Node{
  double x, y, z, u, u_clip;
  int intersect; /* indicate if this point is an intersection, 0 = no, 1= yes, 2=both intersect and vertices */
  int inbound;      /* -1 uninitialized, 0 coincident, 1 outbound, 2 inbound */
  int initialized; /* = 0 means empty list */
  int isInside;   /* = 1 means one point is inside the other polygon, 0 is not, -1 undecided. */
  int subj_index; /* the index of subject point that an intersection follow. */
  int clip_index; /* the index of clip point that an intersection follow */
  struct Node *Next;
};

#pragma acc routine seq
int nearest_index_acc(double value, const double *array, int ia);

#pragma acc routine seq
int lon_fix(double *x, double *y, int n_in, double tlon);

#pragma acc routine seq
double minval_double(int size, const double *data);

#pragma acc routine seq
double maxval_double(int size, const double *data);

#pragma acc routine seq
double avgval_double(int size, const double *data);

#pragma acc routine seq
void latlon2xyz(int size, const double *lon, const double *lat, double *x, double *y, double *z);

#pragma acc routine seq
void xyz2latlon(int size, const double *x, const double *y, const double *z, double *lon, double *lat);

#pragma acc routine seq
double box_area(double ll_lon, double ll_lat, double ur_lon, double ur_lat);

#pragma acc routine seq
double poly_area(const double lon[], const double lat[], int n);

#pragma acc routine seq
double poly_area_dimensionless(const double lon[], const double lat[], int n);

#pragma acc routine seq
double poly_area_no_adjust(const double x[], const double y[], int n);

#pragma acc routine seq
int fix_lon(double lon[], double lat[], int n, double tlon);

#pragma acc routine seq
double great_circle_distance(double *p1, double *p2);

#pragma acc routine seq
double spherical_excess_area(const double* p_ll, const double* p_ul,
                             const double* p_lr, const double* p_ur, double radius);

#pragma acc routine seq
void vect_cross(const double *p1, const double *p2, double *e );

#pragma acc routine seq
double spherical_angle(const double *v1, const double *v2, const double *v3);

#pragma acc routine seq
void normalize_vect(double *e);

#pragma acc routine seq
void unit_vect_latlon(int size, const double *lon, const double *lat, double *vlon, double *vlat);

#pragma acc routine seq
double great_circle_area(int n, const double *x, const double *y, const double *z);

#pragma acc routine seq
double * cross(const double *p1, const double *p2);

#pragma acc routine seq
double dot(const double *p1, const double *p2);

#pragma acc routine seq
int intersect_tri_with_line(const double *plane, const double *l1, const double *l2, double *p, double *t);

#pragma acc routine seq
int invert_matrix_3x3(double m[], double m_inv[]);

#pragma acc routine seq
void mult(double m[], double v[], double out_v[]);

#pragma acc routine seq
double metric(const double *p);

#pragma acc routine seq
int insidePolygon(struct Node *node, struct Node *list );

#pragma acc routine seq
int inside_a_polygon( double *lon1, double *lat1, int *npts, double *lon2, double *lat2);

#pragma acc routine seq
void rewindList(void);

#pragma acc routine seq
struct Node *getNext();

#pragma acc routine seq
void initNode(struct Node *node);

#pragma acc routine seq
void addEnd(struct Node *list, double x, double y, double z, int intersect, double u, int inbound, int inside);

#pragma acc routine seq
int addIntersect(struct Node *list, double x, double y, double z, int intersect, double u1, double u2,
                int inbound, int is1, int ie1, int is2, int ie2);

#pragma acc routine seq
void insertIntersect(struct Node *list, double x, double y, double z, double u1, double u2, int inbound,
                     double x2, double y2, double z2);

#pragma acc routine seq
int length(struct Node *list);

#pragma acc routine seq
int samePoint(double x1, double y1, double z1, double x2, double y2, double z2);

#pragma acc routine seq
int sameNode(struct Node node1, struct Node node2);

#pragma acc routine seq
void addNode(struct Node *list, struct Node nodeIn);

#pragma acc routine seq
struct Node *getNode(struct Node *list, struct Node inNode);

#pragma acc routine seq
struct Node *getNextNode(struct Node *list);

#pragma acc routine seq
void copyNode(struct Node *node_out, struct Node node_in);

#pragma acc routine seq
void printNode(struct Node *list, char *str);

#pragma acc routine seq
int intersectInList(struct Node *list, double x, double y, double z);

#pragma acc routine seq
void insertAfter(struct Node *list, double x, double y, double z, int intersect, double u, int inbound,
                 double x2, double y2, double z2);

#pragma acc routine seq
double gridArea(struct Node *grid);

#pragma acc routine seq
int isIntersect(struct Node node);

#pragma acc routine seq
int getInbound( struct Node node );

#pragma acc routine seq
struct Node *getLast(struct Node *list);

#pragma acc routine seq
int getFirstInbound( struct Node *list, struct Node *nodeOut);

#pragma acc routine seq
void getCoordinate(struct Node node, double *x, double *y, double *z);

#pragma acc routine seq
void getCoordinates(struct Node *node, double *p);

#pragma acc routine seq
void setCoordinate(struct Node *node, double x, double y, double z);

#pragma acc routine seq
void setInbound(struct Node *interList, struct Node *list);

#pragma acc routine seq
int isInside(struct Node *node);

#pragma acc routine seq
void set_reproduce_siena_true(void);

#pragma acc routine seq
void set_rotate_poly_true(void);

#pragma acc routine seq
int is_near_pole(const double y[], int n);

#pragma acc routine seq
int crosses_pole(const double x[], int n);

#pragma acc routine seq
void rotate_point( double rv[], double rmat [][3]);

#pragma acc routine seq
void rotate_poly(const double x[], const double y[], const int n, double xr[], double yr[]);

#pragma acc routine seq
void set_the_rotation_matrix();
