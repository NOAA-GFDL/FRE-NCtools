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

#define min(a,b) (a<b ? a:b)
#define max(a,b) (a>b ? a:b)

struct Node_gpu{
  double x, y, z, u, u_clip;
  int intersect; /* indicate if this point is an intersection, 0 = no, 1= yes, 2=both intersect and vertices */
  int inbound;      /* -1 uninitialized, 0 coincident, 1 outbound, 2 inbound */
  int initialized; /* = 0 means empty list */
  int isInside;   /* = 1 means one point is inside the other polygon, 0 is not, -1 undecided. */
  int subj_index; /* the index of subject point that an intersection follow. */
  int clip_index; /* the index of clip point that an intersection follow */
  struct Node_gpu *Next_gpu;
};

#pragma acc routine seq
int nearest_index_gpu(double value, const double *array, int ia);

#pragma acc routine seq
int lon_fix_gpu(double *x, double *y, int n_in, double tlon);

#pragma acc routine seq
double minval_double_gpu(int size, const double *data);

#pragma acc routine seq
double maxval_double_gpu(int size, const double *data);

#pragma acc routine seq
double avgval_double_gpu(int size, const double *data);

#pragma acc routine seq
void latlon2xyz_gpu(int size, const double *lon, const double *lat, double *x, double *y, double *z);

#pragma acc routine seq
void xyz2latlon_gpu(int size, const double *x, const double *y, const double *z, double *lon, double *lat);

#pragma acc routine seq
double poly_area_gpu(const double lon[], const double lat[], int n);

#pragma acc routine seq
int fix_lon_gpu(double lon[], double lat[], int n, double tlon);

#pragma acc routine seq
double spherical_angle_gpu(const double *v1, const double *v2, const double *v3);

#pragma acc routine seq
void vect_cross_gpu(const double *p1, const double *p2, double *e );

#pragma acc routine seq
double dot_gpu(const double *p1, const double *p2);

#pragma acc routine seq
double metric_gpu(const double *p) ;

#pragma acc routine seq
int intersect_tri_with_line_gpu(const double *plane, const double *l1, const double *l2, double *p, double *t);

#pragma acc routine seq
void mult_gpu(double m[], double v[], double out_v[]);

#pragma acc routine seq
int invert_matrix_3x3_gpu(double m[], double m_inv[]);

#pragma acc routine seq
double great_circle_area_gpu(int n, const double *x, const double *y, const double *z);

#pragma acc routine seq
int insidePolygon_gpu(struct Node_gpu *node, struct Node_gpu *list );

#pragma acc routine seq
int inside_a_polygon_gpu( double *lon1, double *lat1, int *npts, double *lon2, double *lat2);

#pragma acc routine seq
void rewindList_gpu(void);

#pragma acc routine seq
struct Node_gpu *getNext_gpu();

#pragma acc routine seq
void initNode_gpu(struct Node_gpu *node);

#pragma acc routine seq
void addEnd_gpu(struct Node_gpu *list, double x, double y, double z, int intersect, double u, int inbound, int inside);

#pragma acc routine seq
int addIntersect_gpu(struct Node_gpu *list, double x, double y, double z, int intersect, double u1, double u2,
                int inbound, int is1, int ie1, int is2, int ie2);

#pragma acc routine seq
void insertIntersect_gpu(struct Node_gpu *list, double x, double y, double z, double u1, double u2, int inbound,
                     double x2, double y2, double z2);

#pragma acc routine seq
int length_gpu(struct Node_gpu *list);

#pragma acc routine seq
int samePoint_gpu(double x1, double y1, double z1, double x2, double y2, double z2);

#pragma acc routine seq
int sameNode_gpu(struct Node_gpu node1, struct Node_gpu node2);

#pragma acc routine seq
void addNode_gpu(struct Node_gpu *list, struct Node_gpu nodeIn);

#pragma acc routine seq
struct Node_gpu *getNode_gpu(struct Node_gpu *list, struct Node_gpu inNode_gpu);

#pragma acc routine seq
struct Node_gpu *getNextNode_gpu(struct Node_gpu *list);

#pragma acc routine seq
void copyNode_gpu(struct Node_gpu *node_out, struct Node_gpu node_in);

#pragma acc routine seq
void printNode_gpu(struct Node_gpu *list, char *str);

#pragma acc routine seq
int intersectInList_gpu(struct Node_gpu *list, double x, double y, double z);

#pragma acc routine seq
void insertAfter_gpu(struct Node_gpu *list, double x, double y, double z, int intersect, double u, int inbound,
                 double x2, double y2, double z2);

#pragma acc routine seq
double gridArea_gpu(struct Node_gpu *grid);

#pragma acc routine seq
int isIntersect_gpu(struct Node_gpu node);

#pragma acc routine seq
int getInbound_gpu( struct Node_gpu node );

#pragma acc routine seq
struct Node_gpu *getLast_gpu(struct Node_gpu *list);

#pragma acc routine seq
int getFirstInbound_gpu( struct Node_gpu *list, struct Node_gpu *nodeOut);

#pragma acc routine seq
void getCoordinate_gpu(struct Node_gpu node, double *x, double *y, double *z);

#pragma acc routine seq
void getCoordinates_gpu(struct Node_gpu *node, double *p);

#pragma acc routine seq
void setCoordinate_gpu(struct Node_gpu *node, double x, double y, double z);

#pragma acc routine seq
void setInbound_gpu(struct Node_gpu *interList, struct Node_gpu *list);

#pragma acc routine seq
int isInside_gpu(struct Node_gpu *node);

#pragma acc routine seq
void set_rotate_poly_true_gpu(void);

#pragma acc routine seq
int is_near_pole_gpu(const double y[], int n);

#pragma acc routine seq
int crosses_pole_gpu(const double x[], int n);

#pragma acc routine seq
void rotate_point_gpu( double rv[], double rmat [][3]);

#pragma acc routine seq
void rotate_poly_gpu(const double x[], const double y[], const int n, double xr[], double yr[]);

#pragma acc routine seq
void set_the_rotation_matrix_gpu();

#pragma acc routine seq
void pimod_gpu(double x[],int nn);

#pragma acc routine seq
int inside_edge_gpu(double x0, double y0, double x1, double y1, double x, double y);

#pragma acc routine seq
int line_intersect_2D_3D_gpu(double *a1, double *a2, double *q1, double *q2, double *q3,
                         double *intersect, double *u_a, double *u_q, int *inbound);
