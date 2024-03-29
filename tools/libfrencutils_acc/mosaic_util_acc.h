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
#ifndef RANGE_CHECK_CRITERIA
#define RANGE_CHECK_CRITERIA 0.05
#endif

/* override the `fabs` function based on the type */
#define fabs(X) _Generic((X), \
     long double: fabsl, \
          double: fabs)(X)

#define min(a,b) (a<b ? a:b)
#define max(a,b) (a>b ? a:b)
#define SMALL_VALUE ( 1.e-10 )
struct Node_acc{
  double x, y, z, u, u_clip;
  int intersect; /* indicate if this point is an intersection, 0 = no, 1= yes, 2=both intersect and vertices */
  int inbound;      /* -1 uninitialized, 0 coincident, 1 outbound, 2 inbound */
  int initialized; /* = 0 means empty list */
  int isInside;   /* = 1 means one point is inside the other polygon, 0 is not, -1 undecided. */
  int subj_index; /* the index of subject point that an intersection follow. */
  int clip_index; /* the index of clip point that an intersection follow */
  struct Node_acc *Next;
};


void error_handler_acc(const char *msg);

int nearest_index_acc(double value, const double *array, int ia);

int lon_fix_acc(double *x, double *y, int n_in, double tlon);

double minval_double_acc(int size, const double *data);

double maxval_double_acc(int size, const double *data);

double avgval_double_acc(int size, const double *data);

void latlon2xyz_acc(int size, const double *lon, const double *lat, double *x, double *y, double *z);

void xyz2latlon_acc(int size, const double *x, const double *y, const double *z, double *lon, double *lat);

double box_area_acc(double ll_lon, double ll_lat, double ur_lon, double ur_lat);

double poly_area_acc(const double lon[], const double lat[], int n);

double poly_area_dimensionless_acc(const double lon[], const double lat[], int n);

double poly_area_no_adjust_acc(const double x[], const double y[], int n);

int fix_lon_acc(double lon[], double lat[], int n, double tlon);

void tokenize_acc(const char * const string, const char *tokens, unsigned int varlen,
                  unsigned int maxvar, char * pstring, unsigned int * const nstr);

double great_circle_distance_acc(double *p1, double *p2);

double spherical_excess_area_acc(const double* p_ll, const double* p_ul,
                  const double* p_lr, const double* p_ur, double radius);

void vect_cross_acc(const double *p1, const double *p2, double *e );

double spherical_angle_acc(const double *v1, const double *v2, const double *v3);

void normalize_vect_acc(double *e);

void unit_vect_latlon_acc(int size, const double *lon, const double *lat, double *vlon, double *vlat);

double great_circle_area_acc(int n, const double *x, const double *y, const double *z);

double * cross_acc(const double *p1, const double *p2);

double dot_acc(const double *p1, const double *p2);

int intersect_tri_with_line_acc(const double *plane, const double *l1, const double *l2, double *p,
                                double *t);

int invert_matrix_3x3_acc(long double m[], long double m_inv[]);

void mult_acc(long double m[], long double v[], long double out_v[]);

double metric_acc(const double *p);

int insidePolygon_acc(struct Node_acc *node, struct Node_acc *list );

int inside_a_polygon_acc( double *lon1, double *lat1, int *npts, double *lon2, double *lat2);

void rewindList_acc(void);

struct Node_acc *getNext_acc();

void initNode_acc(struct Node_acc *node);

void addEnd_acc(struct Node_acc *list, double x, double y, double z, int intersect, double u, int inbound, int inside);

int addIntersect_acc(struct Node_acc *list, double x, double y, double z, int intersect, double u1, double u2,
                     int inbound, int is1, int ie1, int is2, int ie2);

void insertIntersect_acc(struct Node_acc *list, double x, double y, double z, double u1, double u2, int inbound,
                         double x2, double y2, double z2);

int length_acc(struct Node_acc *list);

int samePoint_acc(double x1, double y1, double z1, double x2, double y2, double z2);

int sameNode_acc_acc(struct Node_acc node1, struct Node_acc node2);

void addNode_acc_acc(struct Node_acc *list, struct Node_acc nodeIn);

struct Node_acc *getNode_acc_acc(struct Node_acc *list, struct Node_acc inNode_acc);

struct Node_acc *getNextNode_acc_acc(struct Node_acc *list);

void copyNode_acc_acc(struct Node_acc *node_out, struct Node_acc node_in);

void printNode_acc_acc(struct Node_acc *list, char *str);

int intersectInList_acc(struct Node_acc *list, double x, double y, double z);

void insertAfter_acc(struct Node_acc *list, double x, double y, double z, int intersect, double u, int inbound,
       double x2, double y2, double z2);

double gridArea_acc(struct Node_acc *grid);

int isIntersect_acc(struct Node_acc node);

int getInbound_acc( struct Node_acc node );

struct Node_acc *getLast_acc(struct Node_acc *list);

int getFirstInbound_acc( struct Node_acc *list, struct Node_acc *nodeOut);

void getCoordinate_acc(struct Node_acc node, double *x, double *y, double *z);

void getCoordinates_acc(struct Node_acc *node, double *p);

void setCoordinate_acc(struct Node_acc *node, double x, double y, double z);

void setInbound_acc(struct Node_acc *interList, struct Node_acc *list);

int isInside_acc(struct Node_acc *node);

void set_reproduce_siena_true_acc(void);

void set_rotate_poly_true_acc(void);

int is_near_pole_acc(const double y[], int n);

int crosses_pole_acc(const double x[], int n);

void rotate_point_acc( double rv[], double rmat [][3]);

void rotate_poly_acc(const double x[], const double y[], const int n,
                     double xr[], double yr[]);

void set_the_rotation_matrix_acc();
