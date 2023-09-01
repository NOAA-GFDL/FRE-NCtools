
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <execution>
#include <ranges>
#include <set>
#include <span>
#include <string>
#include <vector>
#include <format>

#include <array>
#include <algorithm>
#include "BBox3D.h"
#include "Polygon.h"
#include "BoxedObj.h"
//#include "mosaic_util.h"
#include "create_xgrid.h"
#include "create_xgrid.h"

#include "cartesian_product.hpp"

using std::sin;

#include "constant.h"
#ifndef MAXXGRID
#define MAXXGRID (1e7)
#endif
#define MV 50

#define AREA_RATIO_THRESH (1.e-6)
#define MASK_THRESH (0.5)
// #define EPSLN8            (1.e-8)
#define EPSLN30 (1.0e-30)
// #define EPSLN10           (1.0e-10)

constexpr double MU_TOLORENCE{1.e-6};
constexpr double SMALL_VALUE_PA{1.0E-10};

//TODO: place elsewhere
//Note three of these below are negative, and one is positive
constexpr float FILL_VALUE_FLOAT{ -std::numeric_limits<float>::max()};
constexpr double FILL_VALUE_DOUBLE{ -std::numeric_limits<double>::max()};
constexpr int FILL_VALUE_INT{ std::numeric_limits<int>::min() };
constexpr size_t FILL_VALUE_SIZE_T{ std::numeric_limits<size_t>::max() };

//constexpr unsigned int MAX_NN{10}; //MAX_near neighbors per cell

using std::vector;
using BBox_t = nct::BBox3D;
using BPair_t = nct::BoxAndId;
using Poly_t = nct::MeshPolygon<double>;
using Point_t = nct::Point3D<double>;

// These functions are defined here as they were moved from the .C file to the the h
// file in order to use the std parallelization for GPUs with nvc++ 23.x compiers:

extern double maxval_double_gpu(int size, const double *data) {
  int n;
  double maxval;

  maxval = data[0];
  for (n = 1; n < size; n++) {
    if (data[n] > maxval) maxval = data[n];
  }

  return maxval;
}

/*******************************************************************************
  double minval_double(int size, double *data)
  get the minimum value of double array
*******************************************************************************/
extern double minval_double_gpu(int size, const double *data) {
  int n;
  double minval;

  minval = data[0];
  for (n = 1; n < size; n++) {
    if (data[n] < minval) minval = data[n];
  }

  return minval;
}

extern double avgval_double_gpu(int size, const double *data) {
  int n;
  double avgval;

  avgval = 0;
  for (n = 0; n < size; n++) avgval += data[n];
  avgval /= size;

  return avgval;
}

extern void error_handler_gpu(const char *msg) {
#ifdef use_libMPI
  MPI_Abort(MPI_COMM_WORLD, -1);
#else
  fprintf(stderr, "FATAL Error: %s\n", msg);
  // exit(1); exit not usable with nvc++ par execution polify on GPU
#endif
}

extern void error_handler_gpu(const std::string &msg) {
  error_handler_gpu(msg.c_str());
}

void v_print_gpu(double x[], double y[], int n) {
  for (int i = 0; i < n; i++) printf(" %20g   %20g\n", x[i] * R2D, y[i] * R2D);
}

int delete_vtx_gpu(double x[], double y[], int n, int n_del) {
  for (; n_del < n - 1; n_del++) {
    x[n_del] = x[n_del + 1];
    y[n_del] = y[n_del + 1];
  }

  return (n - 1);
}

int insert_vtx_gpu(double x[], double y[], int n, int n_ins, double lon_in, double lat_in) {
  int i;

  for (i = n - 1; i >= n_ins; i--) {
    x[i + 1] = x[i];
    y[i + 1] = y[i];
  }

  x[n_ins] = lon_in;
  y[n_ins] = lat_in;
  return (n + 1);
}

double poly_area_gpu(const double x[], const double y[], int n) {
  double area = 0.0;
  int i;

  for (i = 0; i < n; i++) {
    int ip = (i + 1) % n;
    double dx = (x[ip] - x[i]);
    double lat1, lat2;
    double dy, dat;

    lat1 = y[ip];
    lat2 = y[i];
    if (dx > M_PI) dx = dx - 2.0 * M_PI;
    if (dx < -M_PI) dx = dx + 2.0 * M_PI;
    /*Sides that go through a pole contribute PI regardless of their direction and extent.*/
    if (std::abs(dx + M_PI) < SMALL_VALUE_PA || std::abs(dx - M_PI) < SMALL_VALUE_PA) {
      area += M_PI;
      continue;  // next i
    }
    /*
     We have to be careful in implementing sin(0.5*(lat2-lat1))/(0.5*(lat2-lat1))
     in the limit of lat1=lat2 where it becomes 1. Note that in that limit we get the well known formula
     for the area of a section of sphere bounded by two longitudes and two latitude circles.
    */
    if (std::abs(lat1 - lat2) < SMALL_VALUE_PA) /* cheap area calculation along latitude */
      area -= dx * std::sin(0.5 * (lat1 + lat2));
    else {
      // NOTE: reproduced_sienna option removed
      // This expression is a trig identity with the above reproduce_siena case
      dy = 0.5 * (lat1 - lat2);
      dat = std::sin(dy) / dy;
      area -= dx * std::sin(0.5 * (lat1 + lat2)) * dat;
    }
  }

  if (std::abs(area) > M_PI_2) printf("WARNING poly_area: Large values for area: %19.15f\n", area);
  if (area < 0)
    return -area * RADIUS * RADIUS;
  else
    return area * RADIUS * RADIUS;
}

int fix_lon_gpu(double x[], double y[], int n, double tlon) {
  double x_sum, dx;
  int i, nn = n, pole = 0;

  for (i = 0; i < nn; i++)
    if (std::abs(y[i]) >= M_PI_2 - MU_TOLORENCE) pole = 1;
  if (0 && pole) {
    printf("fixing pole cell\n");
    v_print_gpu(x, y, nn);
    printf("---------");
  }

  /* all pole points must be paired */
  /* The reason is poly_area() function needs a contribution equal to the angle (in radians)
     between the sides that connect to the pole. */
  for (i = 0; i < nn; i++)
    if (std::abs(y[i]) >= M_PI_2 - MU_TOLORENCE) {
      int im = (i + nn - 1) % nn, ip = (i + 1) % nn;

      if (y[im] == y[i] && y[ip] == y[i]) {
        nn = delete_vtx_gpu(x, y, nn, i);
        i--;
      } else if (y[im] != y[i] && y[ip] != y[i]) {
        nn = insert_vtx_gpu(x, y, nn, i, x[i], y[i]);
        i++;
      }
    }
  /* first of pole pair has longitude of previous vertex */
  /* second of pole pair has longitude of subsequent vertex */
  for (i = 0; i < nn; i++)
    if (std::abs(y[i]) >= M_PI_2 - MU_TOLORENCE) {
      int im = (i + nn - 1) % nn, ip = (i + 1) % nn;

      if (y[im] != y[i]) x[i] = x[im];
      if (y[ip] != y[i]) x[i] = x[ip];
    }

  /*If a polygon side passes through a Pole insert twin vertices at the Pole*/
  /*A fix is also directly applied to poly_area to handle this case.*/
  for (i = 0; i < nn; i++) {
    int im = (i + nn - 1) % nn;
    // int ip=(i+1)%nn;
    double dx = x[i] - x[im];
    if (std::abs(dx + M_PI) < SMALL_VALUE_PA || std::abs(dx - M_PI) < SMALL_VALUE_PA) {
      double x1 = x[im];
      double x2 = x[i];
      double ypole = M_PI_2;
      if (y[i] < 0.0) ypole = -M_PI_2;
      nn = insert_vtx_gpu(x, y, nn, i, x2, ypole);
      nn = insert_vtx_gpu(x, y, nn, i, x1, ypole);
      break;
    }
  }
  if (nn)
    x_sum = x[0];
  else
    return (0);
  for (i = 1; i < nn; i++) {
    double dx = x[i] - x[i - 1];

    if (dx < -M_PI)
      dx = dx + TPI;
    else if (dx > M_PI)
      dx = dx - TPI;
    x_sum += (x[i] = x[i - 1] + dx);
  }

  dx = (x_sum / nn) - tlon;
  if (dx < -M_PI)
    for (i = 0; i < nn; i++) x[i] += TPI;
  else if (dx > M_PI)
    for (i = 0; i < nn; i++) x[i] -= TPI;

  if (0 && pole) {
    printf("area=%g\n", poly_area_gpu(x, y, nn));
    v_print_gpu(x, y, nn);
    printf("---------");
  }

  return (nn);
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

double poly_ctrlat_gpu(const double x[], const double y[], int n) {
  double ctrlat = 0.0;
  int i;
  for (i = 0; i < n; i++) {
    int ip = (i + 1) % n;
    double dx = (x[ip] - x[i]);
    double dy, avg_y, hdy;
    double lat1, lat2;
    lat1 = y[ip];
    lat2 = y[i];
    dy = lat2 - lat1;
    hdy = dy * 0.5;
    avg_y = (lat1 + lat2) * 0.5;
    if (dx == 0.0) continue;
    if (dx > M_PI) dx = dx - 2.0 * M_PI;
    if (dx <= -M_PI) dx = dx + 2.0 * M_PI;  // flip sign for dx=-pi to fix huge value see comments in function poly_area

    if (fabs(hdy) < SMALL_VALUE_PA) /* cheap area calculation along latitude */
      ctrlat -= dx * (2 * cos(avg_y) + lat2 * sin(avg_y) - cos(lat1));
    else
      ctrlat -= dx * ((sin(hdy) / hdy) * (2 * cos(avg_y) + lat2 * sin(avg_y)) - cos(lat1));
  }
  if (fabs(ctrlat) > M_PI_2) printf("WARNING poly_ctrlat: Large values for ctrlat: %19.15f\n", ctrlat);
  return (ctrlat * RADIUS * RADIUS);
}
/*------------------------------------------------------------------------------
  double poly_ctrlon(const double x[], const double y[], int n, double clon)
  This routine is used to calculate the lontitude of the centroid.
   ---------------------------------------------------------------------------*/
double poly_ctrlon_gpu(const double x[], const double y[], int n, double clon) {
  double ctrlon = 0.0;
  int i;

  clon = clon;
  for (i = 0; i < n; i++) {
    int ip = (i + 1) % n;
    double phi1, phi2, dphi, lat1, lat2, dphi1, dphi2;
    double f1, f2, fac, fint;
    phi1 = x[ip];
    phi2 = x[i];
    lat1 = y[ip];
    lat2 = y[i];
    dphi = phi1 - phi2;

    if (dphi == 0.0) continue;

    f1 = 0.5 * (cos(lat1) * sin(lat1) + lat1);
    f2 = 0.5 * (cos(lat2) * sin(lat2) + lat2);

    /* this will make sure longitude of centroid is at
       the same interval as the center of any grid */
    if (dphi > M_PI) dphi = dphi - 2.0 * M_PI;
    if (dphi < -M_PI) dphi = dphi + 2.0 * M_PI;
    dphi1 = phi1 - clon;
    if (dphi1 > M_PI) dphi1 -= 2.0 * M_PI;
    if (dphi1 < -M_PI) dphi1 += 2.0 * M_PI;
    dphi2 = phi2 - clon;
    if (dphi2 > M_PI) dphi2 -= 2.0 * M_PI;
    if (dphi2 < -M_PI) dphi2 += 2.0 * M_PI;

    if (fabs(dphi2 - dphi1) < M_PI) {
      ctrlon -= dphi * (dphi1 * f1 + dphi2 * f2) / 2.0;
    } else {
      if (dphi1 > 0.0)
        fac = M_PI;
      else
        fac = -M_PI;
      fint = f1 + (f2 - f1) * (fac - dphi1) / fabs(dphi);
      ctrlon -= 0.5 * dphi1 * (dphi1 - fac) * f1 - 0.5 * dphi2 * (dphi2 + fac) * f2 + 0.5 * fac * (dphi1 + dphi2) * fint;
    }
  }
  return (ctrlon * RADIUS * RADIUS);
}

void pimod_gpu(double x[], int nn) {
  for (int i = 0; i < nn; i++) {
    if (x[i] < -M_PI)
      x[i] += TPI;
    else if (x[i] > M_PI)
      x[i] -= TPI;
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
int inside_edge_gpu(double x0, double y0, double x1, double y1, double x, double y) {
  const double SMALL = 1.e-12;
  double product;

  product = (x - x0) * (y1 - y0) + (x0 - x1) * (y - y0);
  return (product <= SMALL) ? 1 : 0;
}

/*******************************************************************************
   Revise Sutherland-Hodgeman algorithm to find the vertices of the overlapping
   between any two grid boxes. It return the number of vertices for the exchange grid.
*******************************************************************************/

int clip_2dx2d_gpu(const double lon1_in[], const double lat1_in[], int n1_in,
                   const double lon2_in[], const double lat2_in[], int n2_in,
                   double lon_out[], double lat_out[]) {
  double lon_tmp[MV], lat_tmp[MV];
  double lon2_tmp[MV], lat2_tmp[MV];
  double x1_0, y1_0, x1_1, y1_1, x2_0, y2_0, x2_1, y2_1;
  double dx1, dy1, dx2, dy2, determ, ds1, ds2;
  int i_out, n_out, inside_last, inside, i1, i2;
  int gttwopi = 0;
  /* clip polygon with each boundary of the polygon */
  /* We treat lon1_in/lat1_in as clip polygon and lon2_in/lat2_in as subject polygon */
  n_out = n1_in;
  for (i1 = 0; i1 < n1_in; i1++) {
    lon_tmp[i1] = lon1_in[i1];
    lat_tmp[i1] = lat1_in[i1];
    if (lon_tmp[i1] > TPI || lon_tmp[i1] < 0.0) gttwopi = 1;
  }
  for (i2 = 0; i2 < n2_in; i2++) {
    lon2_tmp[i2] = lon2_in[i2];
    lat2_tmp[i2] = lat2_in[i2];
  }
  // Some grid boxes near North Pole are clipped wrong (issue #42 )
  // The following heuristic fix seems to work. Why?
  if (gttwopi) {
    pimod_gpu(lon_tmp, n1_in);
    pimod_gpu(lon2_tmp, n2_in);
  }

  x2_0 = lon2_tmp[n2_in - 1];
  y2_0 = lat2_tmp[n2_in - 1];
  for (i2 = 0; i2 < n2_in; i2++) {
    x2_1 = lon2_tmp[i2];
    y2_1 = lat2_tmp[i2];
    x1_0 = lon_tmp[n_out - 1];
    y1_0 = lat_tmp[n_out - 1];
    inside_last = inside_edge_gpu(x2_0, y2_0, x2_1, y2_1, x1_0, y1_0);
    for (i1 = 0, i_out = 0; i1 < n_out; i1++) {
      x1_1 = lon_tmp[i1];
      y1_1 = lat_tmp[i1];
      if ((inside = inside_edge_gpu(x2_0, y2_0, x2_1, y2_1, x1_1, y1_1)) != inside_last) {
        /* there is intersection, the line between <x1_0,y1_0> and  <x1_1,y1_1>
           should not parallel to the line between <x2_0,y2_0> and  <x2_1,y2_1>
           may need to consider truncation error */
        dy1 = y1_1 - y1_0;
        dy2 = y2_1 - y2_0;
        dx1 = x1_1 - x1_0;
        dx2 = x2_1 - x2_0;
        ds1 = y1_0 * x1_1 - y1_1 * x1_0;
        ds2 = y2_0 * x2_1 - y2_1 * x2_0;
        determ = dy2 * dx1 - dy1 * dx2;
        if (fabs(determ) < EPSLN30) {
          error_handler_gpu(
              "the line between <x1_0,y1_0> and  <x1_1,y1_1> should not parallel to "
              "the line between <x2_0,y2_0> and  <x2_1,y2_1>");
        }
        lon_out[i_out] = (dx2 * ds1 - dx1 * ds2) / determ;
        lat_out[i_out++] = (dy2 * ds1 - dy1 * ds2) / determ;
      }
      if (inside) {
        lon_out[i_out] = x1_1;
        lat_out[i_out++] = y1_1;
      }
      x1_0 = x1_1;
      y1_0 = y1_1;
      inside_last = inside;
    }
    if (!(n_out = i_out)) return 0;
    for (i1 = 0; i1 < n_out; i1++) {
      lon_tmp[i1] = lon_out[i1];
      lat_tmp[i1] = lat_out[i1];
    }
    /* shift the starting point */
    x2_0 = x2_1;
    y2_0 = y2_1;
  }
  return (n_out);
}

int create_xgrid_2dx2d_order2_legacy_gpu(const int nlon_in, const int nlat_in, const int nlon_out, const int nlat_out,
                                         const double *lon_in, const double *lat_in, const double *lon_out,
                                         const double *lat_out,
                                         const double *mask_in, int *i_in, int *j_in, int *i_out, int *j_out,
                                         double *xgrid_area, double *xgrid_clon, double *xgrid_clat) {
#define MAX_V 8
  int nx1, nx2, ny1, ny2, nx1p, nx2p, nxgrid;
  double *area_in, *area_out;
  int nblocks = 1;
  int *istart2 = NULL, *iend2 = NULL;
  int npts_left, nblks_left, pos, m, npts_my, ij;
  double *lon_out_min_list, *lon_out_max_list, *lon_out_avg, *lat_out_min_list, *lat_out_max_list;
  double *lon_out_list, *lat_out_list;
  int *pnxgrid = NULL, *pstart;
  int *pi_in = NULL, *pj_in = NULL, *pi_out = NULL, *pj_out = NULL;
  double *pxgrid_area = NULL, *pxgrid_clon = NULL, *pxgrid_clat = NULL;
  int *n2_list;
  int nthreads, nxgrid_block_max;

  nx1 = nlon_in;
  ny1 = nlat_in;
  nx2 = nlon_out;
  ny2 = nlat_out;
  nx1p = nx1 + 1;
  nx2p = nx2 + 1;

  area_in = (double *)malloc(nx1 * ny1 * sizeof(double));
  area_out = (double *)malloc(nx2 * ny2 * sizeof(double));
  get_grid_area(nlon_in, nlat_in, lon_in, lat_in, area_in);
  get_grid_area(nlon_out, nlat_out, lon_out, lat_out, area_out);

  nthreads = 1;

  nblocks = nthreads;

  istart2 = (int *)malloc(nblocks * sizeof(int));
  iend2 = (int *)malloc(nblocks * sizeof(int));

  pstart = (int *)malloc(nblocks * sizeof(int));
  pnxgrid = (int *)malloc(nblocks * sizeof(int));

  nxgrid_block_max = MAXXGRID / nblocks;

  for (m = 0; m < nblocks; m++) {
    pnxgrid[m] = 0;
    pstart[m] = m * nxgrid_block_max;
  }

  if (nblocks == 1) {
    pi_in = i_in;
    pj_in = j_in;
    pi_out = i_out;
    pj_out = j_out;
    pxgrid_area = xgrid_area;
    pxgrid_clon = xgrid_clon;
    pxgrid_clat = xgrid_clat;
  } else {
    pi_in = (int *)malloc(MAXXGRID * sizeof(int));
    pj_in = (int *)malloc(MAXXGRID * sizeof(int));
    pi_out = (int *)malloc(MAXXGRID * sizeof(int));
    pj_out = (int *)malloc(MAXXGRID * sizeof(int));
    pxgrid_area = (double *)malloc(MAXXGRID * sizeof(double));
    pxgrid_clon = (double *)malloc(MAXXGRID * sizeof(double));
    pxgrid_clat = (double *)malloc(MAXXGRID * sizeof(double));
  }

  npts_left = nx2 * ny2;
  nblks_left = nblocks;
  pos = 0;
  for (m = 0; m < nblocks; m++) {
    istart2[m] = pos;
    npts_my = npts_left / nblks_left;
    iend2[m] = istart2[m] + npts_my - 1;
    pos = iend2[m] + 1;
    npts_left -= npts_my;
    nblks_left--;
  }

  lon_out_min_list = (double *)malloc(nx2 * ny2 * sizeof(double));
  lon_out_max_list = (double *)malloc(nx2 * ny2 * sizeof(double));
  lat_out_min_list = (double *)malloc(nx2 * ny2 * sizeof(double));
  lat_out_max_list = (double *)malloc(nx2 * ny2 * sizeof(double));
  lon_out_avg = (double *)malloc(nx2 * ny2 * sizeof(double));
  n2_list = (int *)malloc(nx2 * ny2 * sizeof(int));
  lon_out_list = (double *)malloc(MAX_V * nx2 * ny2 * sizeof(double));
  lat_out_list = (double *)malloc(MAX_V * nx2 * ny2 * sizeof(double));

  auto ids = std::views::iota(0, (nx2 * ny2));
  std::for_each_n(std::execution::par,
                  ids.begin(), (nx2 * ny2),
                  [=](int ij) {
                    // for (ij = 0; ij < nx2 * ny2; ij++) {
                    int i2, j2, n, n0, n1, n2, n3, n2_in, l;
                    double x2_in[MV], y2_in[MV];
                    i2 = ij % nx2;
                    j2 = ij / nx2;
                    n = j2 * nx2 + i2;
                    n0 = j2 * nx2p + i2;
                    n1 = j2 * nx2p + i2 + 1;
                    n2 = (j2 + 1) * nx2p + i2 + 1;
                    n3 = (j2 + 1) * nx2p + i2;
                    x2_in[0] = lon_out[n0];
                    y2_in[0] = lat_out[n0];
                    x2_in[1] = lon_out[n1];
                    y2_in[1] = lat_out[n1];
                    x2_in[2] = lon_out[n2];
                    y2_in[2] = lat_out[n2];
                    x2_in[3] = lon_out[n3];
                    y2_in[3] = lat_out[n3];

                    lat_out_min_list[n] = minval_double_gpu(4, y2_in);
                    lat_out_max_list[n] = maxval_double_gpu(4, y2_in);
                    n2_in = fix_lon_gpu(x2_in, y2_in, 4, M_PI);
                    if (n2_in > MAX_V) error_handler_gpu("create_xgrid.c: n2_in is greater than MAX_V");
                    lon_out_min_list[n] = minval_double_gpu(n2_in, x2_in);
                    lon_out_max_list[n] = maxval_double_gpu(n2_in, x2_in);
                    lon_out_avg[n] = avgval_double_gpu(n2_in, x2_in);
                    n2_list[n] = n2_in;
                    for (l = 0; l < n2_in; l++) {
                      lon_out_list[n * MAX_V + l] = x2_in[l];
                      lat_out_list[n * MAX_V + l] = y2_in[l];
                    }
                  });

  nxgrid = 0;

  for (m = 0; m < nblocks; m++) {
    //NOTE: outer loop is 1/in
    auto ids = std::views::iota(0, ny1);
    // for(j1=0; j1<ny1; j1++)
    std::for_each_n(std::execution::par, ids.begin(), ny1, [=](int j1) {
      for (int i1 = 0; i1 < nx1; i1++) {
        if (mask_in[j1 * nx1 + i1] > MASK_THRESH) {
          int n0, n1, n2, n3, l, n1_in;
          double lat_in_min, lat_in_max, lon_in_min, lon_in_max, lon_in_avg;
          double x1_in[MV], y1_in[MV], x_out[MV], y_out[MV];

          n0 = j1 * nx1p + i1;
          n1 = j1 * nx1p + i1 + 1;
          n2 = (j1 + 1) * nx1p + i1 + 1;
          n3 = (j1 + 1) * nx1p + i1;
          x1_in[0] = lon_in[n0];
          y1_in[0] = lat_in[n0];
          x1_in[1] = lon_in[n1];
          y1_in[1] = lat_in[n1];
          x1_in[2] = lon_in[n2];
          y1_in[2] = lat_in[n2];
          x1_in[3] = lon_in[n3];
          y1_in[3] = lat_in[n3];
          lat_in_min = minval_double_gpu(4, y1_in);
          lat_in_max = maxval_double_gpu(4, y1_in);
          n1_in = fix_lon_gpu(x1_in, y1_in, 4, M_PI);
          lon_in_min = minval_double_gpu(n1_in, x1_in);
          lon_in_max = maxval_double_gpu(n1_in, x1_in);
          lon_in_avg = avgval_double_gpu(n1_in, x1_in);
          //NOTE: Inner loop in 2 our out
          for (int ij = istart2[m]; ij <= iend2[m]; ij++) {
            int n_out, i2, j2, n2_in;
            double xarea, dx, lon_out_min, lon_out_max;
            double x2_in[MAX_V], y2_in[MAX_V];

            i2 = ij % nx2;
            j2 = ij / nx2;

            if (lat_out_min_list[ij] >= lat_in_max || lat_out_max_list[ij] <= lat_in_min) continue;
            /* adjust x2_in according to lon_in_avg*/
            n2_in = n2_list[ij];
            for (l = 0; l < n2_in; l++) {
              x2_in[l] = lon_out_list[ij * MAX_V + l];
              y2_in[l] = lat_out_list[ij * MAX_V + l];
            }
            lon_out_min = lon_out_min_list[ij];
            lon_out_max = lon_out_max_list[ij];
            dx = lon_out_avg[ij] - lon_in_avg;
            if (dx < -M_PI) {
              lon_out_min += TPI;
              lon_out_max += TPI;
              for (l = 0; l < n2_in; l++) x2_in[l] += TPI;
            } else if (dx > M_PI) {
              lon_out_min -= TPI;
              lon_out_max -= TPI;
              for (l = 0; l < n2_in; l++) x2_in[l] -= TPI;
            }

            /* x2_in should in the same range as x1_in after lon_fix, so no need to
               consider cyclic condition
            */
            if (lon_out_min >= lon_in_max || lon_out_max <= lon_in_min) continue;
            if ((n_out = clip_2dx2d_gpu(x1_in, y1_in, n1_in, x2_in, y2_in, n2_in, x_out, y_out)) > 0) {
              double min_area;
              int nn;
              xarea = poly_area_gpu(x_out, y_out, n_out) * mask_in[j1 * nx1 + i1];
              min_area = std::min(area_in[j1 * nx1 + i1], area_out[j2 * nx2 + i2]);
              if (xarea / min_area > AREA_RATIO_THRESH) {
                pnxgrid[m]++;
                if (pnxgrid[m] >= MAXXGRID / nthreads)
                  error_handler_gpu(
                      "nxgrid is greater than MAXXGRID/nthreads, increase MAXXGRID, decrease nthreads, or increase number of MPI ranks");
                nn = pstart[m] + pnxgrid[m] - 1;
                pxgrid_area[nn] = xarea;
                pxgrid_clon[nn] = poly_ctrlon_gpu(x_out, y_out, n_out, lon_in_avg);
                pxgrid_clat[nn] = poly_ctrlat_gpu(x_out, y_out, n_out);
                pi_in[nn] = i1;
                pj_in[nn] = j1;
                pi_out[nn] = i2;
                pj_out[nn] = j2;
              }
            }
          }
        }
      }
    });
  }

  /*copy data if nblocks > 1 */
  if (nblocks == 1) {
    nxgrid = pnxgrid[0];
    pi_in = NULL;
    pj_in = NULL;
    pi_out = NULL;
    pj_out = NULL;
    pxgrid_area = NULL;
    pxgrid_clon = NULL;
    pxgrid_clat = NULL;
  } else {
    int nn, i;
    nxgrid = 0;
    for (m = 0; m < nblocks; m++) {
      for (i = 0; i < pnxgrid[m]; i++) {
        nn = pstart[m] + i;
        i_in[nxgrid] = pi_in[nn];
        j_in[nxgrid] = pj_in[nn];
        i_out[nxgrid] = pi_out[nn];
        j_out[nxgrid] = pj_out[nn];
        xgrid_area[nxgrid] = pxgrid_area[nn];
        xgrid_clon[nxgrid] = pxgrid_clon[nn];
        xgrid_clat[nxgrid] = pxgrid_clat[nn];
        nxgrid++;
      }
    }
    free(pi_in);
    free(pj_in);
    free(pi_out);
    free(pj_out);
    free(pxgrid_area);
    free(pxgrid_clon);
    free(pxgrid_clat);
  }

  free(area_in);
  free(area_out);
  free(lon_out_min_list);
  free(lon_out_max_list);
  free(lat_out_min_list);
  free(lat_out_max_list);
  free(lon_out_avg);
  free(n2_list);
  free(lon_out_list);
  free(lat_out_list);

  return nxgrid;
}

// grid suffix synonyms: '2'/out/target
// grid suffix synonyms: '1'/in/source
int create_xgrid_2dx2d_order2_noahack(const int nlon_in, const int nlat_in, const int nlon_out, const int nlat_out,
                              const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                              const double *mask_in, int *i_in, int *j_in, int *i_out, int *j_out,
                              double *xgrid_area, double *xgrid_clon, double *xgrid_clat)
{

#define MAX_V 8
  double *area_in, *area_out;
  double *lon_out_min_list,*lon_out_max_list,*lon_out_avg,*lat_out_min_list,*lat_out_max_list;
  double *lon_out_list, *lat_out_list;
  int    *n2_list;
  int mxxgrid;

  const int nx1 {nlon_in};
  const int ny1 {nlat_in};
  const int nx2 {nlon_out};
  const int ny2 {nlat_out};
  const int nx1p {nx1 + 1};
  const int nx2p {nx2 + 1};
  mxxgrid = MAXXGRID;

  area_in  = (double *)malloc(nx1*ny1*sizeof(double));
  area_out = (double *)malloc(nx2*ny2*sizeof(double));
  get_grid_area(nlon_in, nlat_in, lon_in, lat_in, area_in);
  get_grid_area(nlon_out, nlat_out, lon_out, lat_out, area_out);

  lon_out_min_list = (double *)malloc(nx2*ny2*sizeof(double));
  lon_out_max_list = (double *)malloc(nx2*ny2*sizeof(double));
  lat_out_min_list = (double *)malloc(nx2*ny2*sizeof(double));
  lat_out_max_list = (double *)malloc(nx2*ny2*sizeof(double));
  lon_out_avg = (double *)malloc(nx2*ny2*sizeof(double));
  n2_list     = (int *)malloc(nx2*ny2*sizeof(int));
  lon_out_list = (double *)malloc(MAX_V*nx2*ny2*sizeof(double));
  lat_out_list = (double *)malloc(MAX_V*nx2*ny2*sizeof(double));

  auto ids = std::views::iota(0, (nx2 * ny2));
  std::for_each_n(std::execution::par, ids.begin(), (nx2 * ny2), [=](int ij) {
    //originals equivalent had been labeled 'acc loop independent'
    int i2, j2, n, n0, n1, n2, n3, n2_in, l;
    double x2_in[MV], y2_in[MV];
    i2 = ij % nx2;
    j2 = ij / nx2;
    n = j2 * nx2 + i2;
    n0 = j2 * nx2p + i2;
    n1 = j2 * nx2p + i2 + 1;
    n2 = (j2 + 1) * nx2p + i2 + 1;
    n3 = (j2 + 1) * nx2p + i2;

    x2_in[0] = lon_out[n0]; y2_in[0] = lat_out[n0];
    x2_in[1] = lon_out[n1]; y2_in[1] = lat_out[n1];
    x2_in[2] = lon_out[n2]; y2_in[2] = lat_out[n2];
    x2_in[3] = lon_out[n3]; y2_in[3] = lat_out[n3];

    lat_out_min_list[n] = minval_double_gpu(4, y2_in);
    lat_out_max_list[n] = maxval_double_gpu(4, y2_in);
    n2_in = fix_lon_gpu(x2_in, y2_in, 4, M_PI);
    //if(n2_in > MAX_V) error_handler("create_xgrid.c: n2_in is greater than MAX_V");
    lon_out_min_list[n] = minval_double_gpu(n2_in, x2_in);
    lon_out_max_list[n] = maxval_double_gpu(n2_in, x2_in);
    lon_out_avg[n] = avgval_double_gpu(n2_in, x2_in);
    n2_list[n] = n2_in;

    //original had been labeled `acc loop independent`
    //consider moving in its own for_each?
    for (l = 0; l < n2_in; l++) {
      lon_out_list[n * MAX_V + l] = x2_in[l];
      lat_out_list[n * MAX_V + l] = y2_in[l];
    }
  });

    // original used acc shared construct here
  // acc loop independent reduction(+:nxgrid) collapse(2)
    //These are the loops over source cells
    int nxgrid;
    ids = std::views::iota(0, ny1);
    std::for_each_n(std::execution::par, ids.begin(), ny1, [=](int j1) {
      for(int  i1=0; i1<nx1; i1++)  //TODO: place i1 in 2D loop (of i1 x j1)
        if( mask_in[j1*nx1+i1] > MASK_THRESH ) {
          int n0, n1, n2, n3, l,n1_in;
          double lat_in_min,lat_in_max,lon_in_min,lon_in_max,lon_in_avg;
          double x1_in[MV], y1_in[MV], x_out[MV], y_out[MV];

          n0 = j1*nx1p+i1;       n1 = j1*nx1p+i1+1;
          n2 = (j1+1)*nx1p+i1+1; n3 = (j1+1)*nx1p+i1;
          x1_in[0] = lon_in[n0]; y1_in[0] = lat_in[n0];
          x1_in[1] = lon_in[n1]; y1_in[1] = lat_in[n1];
          x1_in[2] = lon_in[n2]; y1_in[2] = lat_in[n2];
          x1_in[3] = lon_in[n3]; y1_in[3] = lat_in[n3];
          lat_in_min = minval_double_gpu(4, y1_in);
          lat_in_max = maxval_double_gpu(4, y1_in);
          n1_in = fix_lon_gpu(x1_in, y1_in, 4, M_PI);
          lon_in_min = minval_double_gpu(n1_in, x1_in);
          lon_in_max = maxval_double_gpu(n1_in, x1_in);
          lon_in_avg = avgval_double_gpu(n1_in, x1_in);
          // original had acc loop independent reduction(+:nxgrid)
          for(int ij=0; ij<=nx2*ny2; ij++) { //and the loop over target cells
            int n_out, i2, j2, n2_in;
            double dx, lon_out_min, lon_out_max;
            double x2_in[MAX_V], y2_in[MAX_V];

            i2 = ij%nx2;
            j2 = ij/nx2;

            //Note: longitude gate considered a few lines below this.
            if(lat_out_min_list[ij] >= lat_in_max || lat_out_max_list[ij] <= lat_in_min ) continue;
            /* adjust x2_in according to lon_in_avg*/
            n2_in = n2_list[ij];
            // #pragma acc loop seq
            for(l=0; l<n2_in; l++) {
              x2_in[l] = lon_out_list[ij*MAX_V+l];
              y2_in[l] = lat_out_list[ij*MAX_V+l];
            }
            lon_out_min = lon_out_min_list[ij];
            lon_out_max = lon_out_max_list[ij];
            dx = lon_out_avg[ij] - lon_in_avg;
            if(dx < -M_PI ) {
              lon_out_min += TPI;
              lon_out_max += TPI;
              //#pragma acc loop seq
              for (l=0; l<n2_in; l++) x2_in[l] += TPI;
            }else if (dx >  M_PI) {
              lon_out_min -= TPI;
              lon_out_max -= TPI;
              //#pragma acc loop seq
              for (l=0; l<n2_in; l++) x2_in[l] -= TPI;
            }

            //NOTE: latitudes gate considered a few lines above this.
            if(lon_out_min >= lon_in_max || lon_out_max <= lon_in_min ) continue;
            n_out = 1;
            if (  (n_out = clip_2dx2d( x1_in, y1_in, n1_in, x2_in, y2_in, n2_in, x_out, y_out )) > 0) {
              double xarea = poly_area_gpu (x_out, y_out, n_out ) * mask_in[j1*nx1+i1];
              double min_area = std::min(area_in[j1*nx1+i1], area_out[j2*nx2+i2]);
              if( xarea/min_area > AREA_RATIO_THRESH ) {
//            if(nxgrid>= MAXXGRID/nthreads)
//	      error_handler("nxgrid is greater than MAXXGRID/nthreads, increase MAXXGRID, decrease nthreads, or increase number of MPI ranks");
                xgrid_area[nxgrid] = xarea;
                xgrid_clon[nxgrid] = poly_ctrlon(x_out, y_out, n_out, lon_in_avg);
                xgrid_clat[nxgrid] = poly_ctrlat (x_out, y_out, n_out );
                i_in[nxgrid]       = i1;
                j_in[nxgrid]       = j1;
                i_out[nxgrid]      = i2;
                j_out[nxgrid]      = j2;
                //nxgrid++; TODO: needs to accumulate
              }
            }
          }
        }
  });

  free(area_in);
  free(area_out);
  free(lon_out_min_list);
  free(lon_out_max_list);
  free(lat_out_min_list);
  free(lat_out_max_list);
  free(lon_out_avg);
  free(n2_list);
  free(lon_out_list);
  free(lat_out_list);

  return nxgrid;

}

bool is_near(float x , float y){
    return ((x-y) * (x-y) < 1.5);
}

//TODO: This formula for the maximum number of cells in one grid to any other in another grid is
// merely a heuristic. A better one might involve max and min cell areas. Consider how many
//cells from one grid can overlap the largest cell from another. The factor of 1000 below
//should take care of really wierd cells.
size_t get_max_cell_nns(const size_t nx1, const size_t nx2, const size_t ny1, const size_t ny2)
{ return 1000 * std::max(nx1 * ny1, nx2 * ny2) / std::min(nx1 * ny1, nx2 * ny2); }

//Return the max number of near neighbors for a grid.
size_t
get_max_grid_nns(const size_t nx1, const size_t nx2, const size_t ny1, const size_t ny2)
{ return 10 * std::max(nx1 * ny1, nx2 * ny2); }

/*
 * create_xgrid_2dx2d_order2 - brute force with bounding box usage
 */
void  create_xgrid_2dx2d_order2_bfwbb(const int nlon_in, const int nlat_in, const int nlon_out, const int nlat_out,
                                   const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                                   const double *mask_in, vector<size_t>& i_in, vector<size_t>& j_in,
                                   vector<size_t>& i_out, vector<size_t>& j_out, vector<double>& xgrid_area,
                                   vector<double>& xgrid_clon, vector<double>& xgrid_clat) {
#define MAX_V 8
  const size_t nx1 {(size_t) nlon_in}, nx2{(size_t)nlon_out}, ny1{(size_t)nlat_in}, ny2{(size_t)nlat_out};
  const size_t nx1p{nx1 + 1};
  const size_t nx2p{nx2 + 1};

  std::vector<double> area_in(nx1 * ny1);
  std::vector<double> area_out(nx2 * ny2);
  //TODO: note get_grid_area can be parallelized.
  get_grid_area(nlon_in, nlat_in, lon_in, lat_in, area_in.data());
  get_grid_area(nlon_out, nlat_out, lon_out, lat_out, area_out.data());

  std::cout << "*** S1 ***" <<std::endl;

  // grid suffix synonyms: '2'/out/target
  // grid suffix synonyms: '1'/in/source
  //The out our target pairs
  vector<Poly_t> polys_2;
  vector<BBox_t> boxes_2;
  vector<BPair_t> bPairs_2;
  boxes_2.reserve(nx2 * ny2);
  //bPairs_2.reserve(nx2 * ny2);
  //TODO: Parallelize
  for (size_t ij = 0; ij < nx2 * ny2; ij++) {
    auto i2 = ij % nx2;
    auto j2 = ij / nx2;
    auto n = j2 * nx2 + i2;  //NOTE: using nx2; not nx2p
    assert (n == ij);
    auto ip = get_cell_idxs_ccw_4(i2, j2, nx2p);
    boxes_2.emplace_back(getBoxForSphericalPolygon(lat_out, lon_out, ip));
    //bPairs_2.emplace_back(ij, &boxes_2[n] );
  }

  //The "in" or source pairs
  //TODO::parallelize
  vector<Poly_t> polys_1;
  vector<BBox_t> boxes_1;
  vector<BPair_t> bPairs_1;
  boxes_1.reserve(nx1 * ny1);
  //bPairs_1.reserve(nx1 * ny1);
  for (size_t ij = 0; ij < nx1 * ny1; ij++) {
    auto i1 = ij % nx1;
    auto j1 = ij / nx1;
    auto n = j1 * nx1 + i1;
    assert (n == ij);
    auto ip = get_cell_idxs_ccw_4(i1, j1, nx1p);
    boxes_1.emplace_back(getBoxForSphericalPolygon(lat_in, lon_in, ip));
    //bPairs_1.emplace_back(ij, &boxes_1[n] );
  }
  std::cout << "*** S2 ***" <<std::endl;
  std::cout << "boxes_1.size = " << boxes_1.size()<< std::endl;
  std::cout << "boxes_2.size = " << boxes_2.size()<< std::endl;

  const size_t max_cell_nns = get_max_cell_nns(nx1, nx2, ny1, ny2);
  const size_t max_grid_nns  = get_max_grid_nns(nx1, nx2, ny1, ny2);

  std::cout << "*** S3 ***" <<std::endl;

  //vector<int> nn1(max_cell_nns, FILL_VALUE_INT);
  vector<std::tuple<int, int>> nn_pairs1;
  for (int j = 0; j < boxes_1.size(); j++) {//NOTE:outer loop over boxes 1
    vector<int> nn1;
    auto g2_ids = std::views::iota(0, (int) (nx2 * ny2));
    std::copy_if(std::execution::seq,  //TODO: change 'seq' to 'par'
            g2_ids.begin(), g2_ids.end(), //The range of pairs from ids set
            std::back_inserter(nn1),
            [=, boxes_1 = boxes_1.data(), boxes_2 = boxes_2.data()](auto idx) {
             return (nct::BBox3D::intersect(boxes_1[j], boxes_2[idx]));
            });

    //For reproducibility, results are ordered by increasing ij value:
    vector<int> nn1b;
    for (size_t k = 0; k <nn1.size(); k++) {
      if (nn1[k] != FILL_VALUE_INT) {
        nn1b.push_back(nn1[k]);
        nn1[k] = FILL_VALUE_INT;
      }else{
        break;//TODO: and NOTE: Nowhere in the copy_if documentation does it say that
        //the inserts will be at the front of array if fixed length array is used; but it seems to be true.
      }
    }
    sort(nn1b.begin(), nn1b.end(), [](auto x, auto y) {return (x<y);});
    //std::cout << " j nn1.size nns_pairs.size: " << j << " " <<
   // nn1.size() << " " << nn_pairs1.size() << std::endl;
     for (size_t k = 0; k <nn1b.size(); k++) {
       if (nn1b[k] != FILL_VALUE_INT) {
         nn_pairs1.push_back({j, nn1b[k]});
         nn1b[k] = FILL_VALUE_INT;
       }else{
       break;
       }
     }
     nn1.clear();
     nn1b.clear();


  }

  //TRIM nn_pairs to only hold the real pairs.
  vector<std::tuple<int, int>> nn_pairs;
  for (auto &p: nn_pairs1) {
    auto [ij1, ij2] = p;
    //auto ix = get<0>(p); auto jy = get<1>(p);
    if (ij1 != FILL_VALUE_INT) {
      nn_pairs.push_back(p);
    } else {
    break;
    }
  }
  //std::cout << " nn_pairs.size() : " << nn_pairs.size() << std::endl;

  //Given the initial neighbor pairs, calculate the final pair intersections.
  std::vector<double> xarea_v(max_grid_nns, FILL_VALUE_DOUBLE);
  std::vector<double> clon_v(max_grid_nns, FILL_VALUE_DOUBLE);
  std::vector<double> clat_v(max_grid_nns, FILL_VALUE_DOUBLE);
  vector<size_t>i_in_v(max_grid_nns, FILL_VALUE_INT);
  vector<size_t>j_in_v(max_grid_nns, FILL_VALUE_INT);
  vector<size_t>i_out_v(max_grid_nns, FILL_VALUE_INT);
  vector<size_t>j_out_v(max_grid_nns, FILL_VALUE_INT);
  //TODO: parallelize
  for(int ij = 0; ij < nn_pairs.size(); ij++) {
    auto [ij1, ij2] = nn_pairs[ij];
    auto i1 = ij1 % nx1;
    auto j1 = ij1 / nx1;
    auto i2 = ij2 % nx2;
    auto j2 = ij2 / nx2;

    //TODO: refactor into a function ?
    double x1_in[MAX_V], y1_in[MAX_V], x2_in[MAX_V], y2_in[MAX_V], x_out[MAX_V], y_out[MAX_V];
    auto idxs_1 = get_cell_idxs_ccw_4(i1, j1, nx1p);
    for (int i = 0; i < 4; i++) {
      x1_in[i] = lon_in[idxs_1[i]];
      y1_in[i] = lat_in[idxs_1[i]];
    }

    auto idxs_2 = get_cell_idxs_ccw_4(i2, j2, nx2p);
    //The legacy polygon lat-lon representation:
    for (int i = 0; i < 4; i++) {
      x2_in[i] = lon_out[idxs_2[i]];
      y2_in[i] = lat_out[idxs_2[i]];
    }

    //Some polys require "fixing" before calling area (and clipping?) function.
    auto n1_in = fix_lon_gpu(x1_in, y1_in, 4, M_PI);
    auto n2_in = fix_lon_gpu(x2_in, y2_in, 4, M_PI); //TODO: does this one get fixed
    auto lon_in_avg = avgval_double_gpu(n1_in, x1_in);

    //Call the 2D_by_2D clipping algorithm
    auto n_out = clip_2dx2d(x1_in, y1_in, n1_in, x2_in,
                            y2_in, n2_in, x_out, y_out);
    if (n_out > 0) {
      auto xarea = poly_area_gpu(x_out, y_out, n_out) * mask_in[j1 * nx1 + i1];
      auto min_area = std::min(area_in[j1 * nx1 + i1], area_out[j2 * nx2 + i2]);
      if (xarea / min_area > AREA_RATIO_THRESH) {
        xarea_v[ij] = xarea;
        clon_v[ij] = poly_ctrlon(x_out, y_out, n_out, lon_in_avg);
        clat_v[ij] = poly_ctrlat(x_out, y_out, n_out);
        i_in_v[ij] = i1; //stores the lower corner of poly
        j_in_v[ij] = j1;
        i_out_v[ij] = i2;
        j_out_v[ij] = j2;
      } //else leave in the FILL_VALUE
    }
  }
  std::cout << "*** S5 ***" <<std::endl;

  //Save final answers (though copy_if can be used again); but don't use more memory than needed.
    auto nxgrid = count_if(xarea_v.begin(), xarea_v.end(), [](double x) { return x > 0; });

    xgrid_area.reserve(nxgrid);  //TODO: is assign(capacity(, FILL_VALUE) needed?
    i_in.reserve(nxgrid);
    j_in.reserve(nxgrid);
    i_out.reserve(nxgrid);
    j_out.reserve(nxgrid);
    xgrid_clon.reserve(nxgrid);
    xgrid_clat.reserve(nxgrid);

    //auto ir = 0;
    for(int i = 0; i  < nn_pairs.size(); i ++) {
      if (xarea_v[i] > 0.0) {
        xgrid_area.emplace_back(xarea_v[ i ]);
        xgrid_clon.emplace_back(clon_v[i]);
        xgrid_clat.emplace_back(clat_v[i]);
        i_in.emplace_back(i_in_v[i]);
        j_in.emplace_back(j_in_v[i]);
        i_out .emplace_back(i_out_v[i]);
        j_out.emplace_back(j_out_v[i]);
      ///  ++i;
      }
    }

    std::cout <<"create_xgrid_2dx2d_order2_ws end; xgrid_are.size= " << xgrid_area.size() <<std::endl;
   // return nxgrid;
}
/*
 *
 //vector<std::tuple<int, int>> nn_pairs1(
   //       max_grid_nns,
  //        std::tuple<int, int>{FILL_VALUE_INT, FILL_VALUE_INT});
  //auto g1_ids = std::views::iota(0, (int) (nx1 * ny1));
  //auto g2_ids = std::views::iota(0, (int) (nx2 * ny2));
  //auto gxg_ids = tl::views::cartesian_product(g1_ids, g2_ids);



  std::copy_if(
          std::execution::seq,  //TODO: change to 'par' (seq-quential is the default)
          gxg_ids.begin(), gxg_ids.end(), //The range of pairs from ids set
          nn_pairs1.begin(), //GPU nvc++ requires random_access iter; container resizing not allowed.
          //for this item, an access by iterator idiom
          [=, b1v = boxes_1.data(), b2v = boxes_2.data()](auto idx) { //[] is start of lambda definition.
            auto [ij1, ij2] = idx;
            return (nct::BBox3D::intersect(b1v[ ij1], b2v [ ij2]));
          });
  std::cout << "*** S4 ***" <<std::endl;

  // Make the indexes of cells of each grid, and the cartesian product between the two.
  // Note the cartesian product is actually a lazy generator of indices and does not
  // (fortunately) take up N1 x N2 space. The motivation is to use copy_if with index access idiom.
  // Also by constraints from the nvc++ (for GPUs), container resizing is not allowed (so its space
  // must be pre-allocated) and container must feature a random_access iterator.

  vector<std::tuple<int, int>> nn_pairs1(  //Allocate and initialize 1D array of Pairs of Near Neighbors:
          max_grid_nns,
          std::tuple<int, int>{FILL_VALUE_INT, FILL_VALUE_INT});
  //auto g1_ids = std::views::iota(0, (int) (nx1 * ny1));
  //auto g2_ids = std::views::iota(0, (int) (nx2 * ny2));
  //auto gxg_ids = tl::views::cartesian_product(g1_ids, g2_ids);

  std::cout << "*** S3 ***" <<std::endl;

  // Determine the near-neighbors pairs and copy the indexes of each pair. The near-neighbors in this step
  // are merely a superset of the final near-neighbors, and will later be further filtered by clip-polygon.
  // In this step the near-neighbor intersection test is the box-box intersection test. If the final
  // intersection test (with clip_polygon) were used, only about (k * max(N1, N2)), where k is small (~4).
  // If instead polygon boxes are calculated fairly tight around the polygons, that 'k' should only be
  // slightly bigger.
  /*std::copy_if(
          std::execution::seq,  //TODO: change to 'par' (seq-quential is the default)
          gxg_ids.begin(), gxg_ids.end(), //The range of pairs from ids set
          nn_pairs1.begin(), //GPU nvc++ requires random_access iter; container resizing not allowed.
          //for this item, an access by iterator idiom
          [=, b1v = boxes_1.data(), b2v = boxes_2.data()](auto idx) { //[] is start of lambda definition.
            auto [ij1, ij2] = idx;
            return (nct::BBox3D::intersect(b1v[ ij1], b2v [ ij2]));
          });
  std::cout << "*** S4 ***" <<std::endl;



*/
