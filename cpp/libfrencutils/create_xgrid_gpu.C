
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

#include "create_xgrid.h"

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

int create_xgrid_2dx2d_order2_noahack(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
                              const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                              const double *mask_in, int *i_in, int *j_in, int *i_out, int *j_out,
                              double *xgrid_area, double *xgrid_clon, double *xgrid_clat)
{

#define MAX_V 8
  int nx1, nx2, ny1, ny2, nx1p, nx2p, nxgrid;
  double *area_in, *area_out;
  int ij, i1, j1;
  double *lon_out_min_list,*lon_out_max_list,*lon_out_avg,*lat_out_min_list,*lat_out_max_list;
  double *lon_out_list, *lat_out_list;
  int    *n2_list;
  int mxxgrid;

  nx1 = *nlon_in;
  ny1 = *nlat_in;
  nx2 = *nlon_out;
  ny2 = *nlat_out;
  nx1p = nx1 + 1;
  nx2p = nx2 + 1;
  mxxgrid = get_maxxgrid();

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
#if defined(_OPENMP)
  #pragma omp parallel for default(none) shared(nx2,ny2,nx2p,lon_out,lat_out,lat_out_min_list, \
                                              lat_out_max_list,lon_out_min_list,lon_out_max_list, \
                                              lon_out_avg,n2_list,lon_out_list,lat_out_list)
#endif
  nxgrid = 0;
#pragma acc kernels copyin(lon_out[0:(nx2+1)*(ny2+1)], lat_out[0:(nx2+1)*(ny2+1)], mask_in[0:nx1*ny1], \
                           xgrid_area[0:mxxgrid], xgrid_clon[0:mxxgrid], xgrid_clat[0:mxxgrid], \
                           i_in[0:mxxgrid], j_in[0:mxxgrid], i_out[0:mxxgrid],j_out[0:mxxgrid], \
                           area_in[0:nx1*ny1], area_out[0:nx2*ny2], \
                           lon_in[0:(nx1+1)*(ny1+1)], lat_in[0:(nx1+1)*(ny1+1)]) \
                   copyout(lon_out_list[0:MAX_V*nx2*ny2], lat_out_list[0:MAX_V*nx2*ny2], \
		          lat_out_min_list[0:nx2*ny2], lat_out_max_list[0:nx2*ny2],  \
		          lon_out_min_list[0:nx2*ny2], lon_out_max_list[0:nx2*ny2],  \
		          lon_out_avg[0:nx2*ny2], n2_list[0:nx2*ny2]) \
                   copy (nxgrid)
  {
#pragma acc loop independent
    for(ij=0; ij<nx2*ny2; ij++){
      int i2, j2, n, n0, n1, n2, n3, n2_in, l;
      double x2_in[MV], y2_in[MV];
      i2 = ij%nx2;
      j2 = ij/nx2;
      n = j2*nx2+i2;
      n0 = j2*nx2p+i2; n1 = j2*nx2p+i2+1;
      n2 = (j2+1)*nx2p+i2+1; n3 = (j2+1)*nx2p+i2;
      x2_in[0] = lon_out[n0]; y2_in[0] = lat_out[n0];
      x2_in[1] = lon_out[n1]; y2_in[1] = lat_out[n1];
      x2_in[2] = lon_out[n2]; y2_in[2] = lat_out[n2];
      x2_in[3] = lon_out[n3]; y2_in[3] = lat_out[n3];

      lat_out_min_list[n] = minval_double(4, y2_in);
      lat_out_max_list[n] = maxval_double(4, y2_in);
      n2_in = fix_lon(x2_in, y2_in, 4, M_PI);
//    if(n2_in > MAX_V) error_handler("create_xgrid.c: n2_in is greater than MAX_V");
      lon_out_min_list[n] = minval_double(n2_in, x2_in);
      lon_out_max_list[n] = maxval_double(n2_in, x2_in);
      lon_out_avg[n] = avgval_double(n2_in, x2_in);
      n2_list[n] = n2_in;
#pragma acc loop independent
      for(l=0; l<n2_in; l++) {
        lon_out_list[n*MAX_V+l] = x2_in[l];
        lat_out_list[n*MAX_V+l] = y2_in[l];
      }
    }


#if defined(_OPENMP)
    #pragma omp parallel for default(none) shared(nblocks,nx1,ny1,nx1p,mask_in,lon_in,lat_in, \
                                              istart2,iend2,nx2,lat_out_min_list,lat_out_max_list, \
                                              n2_list,lon_out_list,lat_out_list,lon_out_min_list, \
                                              lon_out_max_list,lon_out_avg,area_in,area_out, \
                                              pxgrid_area,pnxgrid,pxgrid_clon,pxgrid_clat,pi_in, \
                                              pj_in,pi_out,pj_out,pstart,nthreads)
#endif
#pragma acc loop independent reduction(+:nxgrid) collapse(2)
    for(j1=0; j1<ny1; j1++) for(i1=0; i1<nx1; i1++) if( mask_in[j1*nx1+i1] > MASK_THRESH ) {
          int n0, n1, n2, n3, l,n1_in;
          double lat_in_min,lat_in_max,lon_in_min,lon_in_max,lon_in_avg;
          double x1_in[MV], y1_in[MV], x_out[MV], y_out[MV];

          n0 = j1*nx1p+i1;       n1 = j1*nx1p+i1+1;
          n2 = (j1+1)*nx1p+i1+1; n3 = (j1+1)*nx1p+i1;
          x1_in[0] = lon_in[n0]; y1_in[0] = lat_in[n0];
          x1_in[1] = lon_in[n1]; y1_in[1] = lat_in[n1];
          x1_in[2] = lon_in[n2]; y1_in[2] = lat_in[n2];
          x1_in[3] = lon_in[n3]; y1_in[3] = lat_in[n3];
          lat_in_min = minval_double(4, y1_in);
          lat_in_max = maxval_double(4, y1_in);
          n1_in = fix_lon(x1_in, y1_in, 4, M_PI);
          lon_in_min = minval_double(n1_in, x1_in);
          lon_in_max = maxval_double(n1_in, x1_in);
          lon_in_avg = avgval_double(n1_in, x1_in);
#pragma acc loop independent reduction(+:nxgrid)
          for(ij=0; ij<=nx2*ny2; ij++) {
            int n_out, i2, j2, n2_in;
            double xarea, dx, lon_out_min, lon_out_max;
            double x2_in[MAX_V], y2_in[MAX_V];

            i2 = ij%nx2;
            j2 = ij/nx2;

            if(lat_out_min_list[ij] >= lat_in_max || lat_out_max_list[ij] <= lat_in_min ) continue;
            /* adjust x2_in according to lon_in_avg*/
            n2_in = n2_list[ij];
#pragma acc loop seq
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
#pragma acc loop seq
              for (l=0; l<n2_in; l++) x2_in[l] += TPI;
            }
            else if (dx >  M_PI) {
              lon_out_min -= TPI;
              lon_out_max -= TPI;
#pragma acc loop seq
              for (l=0; l<n2_in; l++) x2_in[l] -= TPI;
            }

            /* x2_in should in the same range as x1_in after lon_fix, so no need to
               consider cyclic condition
            */
            if(lon_out_min >= lon_in_max || lon_out_max <= lon_in_min ) continue;
            n_out = 1;
            if (  (n_out = clip_2dx2d( x1_in, y1_in, n1_in, x2_in, y2_in, n2_in, x_out, y_out )) > 0) {
              double min_area;
              xarea = poly_area (x_out, y_out, n_out ) * mask_in[j1*nx1+i1];
              min_area = min(area_in[j1*nx1+i1], area_out[j2*nx2+i2]);
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
                nxgrid++;
              }
            }
          }
        }
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
