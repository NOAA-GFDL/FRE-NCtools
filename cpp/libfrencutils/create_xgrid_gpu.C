
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
#include <array>
#include <algorithm>

#include "mpp.h"
#include "BBox3D.h"
#include "Polygon.h"
#include "BoxedObj.h"
#include "create_xgrid.h"

//#include "cartesian_product.hpp"

/**
 * README  IMPORTANT
 * Most of the functions in this file are direct copies of those from the
 * create_xgrid.C file. The reason is that the nvc++ compiler using C++ std
 * parallelism requires many functions to be defined in the same file as the
 * function targeted for GPU execution (e.g. create_xgrid_2dx2d_order2_bfwbb
 * near the bottom of this file). The amount of source copying within this file
 * will be improved ASAP.
 */

using std::sin;
using std::vector;
using std::string;
using std::array;

using BBox_t = nct::BBox3D;
using BPair_t = nct::BoxAndId;
using Point_t = nct::Point3D<double>;


#include "constant.h"
#ifndef MAXXGRID
#define MAXXGRID (1e7)
#endif
#define MV 50

#define AREA_RATIO_THRESH (1.e-6)
#define MASK_THRESH (0.5)
#define EPSLN30 (1.0e-30)

constexpr double MU_TOLORENCE{1.e-6};
constexpr double SMALL_VALUE_PA{1.0E-10};

//TODO: place elsewhere
//NOTE: some of these are are negative values:
constexpr double FILL_VALUE_DOUBLE{ -std::numeric_limits<double>::max()};
constexpr int FILL_VALUE_INT{ std::numeric_limits<int>::max() };
constexpr size_t FILL_VALUE_SIZE_T{ std::numeric_limits<size_t>::max() };


void gpu_error(const string& str)
{
  gpu_error(str.c_str());
}

//TODO: best way to abort processing?
void gpu_error(const char * str)
{
  fprintf(stderr, "Error from  GPU: %s\n",  str );
  //exit(1);  //TODO:
}


void latlon2xyz_gpu(const double lat, const double lon,  std::array<double,3> &  v){
  v[0] = RADIUS * cos(lat) * cos(lon );
  v[1] = RADIUS * cos(lat) * sin(lon);
  v[2] = RADIUS * sin(lat);
}


size_t  latlons_outside_ccd_domain_gpu(const unsigned int NV4, const double *yv, double *xv) {
  size_t count {0};
  for (unsigned int i = 0; i < NV4; i++) {
    if (xv[i] == 2 * M_PI) xv[i] = 0;
    if (xv[i] >= 2 * M_PI || xv[i] < 0.) {
      count++;
    }
    if (yv[i] < -M_PI_2 || yv[i] > M_PI_2) {
      count++;
    }
  }
  return count;
}
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
  //Unexpected branch type (/home/mzuniga/nct_search/FRE-NCtools/cpp/libfrencutils/create_xgrid_gpu.C: 428)
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


inline
size_t pt_idx_gpu(const size_t i, const size_t j,  const size_t nx) {
  return ( j * nx + i);
}
std::array<size_t, 4>
get_cell_idxs_ccw_4_gpu(const size_t i, const size_t j, const size_t nx) {
  std::array<size_t, 4> idxs;
  idxs[0] = pt_idx_gpu(i, j, nx); //ll
  idxs[1] = pt_idx_gpu(i + 1, j , nx); //lr
  idxs[2] = pt_idx_gpu(i + 1, j + 1, nx); //ur
  idxs[3] = pt_idx_gpu(i, j + 1, nx);//ul
  return idxs;
}

std::tuple<double, double>
 getPolygonMinMax_gpu(const double lat_m[], const double lon_m[],
                      const array<size_t, 4> &is, bool debugf=false) {
  constexpr unsigned int NV4{4};
  auto tempy = lat_m[is[0]];
  double maxl = tempy;
  double minl = tempy;
  for (auto i = 1; i < NV4; i++) {
    tempy = lat_m[is[i]];
    if (tempy > maxl) maxl = tempy;
    if (tempy < minl) minl = tempy;
  }
  return {minl, maxl};
}


BBox_t getBoxForSphericalPolygon_gpu(const double lat_m[], const double lon_m[],
                                 const array<size_t, 4> &is, bool debugf=false) {
  constexpr unsigned int NV4{4};
  // xlons are the longitudes that define the X=0 and Y=0 planes(see ll2xyz function)
  // where its possible to have an extrema in X or Y when an edge crosses them.
  constexpr  array<double, 7> xlons {0, M_PI_2, M_PI, 3. * M_PI_2, 2. * M_PI};
  std::array<double, 3> pt{}; // a 3D point.
  BBox_t box; //the thing we are calculating.
  bool crosses_equator = false;
  double yv[ NV4 ], xv[ NV4 ];  //the latitudes (yv) and longitudes (xv) of one polygon/cell.
  for (auto i = 0; i< NV4 ; i++) {
    xv[i] = lon_m[is[i]];
    yv[i] = lat_m[is[i]];
  }

  auto occd = latlons_outside_ccd_domain_gpu(4, yv, xv);
  if(occd>0){
    gpu_error("There are lats or lons outside CCD");
  }

  //Expand the bounding box in consideration of the max and mins of the lats and lons:
  const auto [miny_it, maxy_it] = std::minmax_element(std::begin(yv), std::end(yv));
  const double miny = *miny_it;
  const double maxy = *maxy_it;
  const auto [minx_it, maxx_it] = std::minmax_element(std::begin(xv), std::end(xv));
  const double minx = *minx_it;
  const double maxx = *maxx_it;
  latlon2xyz_gpu( miny, minx, pt);
  box.expand( pt );
  latlon2xyz_gpu( maxy, minx, pt);
  box.expand( pt );
  latlon2xyz_gpu( miny, maxx, pt);
  box.expand( pt );
  latlon2xyz_gpu( maxy, maxx, pt);
  box.expand( pt );

  //For case where a pole is inside a polygon. Note that
  // any Z=k plane of the bbox calculated from above
  // (or even from the actual vertices)
  // would be sufficient to determine this:
  latlon2xyz_gpu(M_PI_2, 0.0, pt); //the north pole
  if (BBox_t::contains_zk(box, pt) && (yv[0] >=  0)){
    box.expand(pt);
  }
  latlon2xyz_gpu(-M_PI_2, 0.0, pt); //the south pole
  if (BBox_t::contains_zk(box, pt) && (yv[0] < 0)){
    box.expand(pt);
  }

  //Expand the box for coordinate value extrema when edge crosses the equator:
  for (auto i = 0; i < NV4; i++) {
    auto equator_lat = 0.;
    auto t1 = yv[i];
    auto t2 = yv[(i + 1) % NV4];
    if (std::signbit(t1) != std::signbit(t2 )) {
      crosses_equator = true;
      latlon2xyz_gpu(equator_lat, xv[i], pt);
      box.expand(pt);
      latlon2xyz_gpu(equator_lat,  xv[(i + 1) % NV4], pt);
      box.expand(pt);
    }
  }

  //Expand the box for coordinate value extrema when edge crosses XZ or YZ plane
  for (auto i = 0; i < NV4; i++) {
    auto t1 = xv[i];
    auto t2 = xv[(i + 1) % NV4];
    for (const auto &xang: xlons) {
      if (std::signbit(t1 - xang) != std::signbit(t2 - xang)) {
        latlon2xyz_gpu(yv[i], xang, pt);
        box.expand(pt);
        latlon2xyz_gpu(yv[(i + 1) % NV4], xang, pt);
        box.expand(pt);
        if(crosses_equator) {
          latlon2xyz_gpu(0.0, xang, pt);
          box.expand(pt);
        }
      }
    }
  }
  box.expand_for_doubles();
  return box;
}

// grid suffix synonyms: '2'/out/target
// grid suffix synonyms: '1'/in/source

//TODO: This formula for the maximum number of cells in one grid to any other in another grid is
// merely a heuristic. A better one might involve max and min cell areas. Consider how many
//cells from one grid can overlap the largest cell from another. The factor of 1000 below
//should take care of really wierd cells.
//size_t get_max_cell_nns(const size_t nx1, const size_t nx2, const size_t ny1, const size_t ny2)
//{ return 1000 * std::max(nx1 * ny1, nx2 * ny2) / std::min(nx1 * ny1, nx2 * ny2); }

//Return the max number of near neighbors for a grid. The worst case may be two grids of the same size
// and one shifted in lat and long. This formula is a heuristic , and any reasonable upper bound will do.
size_t
get_max_grid_nns(const size_t nx1, const size_t nx2, const size_t ny1, const size_t ny2)
{ return 3 * (nx1 * ny1 +  nx2 * ny2); }


void print(std::tuple<int,int> t)
{
  const auto& [a, b] = t;
  std::cout << '(' << a << ' ' << b <<  ')' << std::endl;
}

std::tuple<int, int>
index_pair_from_combo( const int idx_pair, const int nxy2){
  int ij1 = idx_pair / nxy2; //First index
  int ij2 = idx_pair % nxy2;
  return {ij1, ij2};
}

/*
 * Function create_xgrid_2dx2d_order2_bfwbb
 * This function uses the std::copy_if function to collect the [ij1 , ij2]  pairs of indexes of
 * those polygons that pass the "is near neighbors" filters. The filters are bounding-box intersections
 * (and (possibly TBD) also lat/lon overlap tests). The sole purpose of using copy_if
 * is to parallelize the brute-force (O(n1 x n2) ) computations. Note that copy_if is merely (with the
 * cartesian_product class or the iota class) generating the possible index pairs and copying (or saving)
 * only those pairs that pass the neighbors tests.
 * NOTE: This function was written to be compiled for running on GPUS, and complies with the
 *  various restrictions nvc++ 23.7 places on C++ std parallelism usage for GPUS.
 * NOTE: these two triplets are grid suffix synonyms used since the original code. Keeping them in
 * mind may help in understanding: target/'2'/out and source'1'/in.
 *
 * ISSUES: There have been some trouble with copy_if using cartesian_product class (which at this
 * time is not yet available in the C++23 lib). Additionally, there was a compiler bug in creating
 * iotas templated to size_t; and this has forced all the following:
 * A) A single iota is used (instead of the cartesian_product). A trick was used to combine
 *    two integers into one (e.g. ij = i * MAXJ + j). When later needed the 'i' and 'j' are recouped from
 *    from 'ij' (using modulo and division operators).
 * B) The iota is templated on int. Because of the trick in A) above, it would have been better to
 *     use size_t, but unfortunately that does not compile with nvc++.
 * C) To make sure the combine ints fit into one int, an outer loop was introduced. So now a
 *   combine int is made up of ij1 and i2.
 *   TODO: Change to the C++23 cartesian_product when avaialble.
 *   NOTE: Recall C++ is row-major order, so column index (J) should be outer index.
 */

void  create_xgrid_2dx2d_order2_bfwbb(const int nlon_in, const int nlat_in, const int nlon_out, const int nlat_out,
                                       const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                                       const double *mask_in, vector<size_t>& i_in, vector<size_t>& j_in,
                                       vector<size_t>& i_out, vector<size_t>& j_out, vector<double>& xgrid_area,
                                       vector<double>& xgrid_clon, vector<double>& xgrid_clat) {
#define MAX_V 8
  const int nx1{ nlon_in}, nx2{nlon_out}, ny1{ nlat_in}, ny2{nlat_out};
  const int nx1p{nx1 + 1};
  const int nx2p{nx2 + 1};

  //These two are sized int because their use in copy_if
  const int nxy1{nx1 * ny1};
  const int nxy2{nx2 * ny2};

  if (((size_t) nxy1)* ((size_t)ny2) >= std::numeric_limits<int>::max()) {
    // Grids are too big - but only for the current 23.7 nvc++ limitations on
    // using copy_if with iotas templated with ints (i.e. cant use size_t).
    //gpu_error("Grids too big nxy1 * ny2 ) >= std::numeric_limits<int>::max())");
  }
  const int max_grid_nns = get_max_grid_nns(nx1, nx2, ny1, ny2);
  std::cout << "max_grid_nns estimated at :" << max_grid_nns << std::endl;

  std::vector<double> area_in(nx1 * ny1);
  std::vector<double> area_out(nx2 * ny2);
  //TODO: note get_grid_area can be parallelized.
  get_grid_area(nlon_in, nlat_in, lon_in, lat_in, area_in.data());
  get_grid_area(nlon_out, nlat_out, lon_out, lat_out, area_out.data());


  //The out or target pairs
  vector<BBox_t> boxes_2(nx2 * ny2);
  auto ids2 = std::views::iota(0);
  std::for_each_n(std::execution::par, ids2.begin(), nx2 * ny2,
                  [=, boxes_2 = boxes_2.data()](int ij) {
                    auto i2 = ij % nx2;
                    auto j2 = ij / nx2;
                    auto ip = get_cell_idxs_ccw_4_gpu(i2, j2, nx2p);
                    boxes_2[ij] = getBoxForSphericalPolygon_gpu(lat_out, lon_out, ip);
                  });
  //The "in" or source pairs
  vector<BBox_t> boxes_1(nx1 * ny1);
  auto ids1 = std::views::iota(0);
  std::for_each_n(std::execution::par, ids1.begin(), nx1 * ny1,
                  [=, boxes_1 = boxes_1.data()](int ij) {
                    auto i1 = ij % nx1;
                    auto j1 = ij / nx1;
                    auto ip = get_cell_idxs_ccw_4_gpu(i1, j1, nx1p);
                    boxes_1[ij] = getBoxForSphericalPolygon_gpu(lat_in, lon_in, ip);
                  });

  std::cout << "BBox array sizes: " << boxes_1.size() << " ; " << boxes_2.size() << std::endl;

  //NOTE: using std::tuple in liu of std::pair may be needed when switch
  // to cartesian product.

  vector<std::pair<int,int>> nn_pairs;
  for (int j2 = 0; j2 < ny2; j2++) {
    vector<int> nn_pairs1(max_grid_nns, FILL_VALUE_INT);
    //An iota that generates all the continuous integers in [0, nxy1 * nx2):
    auto g12_ids = std::views::iota((int) 0, (int) (nxy1 * nx2));
    std::copy_if(std::execution::par,
                 g12_ids.begin(), g12_ids.end(),
                 nn_pairs1.begin(),
                 [=, boxes_1 = boxes_1.data(), boxes_2 = boxes_2.data()](auto ij12) {
                   auto [ij1, i2] = index_pair_from_combo(ij12, nx2);
                   auto ij2 = j2 * nx2 + i2;
                   return (nct::BBox3D::intersect(boxes_1[ij1], boxes_2[ij2]));
                 });

    //TRIM nn_pairs to only hold the real pairs.
    for (auto &ij12: nn_pairs1) {
      if (ij12 == FILL_VALUE_INT) break;
      auto [ij1, i2] = index_pair_from_combo(ij12, nx2);
      int ij2 = j2 * nx2 + i2;
      nn_pairs.push_back({(int)ij1, ij2});
    }
  }

  //NOTE: For reproducibility with baseline, results are ordered by increasing : i1; j1; ij2
  auto iidx_cmp = [](const std::pair<int,int>& a, const std::pair<int,int>& b) -> bool {
    if(a.first == b.first){
      return  a.second < b.second;
    }else{
      return a.first < b.first;
    }
  };
  std::sort(nn_pairs.begin(), nn_pairs.end(), iidx_cmp); //TODO: parallelize ?


  //Given the initial neighbor pairs, calculate the final pair intersections.
  std::vector<double> xarea_v( nn_pairs.size(), FILL_VALUE_DOUBLE);
  std::vector<double> clon_v( nn_pairs.size(), FILL_VALUE_DOUBLE);
  std::vector<double> clat_v( nn_pairs.size(), FILL_VALUE_DOUBLE);
  vector<int>i_in_v( nn_pairs.size(), FILL_VALUE_INT);
  vector<int>j_in_v( nn_pairs.size(), FILL_VALUE_INT);
  vector<int>i_out_v( nn_pairs.size(), FILL_VALUE_INT);
  vector<int>j_out_v( nn_pairs.size(), FILL_VALUE_INT);

  //NOTE: Parallelizing the loop below may not buy much since size of nn_pairs should
  // already be linear in the size of the larger grid.
  auto ids12 = std::views::iota(0);
  std::for_each_n(std::execution::par, ids12.begin(), nn_pairs.size(),
    [=,nn_pairs=nn_pairs.data(), area_in=area_in.data(), area_out=area_out.data(),
    xarea_v = xarea_v.data(), clon_v = clon_v.data(), clat_v = clat_v.data(),
    i_in_v=i_in_v.data(), j_in_v=j_in_v.data(), i_out_v=i_out_v.data(),
    j_out_v=j_out_v.data()](int ij) {
      auto [ij1, ij2] = nn_pairs[ij];
      auto i1 = ij1 % nx1;
      auto j1 = ij1 / nx1;
      auto i2 = ij2 % nx2;
      auto j2 = ij2 / nx2;
      assert(i1 >= 0 && j1 >= 0);
      assert(i2 < nx2);
      assert(j2 < ny2);

      if (mask_in[j1 * nx1 + i1] > MASK_THRESH) {//NOTE: continue not allowed by ncv++ on GPU;

        double x1_in[MAX_V], y1_in[MAX_V], x2_in[MAX_V], y2_in[MAX_V], x_out[MAX_V], y_out[MAX_V];
        auto idxs_1 = get_cell_idxs_ccw_4_gpu(i1, j1, nx1p);
        for (int i = 0; i < 4; i++) {
          x1_in[i] = lon_in[idxs_1[i]];
          y1_in[i] = lat_in[idxs_1[i]];
        }

        auto idxs_2 = get_cell_idxs_ccw_4_gpu(i2, j2, nx2p);
        for (int i = 0; i < 4; i++) {
          x2_in[i] = lon_out[idxs_2[i]];
          y2_in[i] = lat_out[idxs_2[i]];
        }

        //Some polys require "fixing" before calling area (and clipping) function.
        auto n1_in = fix_lon_gpu(x1_in, y1_in, 4, M_PI);
        auto lon_in_avg = avgval_double_gpu(n1_in, x1_in);
        auto n2_in = fix_lon_gpu(x2_in, y2_in, 4, M_PI);
        auto lon_out_avg = avgval_double_gpu(n2_in, x2_in);

        auto dx = lon_out_avg - lon_in_avg;
        if (dx < -M_PI) {
          for (auto l = 0; l < n2_in; l++) x2_in[l] += TPI;
        }else if (dx > M_PI) {
          for (auto l = 0; l < n2_in; l++) x2_in[l] -= TPI;
        }

        //Call the 2D_by_2D clipping algorithm
        auto n_out = clip_2dx2d_gpu(x1_in, y1_in, n1_in, x2_in,
                                    y2_in, n2_in, x_out, y_out);
        if (n_out > 0) {
          auto xarea = poly_area_gpu(x_out, y_out, n_out) * mask_in[j1 * nx1 + i1];
          auto min_area = std::min(area_in[j1 * nx1 + i1], area_out[j2 * nx2 + i2]);
          if (xarea / min_area > AREA_RATIO_THRESH) {
            xarea_v[ij] = xarea;
            clon_v[ij] = poly_ctrlon_gpu(x_out, y_out, n_out, lon_in_avg);
            clat_v[ij] = poly_ctrlat_gpu(x_out, y_out, n_out);
            i_in_v[ij] = i1; //stores the lower corner of poly
            j_in_v[ij] = j1;
            i_out_v[ij] = i2;
            j_out_v[ij] = j2;
          } //else leave in the FILL_VALUE
        }
      }
    }
  );

  //Save final answers (though copy_if can be used again); but don't use more memory than needed.
  auto nxgrid = count_if(xarea_v.begin(), xarea_v.end(), [](double x) { return x > 0; });

  xgrid_area.reserve(nxgrid);  //TODO: is assign(capacity(, FILL_VALUE) needed?
  i_in.reserve(nxgrid);
  j_in.reserve(nxgrid);
  i_out.reserve(nxgrid);
  j_out.reserve(nxgrid);
  xgrid_clon.reserve(nxgrid);
  xgrid_clat.reserve(nxgrid);

  for(int i = 0; i  < nn_pairs.size(); i ++) {
    if (xarea_v[i] > 0.0) {
      xgrid_area.emplace_back(xarea_v[ i ]);
      xgrid_clon.emplace_back(clon_v[i]);
      xgrid_clat.emplace_back(clat_v[i]);
      i_in.emplace_back(i_in_v[i]);
      j_in.emplace_back(j_in_v[i]);
      i_out .emplace_back(i_out_v[i]);
      j_out.emplace_back(j_out_v[i]);
    }
  }

  std::cout <<"create_xgrid_2dx2d_order2_bff2 end; xgrid_area.size= " << xgrid_area.size() <<std::endl;
}
