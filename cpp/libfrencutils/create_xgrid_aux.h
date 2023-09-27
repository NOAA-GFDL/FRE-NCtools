
#ifndef FREGRID_CREATE_XGRID_AUX_H
#define FREGRID_CREATE_XGRID_AUX_H

#include <cmath>
#include <vector>
#include <array>
#include <algorithm>
#include <span>
#include <source_location>

#include "constant.h"
#include "mpp.h"
#include "create_xgrid.h"


#include "BBox3D.h"
#include "BoxedObj.h"
#include "Polygon.h"
#include "mosaic_util.h"


static bool checkBBoxViaPolySamples(std::span<double> yv, std::span<double> xv,
                              nct::BBox3D & box, unsigned int npoints1D = 5 ) {
  bool passed = true;
  std::array<double, 3> pt{};
  const auto [miny_it, maxy_it] = std::minmax_element(begin(yv), end(yv));
  const double miny = *miny_it;
  const double maxy = *maxy_it;
  const auto [minx_it, maxx_it] = std::minmax_element(begin(xv), end(xv));
  const double minx = *minx_it;
  const double maxx = *maxx_it;

  const unsigned int NDX{npoints1D};
  const unsigned int NDY{npoints1D};
  const double dy = (maxy - miny) / NDY;
  const double dx = (maxx - minx) / NDX;
  for (int i = 0; i < NDX; ++i) {
    auto ax = minx + i * dx;
    for (int j = 0; j < NDY; ++j) {
      latlon2xyz(miny + j * dy, ax, pt);
      auto contains_point = nct::BBox3D::contains(box, pt);
      if (!contains_point) {
        passed = false;
        //TODO: DEBUG or CTest option ? : saved to strstream the polygon and the box ?
       //  printPolygon<double>(std::cout, xv, yv);
       //  std::cout << box << std::endl;
      //   auto str = std::format("< {:16.10e}, {:16.10e}, {:16.10e}>", pt[0],pt[1],pt[2]);
      //   std::cout << str <<std::endl;
      }
    }
  }
  return passed;
}
#endif //FREGRID_CREATE_XGRID_AUX_H
