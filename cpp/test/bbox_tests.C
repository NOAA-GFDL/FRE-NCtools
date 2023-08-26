#include <array>
#include <sstream>
#include <string>
// #include <fmt/format.h>
//#include <format>

#include <catch2/catch_test_macros.hpp>
//#include "ApprovalTests.hpp"

#include "constant.h"
#include "create_xgrid.h"
#include "create_xgrid_aux.h"

TEST_CASE("TEST SIMPLE_POLY BOX")
{
  std::array<size_t,4> is {0,1,2,3};
  //std::vector<double> lats{88.98,88.98, 88.72, 88.72 };
 // std::vector<double> lons{ 251.7, 148.3,57.81, 342.2 };
  std::vector<double> lats{0,1, 0, 1 };
  std::vector<double> lons{ 44, 45, 46, 45};
  for(int i = 0; i<4; i++){
    lats[i] = lats[i] * D2R;
    lons[i] = lons[i] * D2R;
  }
  auto box = getBoxForSphericalPolygon(lats.data(), lons.data(), is);
  bool passed = checkBBoxViaPolySamples({lats.data(),4}, {lons.data(),4}, box, 5);

  REQUIRE( passed == true );
}


TEST_CASE("TEST NIKI_POLY BOX")
{
  std::array<size_t,4> is {0,1,2,3};
  std::vector<double> lats{88.98,88.98, 88.72, 88.72 };
  std::vector<double> lons{ 251.7, 148.3,57.81, 342.2 };
  for(int i = 0; i<4; i++){
    lats[i] = lats[i] * D2R;
    lons[i] = lons[i] * D2R;
  }
  auto box = getBoxForSphericalPolygon(lats.data(), lons.data(), is);
  bool passed = checkBBoxViaPolySamples({lats.data(),4}, {lons.data(),4}, box, 5);

  REQUIRE( passed == true );
}

//TODO: Many More Tests