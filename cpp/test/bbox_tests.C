#include <array>
#include <sstream>
#include <string>
// #include <fmt/format.h>
//#include <format>

#include <catch2/catch_test_macros.hpp>
//#include "ApprovalTests.hpp"

#include "constant.h"
#include "BBox3D.h"
#include "create_xgrid.h"
#include "create_xgrid_aux.h"

TEST_CASE("TEST SIMPLE_SPOLY BBOX CREATE")
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


TEST_CASE("TEST NIKI_SPOLY BBOX CREATE")
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


TEST_CASE("TEST GCA SPOLY1 BBOX CREATE")
{
std::array<size_t,4> is {0,1,2,3};
std::vector<double> lats{-0.7696676255242643,-0.7679656084524836, -0.7682291223935688, -0.7699311689891934 };
std::vector<double> lons{ 5.87284438368007, 5.87284438368007,5.875049806476964, 5.875049806476964 };
auto box = getBoxForSphericalPolygon(lats.data(), lons.data(), is);
bool passed = checkBBoxViaPolySamples({lats.data(),4}, {lons.data(),4}, box, 500);

REQUIRE( passed == true );
}

TEST_CASE("TEST GCA SPOLY #2 BBOX CREATE")
{
std::array<size_t,4> is {0,1,2,3};
std::vector<double> lats{-0.767944870877505,-0.7330382858376184, -0.7330382858376184, -0.767944870877505};
std::vector<double> lons{ 5.061454830783556,  5.026548245743669,5.061454830783556, 5.096361415823442 };
auto box = getBoxForSphericalPolygon(lats.data(), lons.data(), is);
bool passed = checkBBoxViaPolySamples({lats.data(),4}, {lons.data(),4}, box, 500);

REQUIRE( passed == true );
}


TEST_CASE("TEST GCA SPOLY 1-2 BBOX INTERSECT")
{
  std::array<size_t,4> is {0,1,2,3};
  std::vector<double> lats1{-0.7696676255242643 ,-0.7679656084524836 ,-0.7682291223935688 ,-0.7699311689891934 };
  std::vector<double> lons1{ 5.87284438368007 , 5.87284438368007, 5.875049806476964,5.875049806476964  };

  std::vector<double> lats2{-0.767944870877505,-0.7330382858376184, -0.7330382858376184, -0.767944870877505};
  std::vector<double> lons2{ 5.061454830783556,  5.026548245743669,5.061454830783556, 5.096361415823442 };

  auto box1 = getBoxForSphericalPolygon(lats1.data(), lons1.data(), is);
  auto box2 = getBoxForSphericalPolygon(lats2.data(), lons2.data(), is);

 auto isect = nct::BBox3D::intersect(box1, box2);

  REQUIRE( isect == false );
}
//TODO: Many More Tests