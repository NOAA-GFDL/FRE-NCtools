/*
20111203 rsz Created.
*/

#include "ezTest.hpp"
#include "ezNcCF.hpp"

std::string filename1(ez::ezTestRunner& runner) {
  return runner.tmpdir + "./tmp1.nc";
}

bool createFile1(ez::ezTestRunner& runner) {
  ez::ezNc nc;
  
  nc.filename = filename1(runner);
  ezassert( nc.createClobber() == 0 );

  int n1=10;
  int n2=20;
  int n3=30;
  int n4=40;
  int dx = nc.createDim("x", n1);
  ezassert(dx != -1);
  int dy = nc.createDim("y", n2);
  ezassert(dy != -1);
  int dz = nc.createDim("z", n3);
  ezassert(dz != -1);
  int dt = nc.createDim("t", n1);
  ezassert(dt != -1);

  int dimids1dx[1] = {dx};
  int dimids1dy[1] = {dy};
  int dimids1dz[1] = {dz};
  int dimids1dt[1] = {dt};

  ezassert( nc.createVar("x1", NC_FLOAT, 1, dimids1dx) != -1 );
  ezassert( nc.createVar("x2", NC_FLOAT, 1, dimids1dx) != -1 );
  ezassert( nc.createVar("x3", NC_FLOAT, 1, dimids1dx) != -1 );
  ezassert( nc.createVar("x4", NC_FLOAT, 1, dimids1dx) != -1 );
  ezassert( nc.createVar("x5", NC_FLOAT, 1, dimids1dx) != -1 );
  ezassert( nc.createVar("x6", NC_FLOAT, 1, dimids1dx) != -1 );
  ezassert( nc.createVar("x7", NC_FLOAT, 1, dimids1dx) != -1 );
  ezassert( nc.createVar("x8", NC_FLOAT, 1, dimids1dx) != -1 );
  ezassert( nc.createVar("x9", NC_FLOAT, 1, dimids1dx) != -1 );
  ezassert( nc.createVar("x10", NC_FLOAT, 1, dimids1dx) != -1 );
  ezassert( nc.createVar("x11", NC_FLOAT, 1, dimids1dx) != -1 );

  ezassert( nc.createVar("y1", NC_FLOAT, 1, dimids1dy) != -1 );
  ezassert( nc.createVar("y2", NC_FLOAT, 1, dimids1dy) != -1 );
  ezassert( nc.createVar("y3", NC_FLOAT, 1, dimids1dy) != -1 );
  ezassert( nc.createVar("y4", NC_FLOAT, 1, dimids1dy) != -1 );
  ezassert( nc.createVar("y5", NC_FLOAT, 1, dimids1dy) != -1 );
  ezassert( nc.createVar("y6", NC_FLOAT, 1, dimids1dy) != -1 );
  ezassert( nc.createVar("y7", NC_FLOAT, 1, dimids1dy) != -1 );
  ezassert( nc.createVar("y8", NC_FLOAT, 1, dimids1dy) != -1 );
  ezassert( nc.createVar("y9", NC_FLOAT, 1, dimids1dy) != -1 );
  ezassert( nc.createVar("y10", NC_FLOAT, 1, dimids1dy) != -1 );
  ezassert( nc.createVar("y11", NC_FLOAT, 1, dimids1dy) != -1 );
  
  ezassert( nc.createVar("z1", NC_FLOAT, 1, dimids1dz) != -1 );
  ezassert( nc.createVar("z2", NC_FLOAT, 1, dimids1dz) != -1 );
  ezassert( nc.createVar("z3", NC_FLOAT, 1, dimids1dz) != -1 );
  ezassert( nc.createVar("z4", NC_FLOAT, 1, dimids1dz) != -1 );
  ezassert( nc.createVar("z5", NC_FLOAT, 1, dimids1dz) != -1 );
  ezassert( nc.createVar("z6", NC_FLOAT, 1, dimids1dz) != -1 );
  ezassert( nc.createVar("z7", NC_FLOAT, 1, dimids1dz) != -1 );
  ezassert( nc.createVar("z8", NC_FLOAT, 1, dimids1dz) != -1 );
  ezassert( nc.createVar("z9", NC_FLOAT, 1, dimids1dz) != -1 );
  ezassert( nc.createVar("z10", NC_FLOAT, 1, dimids1dz) != -1 );
  ezassert( nc.createVar("z11", NC_FLOAT, 1, dimids1dz) != -1 );
  ezassert( nc.createVar("z12", NC_FLOAT, 1, dimids1dz) != -1 );
  ezassert( nc.createVar("z13", NC_FLOAT, 1, dimids1dz) != -1 );
  ezassert( nc.createVar("z14", NC_FLOAT, 1, dimids1dz) != -1 );
  ezassert( nc.createVar("z15", NC_FLOAT, 1, dimids1dz) != -1 );
  ezassert( nc.createVar("z16", NC_FLOAT, 1, dimids1dz) != -1 );
  ezassert( nc.createVar("z17", NC_FLOAT, 1, dimids1dz) != -1 );
  ezassert( nc.createVar("z18", NC_FLOAT, 1, dimids1dz) != -1 );
  ezassert( nc.createVar("z19", NC_FLOAT, 1, dimids1dz) != -1 );
  ezassert( nc.createVar("z20", NC_FLOAT, 1, dimids1dz) != -1 );
  ezassert( nc.createVar("z21", NC_FLOAT, 1, dimids1dz) != -1 );
  ezassert( nc.createVar("z22", NC_FLOAT, 1, dimids1dz) != -1 );
  ezassert( nc.createVar("z23", NC_FLOAT, 1, dimids1dz) != -1 );
  ezassert( nc.createVar("z24", NC_FLOAT, 1, dimids1dz) != -1 );
  ezassert( nc.createVar("z25", NC_FLOAT, 1, dimids1dz) != -1 );
  ezassert( nc.createVar("z26", NC_FLOAT, 1, dimids1dz) != -1 );
  ezassert( nc.createVar("z27", NC_FLOAT, 1, dimids1dz) != -1 );
  ezassert( nc.createVar("z28", NC_FLOAT, 1, dimids1dz) != -1 );
  ezassert( nc.createVar("z29", NC_FLOAT, 1, dimids1dz) != -1 );
  ezassert( nc.createVar("z30", NC_FLOAT, 1, dimids1dz) != -1 );
  ezassert( nc.createVar("z31", NC_FLOAT, 1, dimids1dz) != -1 );
  ezassert( nc.createVar("z32", NC_FLOAT, 1, dimids1dz) != -1 );
  ezassert( nc.createVar("z33", NC_FLOAT, 1, dimids1dz) != -1 );
  ezassert( nc.createVar("z34", NC_FLOAT, 1, dimids1dz) != -1 );
  ezassert( nc.createVar("z35", NC_FLOAT, 1, dimids1dz) != -1 );
  ezassert( nc.createVar("z36", NC_FLOAT, 1, dimids1dz) != -1 );
  ezassert( nc.createVar("z37", NC_FLOAT, 1, dimids1dz) != -1 );
  ezassert( nc.createVar("z38", NC_FLOAT, 1, dimids1dz) != -1 );
  ezassert( nc.createVar("z39", NC_FLOAT, 1, dimids1dz) != -1 );

  ezassert( nc.createVar("t1", NC_FLOAT, 1, dimids1dt) != -1 );
  ezassert( nc.createVar("t2", NC_FLOAT, 1, dimids1dt) != -1 );
  ezassert( nc.createVar("t3", NC_FLOAT, 1, dimids1dt) != -1 );
  ezassert( nc.createVar("t4", NC_FLOAT, 1, dimids1dt) != -1 );
  ezassert( nc.createVar("t5", NC_FLOAT, 1, dimids1dt) != -1 );
  ezassert( nc.createVar("t6", NC_FLOAT, 1, dimids1dt) != -1 );
  ezassert( nc.createVar("t7", NC_FLOAT, 1, dimids1dt) != -1 );
  ezassert( nc.createVar("t8", NC_FLOAT, 1, dimids1dt) != -1 );
  ezassert( nc.createVar("t9", NC_FLOAT, 1, dimids1dt) != -1 );
  ezassert( nc.createVar("t10", NC_FLOAT, 1, dimids1dt) != -1 );
  ezassert( nc.createVar("t11", NC_FLOAT, 1, dimids1dt) != -1 );
  ezassert( nc.createVar("t12", NC_FLOAT, 1, dimids1dt) != -1 );
  ezassert( nc.createVar("t13", NC_FLOAT, 1, dimids1dt) != -1 );
  ezassert( nc.createVar("t14", NC_FLOAT, 1, dimids1dt) != -1 );
  ezassert( nc.createVar("t15", NC_FLOAT, 1, dimids1dt) != -1 );
  ezassert( nc.createVar("t16", NC_FLOAT, 1, dimids1dt) != -1 );
  ezassert( nc.createVar("t17", NC_FLOAT, 1, dimids1dt) != -1 );
  ezassert( nc.createVar("t18", NC_FLOAT, 1, dimids1dt) != -1 );
  ezassert( nc.createVar("t19", NC_FLOAT, 1, dimids1dt) != -1 );
  ezassert( nc.createVar("t20", NC_FLOAT, 1, dimids1dt) != -1 );
  ezassert( nc.createVar("t21", NC_FLOAT, 1, dimids1dt) != -1 );
  ezassert( nc.createVar("t22", NC_FLOAT, 1, dimids1dt) != -1 );
  ezassert( nc.createVar("t23", NC_FLOAT, 1, dimids1dt) != -1 );
  ezassert( nc.createVar("t24", NC_FLOAT, 1, dimids1dt) != -1 );
  ezassert( nc.createVar("t25", NC_FLOAT, 1, dimids1dt) != -1 );
  
  std::vector<std::string> strings;
  strings.resize(1);
  strings[0] = "longitude";
  ezassert( nc.putAtt("x1", "standard_name", NC_CHAR, strings, 0, 0) == 0);
  ezassert( nc.putAtt("x2", "long_name", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "x";
  ezassert( nc.putAtt("x3", "axis", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "degrees_east";
  ezassert( nc.putAtt("x4", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "degree_east";
  ezassert( nc.putAtt("x5", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "degree_E";
  ezassert( nc.putAtt("x6", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "degrees_E";
  ezassert( nc.putAtt("x7", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "degreeE";
  ezassert( nc.putAtt("x8", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "degreesE";
  ezassert( nc.putAtt("x9", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "grid_longitude";
  ezassert( nc.putAtt("x10", "standard_name", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "DEGREES_EAST";
  ezassert( nc.putAtt("x11", "units", NC_CHAR, strings, 0, 0) == 0);

  strings[0] = "latitude";
  ezassert( nc.putAtt("y1", "standard_name", NC_CHAR, strings, 0, 0) == 0);
  ezassert( nc.putAtt("y2", "long_name", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "y";
  ezassert( nc.putAtt("y3", "axis", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "degrees_north";
  ezassert( nc.putAtt("y4", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "degree_north";
  ezassert( nc.putAtt("y5", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "degree_N";
  ezassert( nc.putAtt("y6", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "degrees_N";
  ezassert( nc.putAtt("y7", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "degreeN";
  ezassert( nc.putAtt("y8", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "degreesN";
  ezassert( nc.putAtt("y9", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "grid_latitude";
  ezassert( nc.putAtt("y10", "standard_name", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "DEGREES_NORTH";
  ezassert( nc.putAtt("y11", "units", NC_CHAR, strings, 0, 0) == 0);

  strings[0] = "atmosphere_sigma_coordinate";
  ezassert( nc.putAtt("z1", "standard_name", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "atmosphere_hybrid_sigma_pressure_coordinate";
  ezassert( nc.putAtt("z2", "standard_name", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "atmosphere_hybrid_height_coordinate";
  ezassert( nc.putAtt("z3", "standard_name", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "atmosphere_sleve_coordinate";
  ezassert( nc.putAtt("z4", "standard_name", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "ocean_sigma_coordinate";
  ezassert( nc.putAtt("z5", "standard_name", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "ocean_s_coordinate";
  ezassert( nc.putAtt("z6", "standard_name", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "ocean_sigma_z_coordinate";
  ezassert( nc.putAtt("z7", "standard_name", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "ocean_double_sigma_coordinate";
  ezassert( nc.putAtt("z8", "standard_name", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "elevation";
  ezassert( nc.putAtt("z9", "long_name", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "z";
  ezassert( nc.putAtt("z10", "axis", NC_CHAR, strings, 0, 0) == 0);
  
  strings[0] = "bar";
  ezassert( nc.putAtt("z11", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "millibar";
  ezassert( nc.putAtt("z12", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "decibar";
  ezassert( nc.putAtt("z13", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "atmosphere";
  ezassert( nc.putAtt("z14", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "atm";
  ezassert( nc.putAtt("z15", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "pascal";
  ezassert( nc.putAtt("z16", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "Pa";
  ezassert( nc.putAtt("z17", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "hPa";
  ezassert( nc.putAtt("z18", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "meter";
  ezassert( nc.putAtt("z19", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "m";
  ezassert( nc.putAtt("z20", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "metre";
  ezassert( nc.putAtt("z21", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "kilometer";
  ezassert( nc.putAtt("z22", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "km";
  ezassert( nc.putAtt("z23", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "level";
  ezassert( nc.putAtt("z24", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "layer";
  ezassert( nc.putAtt("z25", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "sigma_level";
  ezassert( nc.putAtt("z26", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "bars";
  ezassert( nc.putAtt("z27", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "millibars";
  ezassert( nc.putAtt("z28", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "decibars";
  ezassert( nc.putAtt("z29", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "atmospheres";
  ezassert( nc.putAtt("z30", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "pascals";
  ezassert( nc.putAtt("z31", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "meters";
  ezassert( nc.putAtt("z32", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "metres";
  ezassert( nc.putAtt("z33", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "kilometers";
  ezassert( nc.putAtt("z34", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "levels";
  ezassert( nc.putAtt("z35", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "layers";
  ezassert( nc.putAtt("z36", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "sigma_levels";
  ezassert( nc.putAtt("z37", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "up";
  ezassert( nc.putAtt("z38", "positive", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "down";
  ezassert( nc.putAtt("z39", "positive", NC_CHAR, strings, 0, 0) == 0);

  strings[0] = "time";
  ezassert( nc.putAtt("t1", "standard_name", NC_CHAR, strings, 0, 0) == 0);
  ezassert( nc.putAtt("t2", "long_name", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "T";
  ezassert( nc.putAtt("t3", "axis", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "day since 1980-0-0";
  ezassert( nc.putAtt("t4", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "d since 1980-0-0";
  ezassert( nc.putAtt("t5", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "hour since 1980-0-0";
  ezassert( nc.putAtt("t6", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "hr since 1980-0-0";
  ezassert( nc.putAtt("t7", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "h since 1980-0-0";
  ezassert( nc.putAtt("t8", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "minute since 1980-0-0";
  ezassert( nc.putAtt("t9", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "min since 1980-0-0";
  ezassert( nc.putAtt("t10", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "second since 1980-0-0";
  ezassert( nc.putAtt("t11", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "sec since 1980-0-0";
  ezassert( nc.putAtt("t12", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "s since 1980-0-0";
  ezassert( nc.putAtt("t13", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "year since 1980-0-0";
  ezassert( nc.putAtt("t14", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "month since 1980-0-0";
  ezassert( nc.putAtt("t15", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "days since 1980-0-0";
  ezassert( nc.putAtt("t16", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "hours since 1980-0-0";
  ezassert( nc.putAtt("t17", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "minutes since 1980-0-0";
  ezassert( nc.putAtt("t18", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "seconds since 1980-0-0";
  ezassert( nc.putAtt("t19", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "years since 1980-0-0";
  ezassert( nc.putAtt("t20", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "months since 1980-0-0";
  ezassert( nc.putAtt("t21", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "31 28";
  ezassert( nc.putAtt("t22", "month_lengths", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "none";
  ezassert( nc.putAtt("t23", "calendar", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "1900";
  ezassert( nc.putAtt("t24", "leap_year", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "2";
  ezassert( nc.putAtt("t25", "leap_month", NC_CHAR, strings, 0, 0) == 0);

  ezassert( nc.close() == 0 );

  return true;
}
  
bool test_x(ez::ezTestRunner& runner) {
  ezassert( createFile1(runner) );
  ez::ezNc nc;
  nc.filename = filename1(runner);

  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );

  ezassert( ez::ezNcCF::isX(nc.getVar("x1")) );
  ezassert( ez::ezNcCF::isX(nc.getVar("x2")) );
  ezassert( ez::ezNcCF::isX(nc.getVar("x3")) );
  ezassert( ez::ezNcCF::isX(nc.getVar("x4")) );
  ezassert( ez::ezNcCF::isX(nc.getVar("x5")) );
  ezassert( ez::ezNcCF::isX(nc.getVar("x6")) );
  ezassert( ez::ezNcCF::isX(nc.getVar("x7")) );
  ezassert( ez::ezNcCF::isX(nc.getVar("x8")) );
  ezassert( ez::ezNcCF::isX(nc.getVar("x9")) );
  ezassert( ez::ezNcCF::isX(nc.getVar("x10")) );
  ezassert( ez::ezNcCF::isX(nc.getVar("x11")) );
  ezassert( ez::ezNcCF::isX(nc.getVar("y1")) == false );
  
  ezassert( nc.close() == 0);

  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );

  return true;
}

bool test_y(ez::ezTestRunner& runner) {
  ezassert( createFile1(runner) );
  ez::ezNc nc;
  nc.filename = filename1(runner);

  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );

  ezassert( ez::ezNcCF::isY(nc.getVar("y1")) );
  ezassert( ez::ezNcCF::isY(nc.getVar("y2")) );
  ezassert( ez::ezNcCF::isY(nc.getVar("y3")) );
  ezassert( ez::ezNcCF::isY(nc.getVar("y4")) );
  ezassert( ez::ezNcCF::isY(nc.getVar("y5")) );
  ezassert( ez::ezNcCF::isY(nc.getVar("y6")) );
  ezassert( ez::ezNcCF::isY(nc.getVar("y7")) );
  ezassert( ez::ezNcCF::isY(nc.getVar("y8")) );
  ezassert( ez::ezNcCF::isY(nc.getVar("y9")) );
  ezassert( ez::ezNcCF::isY(nc.getVar("y10")) );
  ezassert( ez::ezNcCF::isY(nc.getVar("y11")) );
  ezassert( ez::ezNcCF::isY(nc.getVar("x1")) == false );
  
  ezassert( nc.close() == 0);

  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );

  return true;
}

bool test_z(ez::ezTestRunner& runner) {
  ezassert( createFile1(runner) );
  ez::ezNc nc;
  nc.filename = filename1(runner);

  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );

  ezassert( ez::ezNcCF::isZ(nc.getVar("z1")) );
  ezassert( ez::ezNcCF::isZ(nc.getVar("z2")) );
  ezassert( ez::ezNcCF::isZ(nc.getVar("z3")) );
  ezassert( ez::ezNcCF::isZ(nc.getVar("z4")) );
  ezassert( ez::ezNcCF::isZ(nc.getVar("z5")) );
  ezassert( ez::ezNcCF::isZ(nc.getVar("z6")) );
  ezassert( ez::ezNcCF::isZ(nc.getVar("z7")) );
  ezassert( ez::ezNcCF::isZ(nc.getVar("z8")) );
  ezassert( ez::ezNcCF::isZ(nc.getVar("z9")) );
  ezassert( ez::ezNcCF::isZ(nc.getVar("z10")) );
  ezassert( ez::ezNcCF::isZ(nc.getVar("z11")) );
  ezassert( ez::ezNcCF::isZ(nc.getVar("z12")) );
  ezassert( ez::ezNcCF::isZ(nc.getVar("z13")) );
  ezassert( ez::ezNcCF::isZ(nc.getVar("z14")) );
  ezassert( ez::ezNcCF::isZ(nc.getVar("z15")) );
  ezassert( ez::ezNcCF::isZ(nc.getVar("z16")) );
  ezassert( ez::ezNcCF::isZ(nc.getVar("z17")) );
  ezassert( ez::ezNcCF::isZ(nc.getVar("z18")) );
  ezassert( ez::ezNcCF::isZ(nc.getVar("z19")) );
  ezassert( ez::ezNcCF::isZ(nc.getVar("z20")) );
  ezassert( ez::ezNcCF::isZ(nc.getVar("z21")) );
  ezassert( ez::ezNcCF::isZ(nc.getVar("z22")) );
  ezassert( ez::ezNcCF::isZ(nc.getVar("z23")) );
  ezassert( ez::ezNcCF::isZ(nc.getVar("z24")) );
  ezassert( ez::ezNcCF::isZ(nc.getVar("z25")) );
  ezassert( ez::ezNcCF::isZ(nc.getVar("z26")) );
  ezassert( ez::ezNcCF::isZ(nc.getVar("z27")) );
  ezassert( ez::ezNcCF::isZ(nc.getVar("z28")) );
  ezassert( ez::ezNcCF::isZ(nc.getVar("z29")) );
  ezassert( ez::ezNcCF::isZ(nc.getVar("z30")) );
  ezassert( ez::ezNcCF::isZ(nc.getVar("z31")) );
  ezassert( ez::ezNcCF::isZ(nc.getVar("z32")) );
  ezassert( ez::ezNcCF::isZ(nc.getVar("z33")) );
  ezassert( ez::ezNcCF::isZ(nc.getVar("z34")) );
  ezassert( ez::ezNcCF::isZ(nc.getVar("z35")) );
  ezassert( ez::ezNcCF::isZ(nc.getVar("z36")) );
  ezassert( ez::ezNcCF::isZ(nc.getVar("z37")) );
  ezassert( ez::ezNcCF::isZ(nc.getVar("z38")) );
  ezassert( ez::ezNcCF::isZ(nc.getVar("z39")) );

  ezassert( ez::ezNcCF::isZ(nc.getVar("x1")) == false );
  ezassert( ez::ezNcCF::isZ(nc.getVar("x8")) == false );
  ezassert( ez::ezNcCF::isZ(nc.getVar("y4")) == false );
  
  ezassert( nc.close() == 0);

  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );

  return true;
}

bool test_t(ez::ezTestRunner& runner) {
  ezassert( createFile1(runner) );
  ez::ezNc nc;
  nc.filename = filename1(runner);

  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );

  ezassert( ez::ezNcCF::isT(nc.getVar("t1")) );
  ezassert( ez::ezNcCF::isT(nc.getVar("t2")) );
  ezassert( ez::ezNcCF::isT(nc.getVar("t3")) );
  ezassert( ez::ezNcCF::isT(nc.getVar("t4")) );
  ezassert( ez::ezNcCF::isT(nc.getVar("t5")) );
  ezassert( ez::ezNcCF::isT(nc.getVar("t6")) );
  ezassert( ez::ezNcCF::isT(nc.getVar("t7")) );
  ezassert( ez::ezNcCF::isT(nc.getVar("t8")) );
  ezassert( ez::ezNcCF::isT(nc.getVar("t9")) );
  ezassert( ez::ezNcCF::isT(nc.getVar("t10")) );
  ezassert( ez::ezNcCF::isT(nc.getVar("t11")) );
  ezassert( ez::ezNcCF::isT(nc.getVar("t12")) );
  ezassert( ez::ezNcCF::isT(nc.getVar("t13")) );
  ezassert( ez::ezNcCF::isT(nc.getVar("t14")) );
  ezassert( ez::ezNcCF::isT(nc.getVar("t15")) );
  ezassert( ez::ezNcCF::isT(nc.getVar("t16")) );
  ezassert( ez::ezNcCF::isT(nc.getVar("t17")) );
  ezassert( ez::ezNcCF::isT(nc.getVar("t18")) );
  ezassert( ez::ezNcCF::isT(nc.getVar("t19")) );
  ezassert( ez::ezNcCF::isT(nc.getVar("t20")) );
  ezassert( ez::ezNcCF::isT(nc.getVar("t21")) );
  ezassert( ez::ezNcCF::isT(nc.getVar("t22")) );
  ezassert( ez::ezNcCF::isT(nc.getVar("t23")) );
  ezassert( ez::ezNcCF::isT(nc.getVar("t24")) );
  ezassert( ez::ezNcCF::isT(nc.getVar("t25")) );

  ezassert( ez::ezNcCF::isT(nc.getVar("x1")) == false );
  ezassert( ez::ezNcCF::isT(nc.getVar("y8")) == false );
  ezassert( ez::ezNcCF::isT(nc.getVar("z4")) == false );
  
  ezassert( nc.close() == 0);

  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );

  return true;
}

int main(int argc, const char* argv[]) {
  ez::ezTestRunner runner;
  runner.getopt(argc, argv);
  #define TEST(N) runner.add(&test_##N, "test_"#N);
  TEST(x);
  TEST(y);
  TEST(z);
  TEST(t);
  
  return runner.run();
}