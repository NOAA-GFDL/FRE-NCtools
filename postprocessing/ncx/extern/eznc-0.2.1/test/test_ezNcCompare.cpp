/*
20111125 rsz Created.
*/

#include "ezTest.hpp"
#include "ezNcCompare.hpp"

std::string filename1(ez::ezTestRunner& runner) {
  return runner.tmpdir + "./tmp1.nc";
}

std::string filename2(ez::ezTestRunner& runner) {
  return runner.tmpdir + "./tmp2.nc";
}

bool createFile1(ez::ezTestRunner& runner) {
  ez::ezNc nc;
  
  nc.filename = filename1(runner);
  ezassert( nc.createClobber() == 0 );

  int n1=10;
  int n2=20;
  int n3=30;
  int n4=40;
  int d1 = nc.createDim("d1", n1);
  ezassert(d1 != -1);
  int d2 = nc.createDim("d2", n2);
  ezassert(d2 != -1);
  int d3 = nc.createDim("d3", n3);
  ezassert(d3 != -1);
  int d4 = nc.createDim("d4", n4);
  ezassert(d4 != -1);
  int dr = nc.createDim("dr", 2, true);
  ezassert(dr != -1);

  int dimids1d[1] = {d1};
  int dimids2d[2] = {d1,d2};
  int dimids3d[3] = {d1,d2,d3};
  int dimids0dt[1] = {dr};
  int dimids1dt[2] = {dr,d1};
  int dimids2dt[3] = {dr,d1,d2};
  int dimids3dt[4] = {dr,d1,d2,d3};
  int v1 = nc.createVar("v1", NC_BYTE, 1, dimids1d);
  ezassert(v1 != -1);
  int v2 = nc.createVar("v2", NC_CHAR, 2, dimids2d);
  ezassert(v2 != -1);
  int v3 = nc.createVar("v3", NC_SHORT, 3, dimids3d);
  ezassert(v3 != -1);
  int v4 = nc.createVar("v4", NC_INT, 1, dimids0dt);
  ezassert(v4 != -1);
  int v5 = nc.createVar("v5", NC_FLOAT, 2, dimids1dt);
  ezassert(v5 != -1);
  int v6 = nc.createVar("v6", NC_DOUBLE, 3, dimids2dt);
  ezassert(v6 != -1);
  int v7 = nc.createVar("v7", NC_DOUBLE, 4, dimids3dt);
  ezassert(v7 != -1);
  
  std::vector<std::string> strings;
  strings.resize(3);
  strings[0] = "Quick Fox";
  ezassert( nc.putAtt(0, "Author", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "2123000";
  ezassert( nc.putAtt(0, "age", NC_INT, strings, 0, 0) == 0);
  strings[0] = "1e-6";
  ezassert( nc.putAtt(0, "microsecond", NC_FLOAT, strings, 0, 0) == 0);
  strings[0] = "1e-9";
  ezassert( nc.putAtt(0, "nanosecond", NC_DOUBLE, strings, 0, 0) == 0);

  strings[0] = "32000";
  ezassert( nc.putAtt("v3", "offset", NC_SHORT, strings, 0, 0) == 0);
  strings[0] = "days since 1980-01-01 12:30:45";
  ezassert( nc.putAtt("v4", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "-2123000";
  ezassert( nc.putAtt("v4", "missing_value", NC_INT, strings, 0, 0) == 0);
  strings[0] = "3.14";
  ezassert( nc.putAtt("v4", "pi", NC_FLOAT, strings, 0, 0) == 0);
  strings[0] = "1.1";
  strings[1] = "2.2";
  strings[2] = "3.3";
  ezassert( nc.putAtt("v4", "double_array", NC_DOUBLE, strings, 0, 2) == 0);
  strings[0] = "-1e30";
  ezassert( nc.putAtt("v5", "_FillValue", NC_FLOAT, strings, 0, 0) == 0);
  strings[0] = "-1e30";
  strings[1] = "1e30";
  ezassert( nc.putAtt("v6", "valid_range", NC_DOUBLE, strings, 0, 1) == 0);
  strings[0] = "quick brown fox jumped";
  ezassert( nc.putAtt("v7", "string", NC_CHAR, strings, 0, 0) == 0);

  ezassert( nc.enddef() == 0 ); 

  ez::ezNcBuffers buf;
  ezassert( nc.reserve(buf) == 0 );
  buf.allocate();

  // Create dummy data.
  ez::seq(buf.uc, buf.nuc);
  ez::seq(buf.c, buf.nc);
  ez::seq(buf.s, buf.ns);
  ez::seq(buf.i, buf.ni);
  ez::seq(buf.f, buf.nf);
  ez::seq(buf.d, buf.nd);
  
  size_t start[4], count[4];
  // Static vars.
  start[0] = start[1] = start[2] = start[3] = 0;
  count[0] = n1; count[1] = n2; count[2] = n3;
  ezassert( nc.write(v1, &buf, start, count) == 0 );
  ezassert( nc.write(v2, &buf, start, count) == 0 );
  ezassert( nc.write(v3, &buf, start, count) == 0 );

  count[0] = 1; count[1] = n1; count[2] = n2; count[3] = n3;
  ezassert( nc.write(v4, &buf, start, count) == 0 );
  ezassert( nc.write(v5, &buf, start, count) == 0 );
  ezassert( nc.write(v6, &buf, start, count) == 0 );
  ezassert( nc.write(v7, &buf, start, count) == 0 );

  start[0] = 1;
  ezassert( nc.write(v4, &buf, start, count) == 0 );
  ezassert( nc.write(v5, &buf, start, count) == 0 );
  ezassert( nc.write(v6, &buf, start, count) == 0 );
  ezassert( nc.write(v7, &buf, start, count) == 0 );
  
  ezassert( nc.close() == 0 );

  return true;
}

bool createFile2(ez::ezTestRunner& runner) {
  ez::ezNc nc;
  
  nc.filename = filename2(runner);
  ezassert( nc.create64Clobber() == 0 );

  int n1=10;
  int n2=20;
  int n3=30;
  int n4=40;
  int d1 = nc.createDim("d1", n1);
  ezassert(d1 != -1);
  int d2 = nc.createDim("d2", n2);
  ezassert(d2 != -1);
  int d3 = nc.createDim("d3", n3);
  ezassert(d3 != -1);
  int dr = nc.createDim("dr", 2, true);
  ezassert(dr != -1);
  
  int dimids1d[1] = {d1};
  int dimids2d[2] = {d1,d2};
  int dimids3d[3] = {d1,d2,d3};
  int dimids0dt[1] = {dr};
  int dimids1dt[2] = {dr,d1};
  int dimids2dt[3] = {dr,d1,d2};
  int dimids3dt[4] = {dr,d1,d2,d3};
  int v1 = nc.createVar("v1", NC_BYTE, 1, dimids1d);
  ezassert(v1 != -1);
  int v2 = nc.createVar("v2", NC_CHAR, 2, dimids2d);
  ezassert(v2 != -1);
  int v3 = nc.createVar("v3", NC_SHORT, 3, dimids3d);
  ezassert(v3 != -1);
  int v4 = nc.createVar("v4", NC_INT, 1, dimids0dt);
  ezassert(v4 != -1);
  int v5 = nc.createVar("v5", NC_FLOAT, 2, dimids1dt);
  ezassert(v5 != -1);
  int v6 = nc.createVar("v6", NC_DOUBLE, 3, dimids2dt);
  ezassert(v6 != -1);
  int v7 = nc.createVar("v7", NC_DOUBLE, 4, dimids3dt);
  ezassert(v7 != -1);
  
  std::vector<std::string> strings;
  strings.resize(3);
  strings[0] = "Quick Fox 64";
  ezassert( nc.putAtt(0, "Author", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "2123000";
  ezassert( nc.putAtt(0, "age", NC_INT, strings, 0, 0) == 0);
  strings[0] = "1e-6";
  ezassert( nc.putAtt(0, "microsecond_64", NC_FLOAT, strings, 0, 0) == 0);
  strings[0] = "1e-9";
  ezassert( nc.putAtt(0, "nanosecond", NC_DOUBLE, strings, 0, 0) == 0);

  strings[0] = "32000";
  ezassert( nc.putAtt("v3", "offset", NC_SHORT, strings, 0, 0) == 0);
  strings[0] = "days since 1980-01-01 12:30:45";
  ezassert( nc.putAtt("v4", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "-2123000";
  ezassert( nc.putAtt("v4", "missing_value", NC_INT, strings, 0, 0) == 0);
  strings[0] = "3.14";
  ezassert( nc.putAtt("v4", "pi", NC_FLOAT, strings, 0, 0) == 0);
  strings[0] = "1.1";
  strings[1] = "2.2";
  strings[2] = "3.3";
  ezassert( nc.putAtt("v4", "double_array", NC_DOUBLE, strings, 0, 2) == 0);
  strings[0] = "-1e30";
  ezassert( nc.putAtt("v5", "_FillValue", NC_FLOAT, strings, 0, 0) == 0);
  strings[0] = "-1e30";
  strings[1] = "1e30";
  ezassert( nc.putAtt("v6", "valid_range", NC_DOUBLE, strings, 0, 1) == 0);
  strings[0] = "quick brown fox jumped";
  ezassert( nc.putAtt("v7", "string", NC_CHAR, strings, 0, 0) == 0);

  ezassert( nc.enddef() == 0 ); 

  ez::ezNcBuffers buf;
  ezassert( nc.reserve(buf) == 0 );
  buf.allocate();

  // Create dummy data.
  ez::seq(buf.uc, buf.nuc);
  ez::seq(buf.c, buf.nc);
  ez::seq(buf.s, buf.ns);
  ez::seq<int>(buf.i, buf.ni, 1);
  ez::seq(buf.f, buf.nf);
  ez::seq(buf.d, buf.nd);
  
  size_t start[4], count[4];
  // Static vars.
  start[0] = start[1] = start[2] = start[3] = 0;
  count[0] = n1; count[1] = n2; count[2] = n3;
  ezassert( nc.write(v1, &buf, start, count) == 0 );
  ezassert( nc.write(v2, &buf, start, count) == 0 );
  ezassert( nc.write(v3, &buf, start, count) == 0 );

  count[0] = 1; count[1] = n1; count[2] = n2; count[3] = n3;
  ezassert( nc.write(v4, &buf, start, count) == 0 );
  ezassert( nc.write(v5, &buf, start, count) == 0 );
  ezassert( nc.write(v6, &buf, start, count) == 0 );
  ezassert( nc.write(v7, &buf, start, count) == 0 );

  start[0] = 1;
  ezassert( nc.write(v4, &buf, start, count) == 0 );
  ezassert( nc.write(v5, &buf, start, count) == 0 );
  ezassert( nc.write(v6, &buf, start, count) == 0 );
  ezassert( nc.write(v7, &buf, start, count) == 0 );
  
  ezassert( nc.close() == 0 );

  return true;
}

bool test_same(ez::ezTestRunner& runner) {
  ezassert( createFile1(runner) );
  std::string fn = filename1(runner);
  ez::ezNcCompare cmp;
  cmp.init(fn, fn);
  
  ezassert( cmp.equals() );
  
  cmp.clear();

  if (runner.clean)
    ezassert( remove(fn.c_str()) == 0 );

  return true;
}

bool test_differ(ez::ezTestRunner& runner) {
  ezassert( createFile1(runner) );
  ezassert( createFile2(runner) );
  std::string fn1 = filename1(runner);
  std::string fn2 = filename2(runner);
  
  ez::ezNcCompare cmp;
  cmp.init(fn1, fn2);
  
  ezassert( cmp.dimEquals() == false );
  ezassert( cmp.equals() == false );
  ezassert( cmp.formatEquals() == false );
  ezassert( cmp.globalAttEquals() == false );
  ezassert( cmp.metaDataEquals() == false );
  ezassert( cmp.varDataEquals() == false );
  
  std::string str1, str2;
  str1 = "d1";
  ezassert( cmp.dimEquals(str1,str1) );
  str1 = "d2";
  ezassert( cmp.dimEquals(str1,str1) );
  str1 = "d2";
  ezassert( cmp.dimEquals(str1,str1) );
  str1 = "dr";
  ezassert( cmp.dimEquals(str1,str1) );
  
  str1 = "Author";
  ezassert( cmp.globalAttEquals(str1,str1) == false );
  str1 = "age";
  ezassert( cmp.globalAttEquals(str1,str1) );
  str1 = "nanosecond";
  ezassert( cmp.globalAttEquals(str1,str1) );
  str1 = "microsecond";
  str2 = "microsecond_64";
  ezassert( cmp.globalAttEquals(str1,str2) );

  ezassert( cmp.varMetaDataEquals() );
  
  str1 = "v1";
  ezassert( cmp.varDataEquals(str1,str1) );
  str1 = "v2";
  ezassert( cmp.varDataEquals(str1,str1) );
  str1 = "v3";
  ezassert( cmp.varDataEquals(str1,str1) );
  str1 = "v4";
  ezassert( cmp.varDataEquals(str1,str1) == false);
  str1 = "v5";
  ezassert( cmp.varDataEquals(str1,str1) );
  str1 = "v6";
  ezassert( cmp.varDataEquals(str1,str1) );
  str1 = "v7";
  ezassert( cmp.varDataEquals(str1,str1) );
  
  cmp.clear();

  if (runner.clean) {
    ezassert( remove(fn1.c_str()) == 0 );
    ezassert( remove(fn2.c_str()) == 0 );
  }

  return true;
}

int main(int argc, const char* argv[]) {
  ez::ezTestRunner runner;
  runner.getopt(argc, argv);
  #define TEST(N) runner.add(&test_##N, "test_"#N);
  TEST(same);
  TEST(differ);
  
  return runner.run();
}