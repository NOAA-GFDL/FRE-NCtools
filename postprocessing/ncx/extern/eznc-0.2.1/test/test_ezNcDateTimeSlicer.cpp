/*
20120103 rsz Created.
*/

#include "ezTest.hpp"
#include "ezNcDateTimeSlicer.hpp"

std::string createFilename(ez::ezTestRunner& runner, const char* id) {
  std::string fname = runner.tmpdir;
  fname.append(id);
  fname.append(".nc");
  return fname;
}

bool createFile(ez::ezTestRunner& runner, std::string& fname, std::string& units, int toffset) {
  ez::ezNc nc;
  nc.filename = fname;  
  ezassert( nc.createClobber() == 0 );

  int n = 40;
  int dr = nc.createDim("time", n, true);
  ezassert(dr != -1);

  int dimids[1] = {dr};
  int v = nc.createVar("time", NC_DOUBLE, 1, dimids);
  ezassert(v != -1);
  
  std::vector<std::string> strings;
  strings.resize(1);
  strings[0] = "time";
  ezassert( nc.putAtt("time", "standard_name", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = units;
  ezassert( nc.putAtt("time", "units", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "noleap";
  ezassert( nc.putAtt("time", "calendar", NC_CHAR, strings, 0, 0) == 0);

  ezassert( nc.enddef() == 0 ); 

  ez::ezNcBuffers buf;
  ezassert( nc.reserve(buf, v) == 0 );
  buf.allocate();

  size_t start, count=1;
  for(int i=0; i < n; ++i) {
    start = i;
    buf.d[0] = i + toffset;
    ezassert( nc.write(v, &buf, &start, &count) == 0 );
  }
  
  ezassert( nc.close() == 0 );
}

// Prints all time axis values as strings.
void print_datetimes(ez::ezNcDateTime &dt) {
  int n = dt.values.size();
  int i, year, month, day, hr, min;
  float sec;
  std::vector<double> rawValues;
  
  dt.rawValues.getValues<double>(dt.var->type, rawValues);
  
  for(i=0; i < n; ++i) {
    ez::ezDateTime::fromDouble(dt.values[i], year, month, day, hr, min, sec, dt.calendar, dt.tzHours, dt.tzMinutes);
    printf("[%d] %g = %04d-%02d-%02d %02d:%02d:%02d UTC\n", i, rawValues[i], year, month, day, hr, min, (int)sec);
  }
}

bool test_build(ez::ezTestRunner& runner) {
  ez::ezNcDateTime dt;
  ez::ezNc nc;
  ez::ezNcDateTimeSlicer<int> slicer;
  std::string fname, units, str;
  int iBegin, iEnd, iStride, n;
  
  fname = createFilename(runner, "tmp1");
  units = "years since 1980-01-01";
  createFile(runner, fname, units, 0);

  nc.filename = fname;
  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );
  ezassert( dt.build(&nc) == 0 );

  //--------------------------------------
  // Slice before valid range.
  slicer.begin = "1800-01-01";
  slicer.end = "1900-12-31";
  slicer.build(&dt);
  
  str = "time";
  n = slicer.getDimSize(str);

  ezassert( slicer.isDimSliced(str) == 0 );

  //--------------------------------------
  // Begin before first value.
  slicer.begin = "1900-01-01";
  slicer.end = "1999-12-31";
  slicer.build(&dt);
  //print_datetimes(dt);
  
  n = slicer.getDimSize(str);
  //printf("%s:%d n = %d\n", __FILE__, __LINE__, n);
  slicer.getDimMinMaxStride(str, iBegin, iEnd, iStride);
  //printf("%s:%d iBegin = %d, iEnd = %d, iStride = %d\n", __FILE__, __LINE__, iBegin, iEnd, iStride);

  ezassert( slicer.isDimSliced(str) );
  ezassert( n == 20 );
  ezassert( iBegin == 0 );
  ezassert( iEnd == 19 );

  //--------------------------------------
  // Slice on begin edge
  slicer.begin = "1980-01-01";
  slicer.end = "1980-12-31";
  slicer.build(&dt);
  
  n = slicer.getDimSize(str);
  slicer.getDimMinMaxStride(str, iBegin, iEnd, iStride);

  ezassert( slicer.isDimSliced(str) );
  ezassert( n == 1 );
  ezassert( iBegin == 0 );
  ezassert( iEnd == 0 );

  //--------------------------------------
  // Slice in the middle.
  slicer.clear();
  slicer.begin = "1999-01-01";
  slicer.end = "2010-1-1";
  slicer.build(&dt);
  
  n = slicer.getDimSize(str);
  slicer.getDimMinMaxStride(str, iBegin, iEnd, iStride);

  ezassert( slicer.isDimSliced(str) );
  ezassert( n == 12 );
  ezassert( iBegin == 19 );
  ezassert( iEnd == 30 );

  //--------------------------------------
  // Slice on end edge.
  slicer.clear();
  slicer.begin = "2019-01-01";
  slicer.end = "2019-1-1";
  slicer.build(&dt);
  
  n = slicer.getDimSize(str);
  slicer.getDimMinMaxStride(str, iBegin, iEnd, iStride);

  ezassert( slicer.isDimSliced(str) );
  ezassert( n == 1 );
  ezassert( iBegin == 39 );
  ezassert( iEnd == 39 );
  
  //--------------------------------------
  // Slice past end edge.
  slicer.clear();
  slicer.begin = "2010-01-01";
  slicer.end = "2100-1-1";
  slicer.build(&dt);
  
  n = slicer.getDimSize(str);
  slicer.getDimMinMaxStride(str, iBegin, iEnd, iStride);

  ezassert( slicer.isDimSliced(str) );
  ezassert( n == 10 );
  ezassert( iBegin == 30 );
  ezassert( iEnd == 39 );
  
  ezassert( nc.close() == 0 );

  if (runner.clean) {
    ezassert( remove(fname.c_str()) == 0 );
  }
  
  return true;
}

int main(int argc, const char* argv[]) {
  ez::ezTestRunner runner;
  runner.getopt(argc, argv);
  #define TEST(N) runner.add(&test_##N, "test_"#N);
  TEST(build);

  return runner.run();
}