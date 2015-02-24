/*
20111111 rsz Created.
*/

#include "ezTest.hpp"
#include "ezNcSlicer.hpp"

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
  int d1 = nc.createDim("d1", n1);
  ezassert(d1 != -1);
  int d2 = nc.createDim("d2", n2);
  ezassert(d2 != -1);
  int d3 = nc.createDim("d3", n3);
  ezassert(d3 != -1);
  int d4 = nc.createDim("d4", n4);
  ezassert(d4 != -1);
  int dr = nc.createDim("dr", NC_UNLIMITED);
  ezassert(dr != -1);
  nc.nrec = 2;

  int dimids1d[1] = {d1};
  int dimids1d2[1] = {d2};
  int dimids1d3[1] = {d3};
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
  int vd1 = nc.createVar("d1", NC_INT, 1, dimids1d);
  ezassert(vd1 != -1);
  int vd2 = nc.createVar("d2", NC_FLOAT, 1, dimids1d2);
  ezassert(vd2 != -1);
  int vd3 = nc.createVar("d3", NC_DOUBLE, 1, dimids1d3);
  ezassert(vd3 != -1);
  int v4 = nc.createVar("v4", NC_INT, 1, dimids0dt);
  ezassert(v4 != -1);
  int v5 = nc.createVar("v5", NC_FLOAT, 2, dimids1dt);
  ezassert(v5 != -1);
  int v6 = nc.createVar("v6", NC_DOUBLE, 3, dimids2dt);
  ezassert(v6 != -1);
  int v7 = nc.createVar("v7", NC_DOUBLE, 4, dimids3dt);
  ezassert(v7 != -1);

  ezassert( nc.enddef() == 0 ); 

  ez::ezNcBuffers buf;
  size_t start[4], count[4];
  
  // Write 3 coord dim data vars.
  ezassert( nc.init(buf, vd1) == 0 );
  buf.allocate();
  ez::seq<int>(buf.i, 0, 9);
  start[0] = 0;
  count[0] = 10;
  ezassert( nc.write(vd1, &buf, start, count) == 0 );
  // Float sequence.
  ezassert( nc.init(buf, vd2) == 0 );
  buf.allocate();
  ez::seq<float>(buf.f, 0.5, 19.5);
  start[0] = 0;
  count[0] = 20;
  ezassert( nc.write(vd2, &buf, start, count) == 0 );
  // Reversed sequence.
  ezassert( nc.init(buf, vd3) == 0 );
  buf.allocate();
  ez::seq<double>(buf.d, 29.5, 0.5);
  start[0] = 0;
  count[0] = 30;
  ezassert( nc.write(vd3, &buf, start, count) == 0 );
  
  ezassert( nc.init(buf, v7) == 0 );
  buf.allocate();
  buf.zeros();
  // Write last record to force other records to be written to file.
  start[0] = 1;
  start[1] = start[2] = start[3] = 0;
  count[0] = 1;
  count[1] = n1;
  count[2] = n2;
  count[3] = n3;
  ezassert( nc.write(v7, &buf, start, count) == 0 );

  ezassert( nc.close() == 0 );

  return true;
}

bool test_build(ez::ezTestRunner& runner) {
  ez::ezNcSlicer<int> slicer;
  ez::ezNc nc;
  std::vector<int> list;
  int min=0, max=0, stride=0;
  std::string str;
  std::vector<std::string> strings;
  
  ezassert( createFile1(runner) );
  nc.filename = filename1(runner);
  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );

  slicer.build(&nc);
  // --------------------------------------------------------------
  // No slicing tests.
  ezassert( slicer.indexSlices.empty() );
  ezassert( slicer.valueSlices.empty() );
  ezassert( slicer.indexLists.empty() );
  ezassert( slicer.valueLists.empty() );
  ezassert( slicer.slicetypes.size() == 5 );
  ezassert( slicer.slicetypes["d1"] == ez::ezNcSlicer<int>::NOSLICE );
  ezassert( slicer.slicetypes["d2"] == ez::ezNcSlicer<int>::NOSLICE );
  ezassert( slicer.slicetypes["d3"] == ez::ezNcSlicer<int>::NOSLICE );
  ezassert( slicer.slicetypes["d4"] == ez::ezNcSlicer<int>::NOSLICE );
  ezassert( slicer.slicetypes["dr"] == ez::ezNcSlicer<int>::NOSLICE );
  
  ezassert( slicer.getDimSize(0) == 10 );
  ezassert( slicer.getDimSize(1) == 20 );
  ezassert( slicer.getDimSize(2) == 30 );
  ezassert( slicer.getDimSize(3) == 40 );
  ezassert( slicer.getDimSize(4) == 2 );
  
  slicer.getDimIndices(0, list);
  ezassert( list.empty() );
  slicer.getDimIndices(1, list);
  ezassert( list.empty() );
  slicer.getDimIndices(2, list);
  ezassert( list.empty() );
  slicer.getDimIndices(3, list);
  ezassert( list.empty() );
  slicer.getDimIndices(4, list);
  ezassert( list.empty() );

  slicer.getDimMinMaxStride(0, min, max, stride);
  ezassert( min == 0 );
  ezassert( max == 0 );
  ezassert( stride == 1 );
  slicer.getDimMinMaxStride(1, min, max, stride);
  ezassert( min == 0 );
  ezassert( max == 0 );
  ezassert( stride == 1 );
  slicer.getDimMinMaxStride(2, min, max, stride);
  ezassert( min == 0 );
  ezassert( max == 0 );
  ezassert( stride == 1 );
  slicer.getDimMinMaxStride(3, min, max, stride);
  ezassert( min == 0 );
  ezassert( max == 0 );
  ezassert( stride == 1 );
  slicer.getDimMinMaxStride(4, min, max, stride);
  ezassert( min == 0 );
  ezassert( max == 0 );
  ezassert( stride == 1 );
  
  slicer.getVarShape(0, list);
  ezassert( list.size() == 1 );
  ezassert( list[0] == 10 );

  slicer.getVarShape(1, list);
  ezassert( list.size() == 2 );
  ezassert( list[0] == 10 );
  ezassert( list[1] == 20 );

  slicer.getVarShape(2, list);
  ezassert( list.size() == 3 );
  ezassert( list[0] == 10 );
  ezassert( list[1] == 20 );
  ezassert( list[2] == 30 );

  slicer.getVarShape(3, list);
  ezassert( list.size() == 1 );
  ezassert( list[0] == 10 );

  str = "v5";
  slicer.getVarShape(str, list);
  ezassert( list.size() == 2 );
  ezassert( list[0] == 2 );
  ezassert( list[1] == 10 );

  str = "v6";
  slicer.getVarShape(str, list);
  ezassert( list.size() == 3 );
  ezassert( list[0] == 2 );
  ezassert( list[1] == 10 );
  ezassert( list[2] == 20 );

  str = "v7";
  slicer.getVarShape(str, list);
  ezassert( list.size() == 4 );
  ezassert( list[0] == 2 );
  ezassert( list[1] == 10 );
  ezassert( list[2] == 20 );
  ezassert( list[3] == 30 );
  
  ezassert( slicer.isDimSliced(0) == 0 );
  str = "d2";
  ezassert( slicer.isDimSliced(str) == 0 );
  ezassert( slicer.isDimSliced(2) == 0 );
  ezassert( slicer.isDimSliced(3) == 0 );
  ezassert( slicer.isDimSliced(4) == 0 );
  
  ezassert( slicer.isVarSliced(0) == 0 );
  str = "v2";
  ezassert( slicer.isVarSliced(str) == 0 );
  ezassert( slicer.isVarSliced(2) == 0 );
  ezassert( slicer.isVarSliced(3) == 0 );
  ezassert( slicer.isVarSliced(4) == 0 );
  ezassert( slicer.isVarSliced(5) == 0 );
  ezassert( slicer.isVarSliced(6) == 0 );
  
  // --------------------------------------------------------------
  // Do slicing with index lists.
  // Try with bad begin,end list limits.
  str = "d1";
  strings.push_back("0");
  strings.push_back("1");
  strings.push_back("9");
  strings.push_back("5");
  slicer.clear();
  slicer.setDimIndices(str, strings, 0, -1);
  slicer.build(&nc);
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::NOSLICE );
  slicer.clear();
  slicer.setDimIndices(str, strings, -2, 0);
  slicer.build(&nc);
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::NOSLICE );
  // Single index.
  slicer.clear();
  slicer.setDimIndices(str, strings, 0, 0);
  slicer.build(&nc);
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::RANGESLICE );
  ezassert( slicer.isDimSliced(str) );
  slicer.getDimMinMaxStride(str, min, max, stride);
  ezassert( min == 0 );
  ezassert( max == 0 );
  ezassert( stride == 1 );
  
  // Multiple indices unordered.
  slicer.clear();
  slicer.setDimIndices(str, strings, 0, strings.size()-1);
  slicer.build(&nc);
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::INDEXSLICE );
  ezassert( slicer.isDimSliced(str) );
  slicer.getDimIndices(str, list);
  ezassert( list.size() == 4 );
  ezassert( list[0] == 0 );
  ezassert( list[1] == 1 );
  ezassert( list[2] == 9 );
  ezassert( list[3] == 5 );

  // Reverse sequence.
  strings.clear();
  strings.push_back("5");
  strings.push_back("4");
  strings.push_back("3");
  slicer.clear();
  slicer.setDimIndices(str, strings, 0, strings.size()-1);
  slicer.build(&nc);
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::RANGESLICE );
  ezassert( slicer.isDimSliced(str) );
  slicer.getDimMinMaxStride(str, min, max, stride);
  ezassert( min == 5 );
  ezassert( max == 3 );
  ezassert( stride == -1 );
 
  // --------------------------------------------------------------
  // Do slicing with values.
  // Empty list.
  str = "d1";
  strings.clear();
  slicer.clear();
  slicer.setDimValues(str, strings, 0, 0);
  slicer.build(&nc);
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::NOSLICE );
  // Bad list indices.
  strings.push_back("0");
  slicer.clear();
  slicer.setDimValues(str, strings, 0, -1);
  slicer.build(&nc);
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::NOSLICE );
  // Bad list indices.
  slicer.clear();
  slicer.setDimValues(str, strings, -1, 0);
  slicer.build(&nc);
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::NOSLICE );
  // Single value.
  slicer.clear();
  slicer.setDimValues(str, strings, 0, 0);
  slicer.build(&nc);
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::RANGESLICE );
  ezassert( slicer.isDimSliced(str) );
  ezassert( slicer.isVarSliced(str) );
  slicer.getDimMinMaxStride(str, min, max, stride);
  ezassert( min == 0 );
  ezassert( max == 0 );
  ezassert( stride == 1 );
  // Multiple values that are sequential.
  strings.push_back("1");
  strings.push_back("2");
  strings.push_back("3");
  slicer.clear();
  slicer.setDimValues(str, strings, 0, 3);
  slicer.build(&nc);
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::RANGESLICE );
  ezassert( slicer.isDimSliced(str) );
  ezassert( slicer.isVarSliced(str) );
  slicer.getDimMinMaxStride(str, min, max, stride);
  ezassert( min == 0 );
  ezassert( max == 3 );
  ezassert( stride == 1 );
  // Multiple values that are strided.
  strings.clear();
  strings.push_back("1");
  strings.push_back("4");
  strings.push_back("7");
  slicer.clear();
  slicer.setDimValues(str, strings, 0, 2);
  slicer.build(&nc);
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::STRIDESLICE );
  ezassert( slicer.isDimSliced(str) );
  ezassert( slicer.isVarSliced(str) );
  slicer.getDimMinMaxStride(str, min, max, stride);
  ezassert( min == 1 );
  ezassert( max == 7 );
  ezassert( stride == 3 );
  // Multiple values that are scattered indices.
  strings.clear();
  strings.push_back("1");
  strings.push_back("4");
  strings.push_back("9");
  strings.push_back("2");
  slicer.clear();
  slicer.setDimValues(str, strings, 0, 3);
  slicer.build(&nc);
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::INDEXSLICE );
  ezassert( slicer.isDimSliced(str) );
  ezassert( slicer.isVarSliced(str) );
  list.clear();
  slicer.getDimIndices(str, list);
  ezassert( list.size() == 4 );
  ezassert( list[0] == 1 );
  ezassert( list[1] == 4 );
  ezassert( list[2] == 9 );
  ezassert( list[3] == 2 );

  // --------------------------------------------------------------
  // Do slicing with float values.
  str = "d2";
  // Single value.
  strings.clear();
  strings.push_back("0.5");
  slicer.clear();
  slicer.setDimValues(str, strings, 0, 0);
  slicer.build(&nc);
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::RANGESLICE );
  ezassert( slicer.isDimSliced(str) );
  ezassert( slicer.isVarSliced(str) );
  slicer.getDimMinMaxStride(str, min, max, stride);
  ezassert( min == 0 );
  ezassert( max == 0 );
  ezassert( stride == 1 );
  // Multiple values that are sequential.
  strings.push_back("1.5");
  strings.push_back("2.5");
  strings.push_back("3.5");
  slicer.clear();
  slicer.setDimValues(str, strings, 0, 3);
  slicer.build(&nc);
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::RANGESLICE );
  ezassert( slicer.isDimSliced(str) );
  ezassert( slicer.isVarSliced(str) );
  slicer.getDimMinMaxStride(str, min, max, stride);
  ezassert( min == 0 );
  ezassert( max == 3 );
  ezassert( stride == 1 );
  // Multiple values that are strided.
  strings.clear();
  strings.push_back("1.5");
  strings.push_back("4.5");
  strings.push_back("7.5");
  slicer.clear();
  slicer.setDimValues(str, strings, 0, 2);
  slicer.build(&nc);
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::STRIDESLICE );
  ezassert( slicer.isDimSliced(str) );
  ezassert( slicer.isVarSliced(str) );
  slicer.getDimMinMaxStride(str, min, max, stride);
  ezassert( min == 1 );
  ezassert( max == 7 );
  ezassert( stride == 3 );
  // Multiple values that are scattered indices.
  strings.clear();
  strings.push_back("1.5");
  strings.push_back("4.5");
  strings.push_back("9.5");
  strings.push_back("2.5");
  slicer.clear();
  slicer.setDimValues(str, strings, 0, 3);
  slicer.build(&nc);
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::INDEXSLICE );
  ezassert( slicer.isDimSliced(str) );
  ezassert( slicer.isVarSliced(str) );
  list.clear();
  slicer.getDimIndices(str, list);
  ezassert( list.size() == 4 );
  ezassert( list[0] == 1 );
  ezassert( list[1] == 4 );
  ezassert( list[2] == 9 );
  ezassert( list[3] == 2 );
  
  // --------------------------------------------------------------
  // Do slicing with stride indices min/max/stride.
  // Try missing dim name.
  str = "";
  strings.clear();
  strings.push_back("");
  strings.push_back("");
  strings.push_back("");
  slicer.clear();
  slicer.setDimStride(str, strings[0], strings[1], strings[2]);
  slicer.build(&nc);
  str = "d1";
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::NOSLICE );
  // Try missing min,max,stride. Results in complete dim range 0:n-1.
  slicer.clear();
  slicer.setDimStride(str, strings[0], strings[1], strings[2]);
  slicer.build(&nc);
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::RANGESLICE );
  slicer.getDimMinMaxStride(str, min, max, stride);
  ezassert( min == 0 );
  ezassert( max == 9 );
  ezassert( stride == 1 );
  // Just min.
  strings[0] = "5";
  slicer.clear();
  slicer.setDimStride(str, strings[0], strings[1], strings[2]);
  slicer.build(&nc);
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::RANGESLICE );
  slicer.getDimMinMaxStride(str, min, max, stride);
  ezassert( min == 5 );
  ezassert( max == 9 );
  ezassert( stride == 1 );
  // Just max.
  strings[0] = "";
  strings[1] = "5";
  slicer.clear();
  slicer.setDimStride(str, strings[0], strings[1], strings[2]);
  slicer.build(&nc);
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::RANGESLICE );
  slicer.getDimMinMaxStride(str, min, max, stride);
  ezassert( min == 0 );
  ezassert( max == 5 );
  ezassert( stride == 1 );
  // Just min,max.
  strings[0] = "2";
  strings[1] = "5";
  slicer.clear();
  slicer.setDimStride(str, strings[0], strings[1], strings[2]);
  slicer.build(&nc);
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::RANGESLICE );
  slicer.getDimMinMaxStride(str, min, max, stride);
  ezassert( min == 2 );
  ezassert( max == 5 );
  ezassert( stride == 1 );
  // Just stride.
  strings[0] = "";
  strings[1] = "";
  strings[2] = "2";
  slicer.clear();
  slicer.setDimStride(str, strings[0], strings[1], strings[2]);
  slicer.build(&nc);
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::STRIDESLICE );
  slicer.getDimMinMaxStride(str, min, max, stride);
  ezassert( min == 0 );
  ezassert( max == 9 );
  ezassert( stride == 2 );
  // Just min and stride.
  strings[0] = "2";
  strings[1] = "";
  strings[2] = "2";
  slicer.clear();
  slicer.setDimStride(str, strings[0], strings[1], strings[2]);
  slicer.build(&nc);
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::STRIDESLICE );
  slicer.getDimMinMaxStride(str, min, max, stride);
  ezassert( min == 2 );
  ezassert( max == 9 );
  ezassert( stride == 2 );
  // Just max and stride.
  strings[0] = "";
  strings[1] = "5";
  strings[2] = "2";
  slicer.clear();
  slicer.setDimStride(str, strings[0], strings[1], strings[2]);
  slicer.build(&nc);
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::STRIDESLICE );
  slicer.getDimMinMaxStride(str, min, max, stride);
  ezassert( min == 0 );
  ezassert( max == 5 );
  ezassert( stride == 2 );
  // min, max and stride.
  strings[0] = "2";
  strings[1] = "7";
  strings[2] = "2";
  slicer.clear();
  slicer.setDimStride(str, strings[0], strings[1], strings[2]);
  slicer.build(&nc);
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::STRIDESLICE );
  slicer.getDimMinMaxStride(str, min, max, stride);
  ezassert( min == 2 );
  ezassert( max == 7 );
  ezassert( stride == 2 );

  // --------------------------------------------------------------
  // Do slicing with values for min/max and a stride, on increasing float var.
  str = "";
  strings.clear();
  strings.push_back("");
  strings.push_back("");
  strings.push_back("");
  slicer.clear();
  slicer.setDimStrideValues(str, strings[0], strings[1], strings[2]);
  slicer.build(&nc);
  str = "d2";
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::NOSLICE );
  // Try missing min,max,stride. Results in complete dim range 0:n-1.
  slicer.clear();
  slicer.setDimStrideValues(str, strings[0], strings[1], strings[2]);
  slicer.build(&nc);
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::RANGESLICE );
  slicer.getDimMinMaxStride(str, min, max, stride);
  ezassert( min == 0 );
  ezassert( max == 19 );
  ezassert( stride == 1 );
  // Just min. 
  strings[0] = "5.1";
  strings[1] = "";
  strings[2] = "";
  slicer.clear();
  slicer.setDimStrideValues(str, strings[0], strings[1], strings[2]);
  slicer.build(&nc);
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::RANGESLICE );
  slicer.getDimMinMaxStride(str, min, max, stride);
  ezassert( min == 5 );
  ezassert( max == 19 );
  ezassert( stride == 1 );
  // Just max. 
  strings[0] = "";
  strings[1] = "5";
  slicer.clear();
  slicer.setDimStrideValues(str, strings[0], strings[1], strings[2]);
  slicer.build(&nc);
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::RANGESLICE );
  slicer.getDimMinMaxStride(str, min, max, stride);
  ezassert( min == 0 );
  ezassert( max == 4 );
  ezassert( stride == 1 );
  // Just min,max.
  strings[0] = "5";
  strings[1] = "10";
  slicer.clear();
  slicer.setDimStrideValues(str, strings[0], strings[1], strings[2]);
  slicer.build(&nc);
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::RANGESLICE );
  slicer.getDimMinMaxStride(str, min, max, stride);
  ezassert( min == 5 );
  ezassert( max == 9 );
  ezassert( stride == 1 );
  // Just stride.
  strings[0] = "";
  strings[1] = "";
  strings[2] = "2";
  slicer.clear();
  slicer.setDimStrideValues(str, strings[0], strings[1], strings[2]);
  slicer.build(&nc);
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::STRIDESLICE );
  slicer.getDimMinMaxStride(str, min, max, stride);
  ezassert( min == 0 );
  ezassert( max == 19 );
  ezassert( stride == 2 );
  // Just min and stride.
  strings[0] = "5";
  strings[1] = "";
  strings[2] = "2";
  slicer.clear();
  slicer.setDimStrideValues(str, strings[0], strings[1], strings[2]);
  slicer.build(&nc);
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::STRIDESLICE );
  slicer.getDimMinMaxStride(str, min, max, stride);
  ezassert( min == 5 );
  ezassert( max == 19 );
  ezassert( stride == 2 );
  // Just max and stride.
  strings[0] = "";
  strings[1] = "5";
  strings[2] = "2";
  slicer.clear();
  slicer.setDimStrideValues(str, strings[0], strings[1], strings[2]);
  slicer.build(&nc);
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::STRIDESLICE );
  slicer.getDimMinMaxStride(str, min, max, stride);
  ezassert( min == 0 );
  ezassert( max == 4 );
  ezassert( stride == 2 );
  // min, max and stride.
  strings[0] = "5";
  strings[1] = "10";
  strings[2] = "2";
  slicer.clear();
  slicer.setDimStrideValues(str, strings[0], strings[1], strings[2]);
  slicer.build(&nc);
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::STRIDESLICE );
  slicer.getDimMinMaxStride(str, min, max, stride);
  ezassert( min == 5 );
  ezassert( max == 9 );
  ezassert( stride == 2 );

  // --------------------------------------------------------------
  // Do slicing with values for min/max and a stride, on decreasing double var.
  str = "";
  strings.clear();
  strings.push_back("");
  strings.push_back("");
  strings.push_back("");
  slicer.clear();
  slicer.setDimStrideValues(str, strings[0], strings[1], strings[2]);
  slicer.build(&nc);
  str = "d3";
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::NOSLICE );
  // Try missing min,max,stride. Results in complete dim range 0:n-1.
  slicer.clear();
  slicer.setDimStrideValues(str, strings[0], strings[1], strings[2]);
  slicer.build(&nc);
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::RANGESLICE );
  slicer.getDimMinMaxStride(str, min, max, stride);
  ezassert( min == 0 );
  ezassert( max == 29 );
  ezassert( stride == 1 );
  // Just min. But data is decreasing, so expect index range for data +inf>x>=5.1 .
  strings[0] = "5.1";
  strings[1] = "";
  strings[2] = "";
  slicer.clear();
  slicer.setDimStrideValues(str, strings[0], strings[1], strings[2]);
  slicer.build(&nc);
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::RANGESLICE );
  slicer.getDimMinMaxStride(str, min, max, stride);
  ezassert( min == 0 );
  ezassert( max == 24 );
  ezassert( stride == 1 );
  // Just max. In decreasing data, this means min from beginning of array.
  strings[0] = "";
  strings[1] = "5";
  slicer.clear();
  slicer.setDimStrideValues(str, strings[0], strings[1], strings[2]);
  slicer.build(&nc);
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::RANGESLICE );
  slicer.getDimMinMaxStride(str, min, max, stride);
  ezassert( min == 25 );
  ezassert( max == 29 );
  ezassert( stride == 1 );
  // Just min,max.
  strings[0] = "5";
  strings[1] = "20";
  slicer.clear();
  slicer.setDimStrideValues(str, strings[0], strings[1], strings[2]);
  slicer.build(&nc);
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::RANGESLICE );
  slicer.getDimMinMaxStride(str, min, max, stride);
  ezassert( min == 10 );
  ezassert( max == 24 );
  ezassert( stride == 1 );
  // Just stride.
  strings[0] = "";
  strings[1] = "";
  strings[2] = "2";
  slicer.clear();
  slicer.setDimStrideValues(str, strings[0], strings[1], strings[2]);
  slicer.build(&nc);
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::STRIDESLICE );
  slicer.getDimMinMaxStride(str, min, max, stride);
  ezassert( min == 0 );
  ezassert( max == 29 );
  ezassert( stride == 2 );
  // Just min and stride.
  strings[0] = "5";
  strings[1] = "";
  strings[2] = "2";
  slicer.clear();
  slicer.setDimStrideValues(str, strings[0], strings[1], strings[2]);
  slicer.build(&nc);
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::STRIDESLICE );
  slicer.getDimMinMaxStride(str, min, max, stride);
  ezassert( min == 0 );
  ezassert( max == 24 );
  ezassert( stride == 2 );
  // Just max and stride.
  strings[0] = "";
  strings[1] = "5";
  strings[2] = "2";
  slicer.clear();
  slicer.setDimStrideValues(str, strings[0], strings[1], strings[2]);
  slicer.build(&nc);
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::STRIDESLICE );
  slicer.getDimMinMaxStride(str, min, max, stride);
  ezassert( min == 25 );
  ezassert( max == 29 );
  ezassert( stride == 2 );
  // min, max and stride.
  strings[0] = "5";
  strings[1] = "20";
  strings[2] = "2";
  slicer.clear();
  slicer.setDimStrideValues(str, strings[0], strings[1], strings[2]);
  slicer.build(&nc);
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::STRIDESLICE );
  slicer.getDimMinMaxStride(str, min, max, stride);
  ezassert( min == 10 );
  ezassert( max == 24 );
  ezassert( stride == 2 );

  // --------------------------------------------------------------
  // Do slicing with fortran 1-based index lists.
  // Try with bad begin,end list limits.
  str = "d1";
  slicer.offset = -1;
  strings.clear();
  strings.push_back("1");
  strings.push_back("2");
  strings.push_back("10");
  strings.push_back("6");
  slicer.clear();
  slicer.setDimIndices(str, strings, 0, -1);
  slicer.build(&nc);
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::NOSLICE );
  slicer.clear();
  slicer.setDimIndices(str, strings, -2, 0);
  slicer.build(&nc);
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::NOSLICE );
  // Single index.
  slicer.clear();
  slicer.setDimIndices(str, strings, 0, 0);
  slicer.build(&nc);
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::RANGESLICE );
  ezassert( slicer.isDimSliced(str) );
  slicer.getDimMinMaxStride(str, min, max, stride);
  ezassert( min == 0 );
  ezassert( max == 0 );
  ezassert( stride == 1 );
  // Multiple indices.
  slicer.clear();
  slicer.setDimIndices(str, strings, 0, strings.size()-1);
  slicer.build(&nc);
  ezassert( slicer.slicetypes[str] == ez::ezNcSlicer<int>::INDEXSLICE );
  ezassert( slicer.isDimSliced(str) );
  slicer.getDimIndices(str, list);
  ezassert( list.size() == 4 );
  // Should be same order as input indices.
  ezassert( list[0] == 0 );
  ezassert( list[1] == 1 );
  ezassert( list[2] == 9 );
  ezassert( list[3] == 5 );

  ezassert( nc.close() == 0 );

  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );
  
  return true;
}

int main(int argc, const char* argv[]) {
  ez::ezTestRunner runner;
  runner.getopt(argc, argv);
  #define TEST(N) runner.add(&test_##N, "test_"#N);
  TEST(build);

  return runner.run();
}