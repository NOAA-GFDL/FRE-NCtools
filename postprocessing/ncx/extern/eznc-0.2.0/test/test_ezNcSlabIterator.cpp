/*
20111115 rsz Created.
*/

#include "ezTest.hpp"
#include "ezNc.hpp"
#include "ezNcSlicer.hpp"
#include "ezNcSlabIterator.hpp"

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

bool test_noslice(ez::ezTestRunner& runner) {
  ez::ezNcSlicer<int> slicer;
  ez::ezNc nc;
  ez::ezNcSlabIterator<int> it;
  ez::ezNcSlab<size_t> slab;
  ez::ezNcVar *var;
  std::vector<std::string> strings;
  std::string str;

  ezassert( createFile1(runner) );
  nc.filename = filename1(runner);
  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );
  slicer.build(&nc);
  
  //-----------------------------------------
  // Create with empty slicer, so slice is for entire var.
  
  var = nc.getVar("v1");
  ezassert( var );
  it.setVar(var, &slicer);
  it.getSlab(&slab); 
  ezassert( slab.start.size() == 1 );
  ezassert( slab.start[0] == 0 );
  ezassert( slab.count.size() == 1 );
  ezassert( slab.count[0] == 10 );
  
  var = nc.getVar("v2");
  ezassert( var );
  it.setVar(var, &slicer);
  it.getSlab(&slab); 
  ezassert( slab.start.size() == 2 );
  ezassert( slab.start[0] == 0 );
  ezassert( slab.start[1] == 0 );
  ezassert( slab.count.size() == 2 );
  ezassert( slab.count[0] == 10 );
  ezassert( slab.count[1] == 20 );

  var = nc.getVar("v3");
  ezassert( var );
  it.setVar(var, &slicer);
  it.getSlab(&slab); 
  ezassert( slab.start.size() == 3 );
  ezassert( slab.start[0] == 0 );
  ezassert( slab.start[1] == 0 );
  ezassert( slab.start[2] == 0 );
  ezassert( slab.count.size() == 3 );
  ezassert( slab.count[0] == 10 );
  ezassert( slab.count[1] == 20 );
  ezassert( slab.count[2] == 30 );
  
  var = nc.getVar("v4");
  ezassert( var );
  it.setVar(var, &slicer);
  it.getSlab(&slab); 
  ezassert( slab.start.size() == 1 );
  ezassert( slab.start[0] == 0 );
  ezassert( slab.count.size() == 1 );
  ezassert( slab.count[0] == 2 );

  var = nc.getVar("v5");
  ezassert( var );
  it.setVar(var, &slicer);
  it.getSlab(&slab); 
  ezassert( slab.start.size() == 2 );
  ezassert( slab.start[0] == 0 );
  ezassert( slab.start[1] == 0 );
  ezassert( slab.count.size() == 2 );
  ezassert( slab.count[0] == 2 );
  ezassert( slab.count[1] == 10 );

  var = nc.getVar("v6");
  ezassert( var );
  it.setVar(var, &slicer);
  it.getSlab(&slab); 
  ezassert( slab.start.size() == 3 );
  ezassert( slab.start[0] == 0 );
  ezassert( slab.start[1] == 0 );
  ezassert( slab.start[2] == 0 );
  ezassert( slab.count.size() == 3 );
  ezassert( slab.count[0] == 2 );
  ezassert( slab.count[1] == 10 );
  ezassert( slab.count[2] == 20 );

  var = nc.getVar("v7");
  ezassert( var );
  it.setVar(var, &slicer);
  it.getSlab(&slab); 
  ezassert( slab.start.size() == 4 );
  ezassert( slab.start[0] == 0 );
  ezassert( slab.start[1] == 0 );
  ezassert( slab.start[2] == 0 );
  ezassert( slab.start[3] == 0 );
  ezassert( slab.count.size() == 4 );
  ezassert( slab.count[0] == 2 );
  ezassert( slab.count[1] == 10 );
  ezassert( slab.count[2] == 20 );
  ezassert( slab.count[3] == 30 );
  
  nc.close();
  
  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );
  
  return true;
}

bool test_incIndexLists(ez::ezTestRunner& runner) {
  ez::ezNcSlicer<int> slicer;
  ez::ezNc nc;
  ez::ezNcSlabIterator<int> it;
  ez::ezNcSlab<size_t> slab;
  ez::ezNcVar *var;
  std::vector<std::string> strings;
  std::string str;

  ezassert( createFile1(runner) );
  nc.filename = filename1(runner);
  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );

  // ---------------------------------------------
  // Test incrementing all dims.
  // Do slicing with index lists 1d.
  str = "d1";
  // Single index.
  strings.clear();
  strings.push_back("0");
  slicer.clear();
  slicer.setDimIndices(str, strings, 0, 0);
  slicer.build(&nc);
  var = nc.getVar("v1");
  ezassert( var );
  it.setVar(var, &slicer);
  it.getSlab(&slab); 
  ezassert( slab.start.size() == 1 );
  ezassert( slab.start[0] == 0 );
  ezassert( slab.count.size() == 1 );
  ezassert( slab.count[0] == 1 );
  ezassert( it.isDone() == 0 );
  it.inc();
  ezassert( it.isDone() );
  // Multiple sequential indices.
  strings.clear();
  strings.push_back("3");
  strings.push_back("2");
  strings.push_back("4");
  slicer.clear();
  slicer.setDimIndices(str, strings, 0, 2);
  slicer.build(&nc);
  var = nc.getVar("v1");
  ezassert( var );
  it.setVar(var, &slicer);
  it.getSlab(&slab); 
  ezassert( slab.start.size() == 1 );
  ezassert( slab.start[0] == 3 );
  ezassert( slab.count.size() == 1 );
  ezassert( slab.count[0] == 1 );
  ezassert( it.isDone() == 0);
  it.inc();
  it.inc();
  it.inc();
  ezassert( it.isDone() );
  // Multiple strided indices.
  strings.clear();
  strings.push_back("7");
  strings.push_back("5");
  strings.push_back("3");
  slicer.clear();
  slicer.setDimIndices(str, strings, 0, 2);
  slicer.build(&nc);
  var = nc.getVar("v1");
  ezassert( var );
  it.setVar(var, &slicer);
  it.getSlab(&slab); 
  ezassert( slab.start.size() == 1 );
  ezassert( slab.start[0] == 7 );
  ezassert( slab.count.size() == 1 );
  ezassert( slab.count[0] == 1 );
  ezassert( it.isDone() == false );
  it.inc();
  it.getSlab(&slab); 
  ezassert( it.isDone() == false );
  ezassert( slab.start[0] == 5 );
  ezassert( slab.count[0] == 1 );
  it.inc();
  ezassert( it.isDone() == 0 );
  it.getSlab(&slab); 
  ezassert( slab.start[0] == 3 );
  ezassert( slab.count[0] == 1 );
  it.inc();
  ezassert( it.isDone() );
  // Multiple scattered indices.
  strings.clear();
  strings.push_back("8");
  strings.push_back("6");
  strings.push_back("2");
  slicer.clear();
  slicer.setDimIndices(str, strings, 0, 2);
  slicer.build(&nc);
  var = nc.getVar("v1");
  ezassert( var );
  it.setVar(var, &slicer);
  it.getSlab(&slab); 
  ezassert( slab.start.size() == 1 );
  ezassert( slab.start[0] == 8 );
  ezassert( slab.count.size() == 1 );
  ezassert( slab.count[0] == 1 );
  ezassert( it.isDone() == false );
  it.inc();
  it.getSlab(&slab); 
  ezassert( it.isDone() == false );
  ezassert( slab.start[0] == 6 );
  ezassert( slab.count[0] == 1 );
  it.inc();
  ezassert( it.isDone() == 0 );
  it.getSlab(&slab); 
  ezassert( slab.start[0] == 2 );
  ezassert( slab.count[0] == 1 );
  it.inc();
  ezassert( it.isDone() );

  // Do slicing with index lists 2d.
  str = "d1";
  // Single index.
  strings.clear();
  strings.push_back("0");
  slicer.clear();
  slicer.setDimIndices(str, strings, 0, 0);
  str = "d2";
  slicer.setDimIndices(str, strings, 0, 0);
  slicer.build(&nc);
  var = nc.getVar("v2");
  ezassert( var );
  it.setVar(var, &slicer);
  it.getSlab(&slab); 
  ezassert( slab.start.size() == 2 );
  ezassert( slab.start[0] == 0 );
  ezassert( slab.start[1] == 0 );
  ezassert( slab.count.size() == 2 );
  ezassert( slab.count[0] == 1 );
  ezassert( slab.count[1] == 1 );
  ezassert( it.isDone() == 0 );
  it.inc();
  ezassert( it.isDone() );
  // Multiple indices.
  strings.clear();
  strings.push_back("3");
  strings.push_back("2");
  strings.push_back("4");
  slicer.clear();
  str = "d1";
  slicer.setDimIndices(str, strings, 0, 2);
  str = "d2";
  slicer.setDimIndices(str, strings, 0, 2);
  slicer.build(&nc);
  var = nc.getVar("v2");
  ezassert( var );
  it.setVar(var, &slicer);
  it.getSlab(&slab); 
  ezassert( slab.start.size() == 2 );
  ezassert( slab.start[0] == 3 );
  ezassert( slab.start[1] == 3 );
  ezassert( slab.count.size() == 2 );
  ezassert( slab.count[0] == 1 );
  ezassert( slab.count[1] == 1 );
  ezassert( it.isDone() == 0 );
  // Multiple strided indices.
  strings.clear();
  strings.push_back("3");
  strings.push_back("5");
  strings.push_back("7");
  slicer.clear();
  str = "d1";
  slicer.setDimIndices(str, strings, 0, 2);
  str = "d2";
  slicer.setDimIndices(str, strings, 0, 2);
  slicer.build(&nc);
  var = nc.getVar("v2");
  ezassert( var );
  it.setVar(var, &slicer);
  it.getSlab(&slab); 
  ezassert( slab.start.size() == 2 );
  ezassert( slab.start[0] == 3 );
  ezassert( slab.start[1] == 3 );
  ezassert( slab.count.size() == 2 );
  ezassert( slab.count[0] == 1 );
  ezassert( slab.count[1] == 1 );
  ezassert( it.isDone() == false );
  it.inc();
  it.getSlab(&slab); 
  ezassert( it.isDone() == false );
  ezassert( slab.start[0] == 3 );
  ezassert( slab.start[1] == 5 );
  it.inc();
  it.getSlab(&slab); 
  ezassert( it.isDone() == false );
  ezassert( slab.start[0] == 3 );
  ezassert( slab.start[1] == 7 );
  it.inc();
  it.getSlab(&slab); 
  ezassert( it.isDone() == false );
  ezassert( slab.start[0] == 5 );
  ezassert( slab.start[1] == 3 );
  it.inc();
  it.getSlab(&slab); 
  ezassert( it.isDone() == false );
  ezassert( slab.start[0] == 5 );
  ezassert( slab.start[1] == 5 );
  while(! it.isDone() ) it.inc();
  it.getSlab(&slab); 
  ezassert( slab.start[0] == 7 );
  ezassert( slab.start[1] == 7 );
  // Multiple scattered indices.
  strings.clear();
  strings.push_back("8");
  strings.push_back("6");
  strings.push_back("2");
  slicer.clear();
  str = "d1";
  slicer.setDimIndices(str, strings, 0, 2);
  str = "d2";
  slicer.setDimIndices(str, strings, 0, 2);
  slicer.build(&nc);
  var = nc.getVar("v2");
  ezassert( var );
  it.setVar(var, &slicer);
  it.getSlab(&slab); 
  ezassert( slab.start.size() == 2 );
  ezassert( slab.start[0] == 8 );
  ezassert( slab.start[1] == 8 );
  ezassert( slab.count.size() == 2 );
  ezassert( slab.count[0] == 1 );
  ezassert( slab.count[1] == 1 );
  ezassert( it.isDone() == false );
  it.inc();
  it.getSlab(&slab); 
  ezassert( it.isDone() == false );
  ezassert( slab.start[0] == 8 );
  ezassert( slab.start[1] == 6 );
  it.inc();
  ezassert( it.isDone() == false );
  it.getSlab(&slab); 
  ezassert( slab.start[0] == 8 );
  ezassert( slab.start[1] == 2 );
  while(! it.isDone() ) it.inc();
  it.getSlab(&slab); 
  ezassert( slab.start[0] == 2 );
  ezassert( slab.start[1] == 2 );

  // Do slicing with index lists 3d.
  // Single index.
  strings.clear();
  strings.push_back("0");
  slicer.clear();
  str = "d1";
  slicer.setDimIndices(str, strings, 0, 0);
  str = "d2";
  slicer.setDimIndices(str, strings, 0, 0);
  str = "d3";
  slicer.setDimIndices(str, strings, 0, 0);
  slicer.build(&nc);
  var = nc.getVar("v3");
  ezassert( var );
  it.setVar(var, &slicer);
  it.getSlab(&slab); 
  ezassert( slab.start.size() == 3 );
  ezassert( slab.start[0] == 0 );
  ezassert( slab.start[1] == 0 );
  ezassert( slab.start[2] == 0 );
  ezassert( slab.count.size() == 3 );
  ezassert( slab.count[0] == 1 );
  ezassert( slab.count[1] == 1 );
  ezassert( slab.count[2] == 1 );
  ezassert( it.isDone() == 0 );
  it.inc();
  ezassert( it.isDone() );
  
  // Multiple sequential indices.
  strings.clear();
  strings.push_back("2");
  strings.push_back("3");
  strings.push_back("4");
  slicer.clear();
  str = "d1";
  slicer.setDimIndices(str, strings, 0, 2);
  str = "d2";
  slicer.setDimIndices(str, strings, 0, 2);
  str = "d3";
  slicer.setDimIndices(str, strings, 0, 2);
  slicer.build(&nc);
  var = nc.getVar("v3");
  ezassert( var );
  it.setVar(var, &slicer);
  it.getSlab(&slab); 
  ezassert( slab.start.size() == 3 );
  ezassert( slab.start[0] == 2 );
  ezassert( slab.start[1] == 2 );
  ezassert( slab.start[2] == 2 );
  ezassert( slab.count.size() == 3 );
  ezassert( slab.count[0] == 3 );
  ezassert( slab.count[1] == 3 );
  ezassert( slab.count[2] == 3 );
  ezassert( it.isDone() == 0 );
  it.inc();
  ezassert( it.isDone() );
  // Multiple scattered indices.
  strings.clear();
  strings.push_back("8");
  strings.push_back("6");
  strings.push_back("2");
  slicer.clear();
  str = "d1";
  slicer.setDimIndices(str, strings, 0, 2);
  str = "d2";
  slicer.setDimIndices(str, strings, 0, 2);
  str = "d3";
  slicer.setDimIndices(str, strings, 0, 2);
  slicer.build(&nc);
  var = nc.getVar("v3");
  ezassert( var );
  it.setVar(var, &slicer);
  it.getSlab(&slab); 
  ezassert( slab.start.size() == 3 );
  ezassert( slab.start[0] == 8 );
  ezassert( slab.start[1] == 8 );
  ezassert( slab.start[2] == 8 );
  ezassert( slab.count.size() == 3 );
  ezassert( slab.count[0] == 1 );
  ezassert( slab.count[1] == 1 );
  ezassert( slab.count[2] == 1 );
  ezassert( it.isDone() == false );
  it.inc();
  it.getSlab(&slab); 
  ezassert( it.isDone() == false );
  ezassert( slab.start[0] == 8 );
  ezassert( slab.start[1] == 8 );
  ezassert( slab.start[2] == 6 );
  it.inc();
  ezassert( it.isDone() == false );
  it.getSlab(&slab); 
  ezassert( slab.start[0] == 8 );
  ezassert( slab.start[1] == 8 );
  ezassert( slab.start[2] == 2 );
  while(! it.isDone() ) it.inc();
  it.getSlab(&slab); 
  ezassert( slab.start[0] == 2 );
  ezassert( slab.start[1] == 2 );
  ezassert( slab.start[2] == 2 );

  // Incremenent all dims. Record vars only.
  // Do slicing with index lists 1d.
  // Single index.
  strings.clear();
  strings.push_back("0");
  slicer.clear();
  str = "d4";
  slicer.setDimIndices(str, strings, 0, 0);
  slicer.build(&nc);
  var = nc.getVar("v4");
  ezassert( var );
  it.setVar(var, &slicer);
  it.getSlab(&slab); 
  ezassert( slab.start.size() == 1 );
  ezassert( slab.start[0] == 0 );
  ezassert( slab.count.size() == 1 );
  ezassert( slab.count[0] == 2 );
  ezassert( it.isDone() == 0 );
  it.inc();
  ezassert( it.isDone() );
  
  // Multiple sequential indices.
  strings.clear();
  strings.push_back("0");
  strings.push_back("1");
  slicer.clear();
  slicer.setDimIndices(str, strings, 0, 2);
  slicer.build(&nc);
  var = nc.getVar("v4");
  ezassert( var );
  it.setVar(var, &slicer);
  it.getSlab(&slab); 
  ezassert( slab.start.size() == 1 );
  ezassert( slab.start[0] == 0 );
  ezassert( slab.count.size() == 1 );
  ezassert( slab.count[0] == 2 );
  ezassert( it.isDone() == 0 );
  it.inc();
  ezassert( it.isDone() );

  // 4d record var with mixed slicing per dim:
  // none,index,range,stride
  strings.clear();
  strings.push_back("3");
  strings.push_back("4");
  strings.push_back("9");
  slicer.clear();
  str = "d1";
  slicer.setDimIndices(str, strings, 0, 2);
  strings.clear();
  strings.push_back("9");
  strings.push_back("10");
  strings.push_back("11");
  str = "d2";
  slicer.setDimIndices(str, strings, 0, 2);
  strings.clear();
  strings.push_back("9");
  strings.push_back("19");
  strings.push_back("29");
  str = "d3";
  slicer.setDimIndices(str, strings, 0, 2);
  slicer.build(&nc);
  var = nc.getVar("v7");
  ezassert( var );
  it.setVar(var, &slicer);
  it.getSlab(&slab); 
  ezassert( slab.start.size() == 4 );
  ezassert( slab.start[0] == 0 );
  ezassert( slab.start[1] == 3 );
  ezassert( slab.start[2] == 9 );
  ezassert( slab.start[3] == 9 );
  ezassert( slab.count.size() == 4 );
  ezassert( slab.count[0] == 2 );
  ezassert( slab.count[1] == 1 );
  ezassert( slab.count[2] == 3 );
  ezassert( slab.count[3] == 1 );
  ezassert( it.isDone() == false );
  while(! it.isDone() ) it.inc();
  it.getSlab(&slab); 
  ezassert( slab.start[0] == 0 );
  ezassert( slab.start[1] == 9 );
  ezassert( slab.start[2] == 9 );
  ezassert( slab.start[3] == 29 );

  nc.close();
  
  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );
  
  return true;
}

bool test_incValueLists(ez::ezTestRunner& runner) {
  ez::ezNcSlicer<int> slicer;
  ez::ezNc nc;
  ez::ezNcSlabIterator<int> it;
  ez::ezNcSlab<size_t> slab;
  ez::ezNcVar *var;
  std::vector<std::string> strings;
  std::string str;

  ezassert( createFile1(runner) );
  nc.filename = filename1(runner);
  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );

  // ---------------------------------------------
  // 4d mixed dim value slicing
  // none,index,range,stride
  strings.clear();
  strings.push_back("3");
  strings.push_back("4");
  strings.push_back("9");
  slicer.clear();
  str = "d1";
  slicer.setDimValues(str, strings, 0, 2);
  strings.clear();
  strings.push_back("9.5");
  strings.push_back("10.5");
  strings.push_back("11.5");
  str = "d2";
  slicer.setDimValues(str, strings, 0, 2);
  strings.clear();
  strings.push_back("19.5");
  strings.push_back("14.5");
  strings.push_back("9.5");
  str = "d3"; // Remember that d3 values are decreasing sequence.
  slicer.setDimValues(str, strings, 0, 2);
  slicer.build(&nc);
  var = nc.getVar("v7");
  ezassert( var );
  it.setVar(var, &slicer);
  it.getSlab(&slab); 
  ezassert( slab.start.size() == 4 );
  ezassert( slab.start[0] == 0 );
  ezassert( slab.start[1] == 3 );
  ezassert( slab.start[2] == 9 );
  ezassert( slab.start[3] == 10 );
  ezassert( slab.count.size() == 4 );
  ezassert( slab.count[0] == 2 );
  ezassert( slab.count[1] == 1 );
  ezassert( slab.count[2] == 3 );
  ezassert( slab.count[3] == 1 );
  ezassert( it.isDone() == false );
  while(! it.isDone() ) it.inc();
  it.getSlab(&slab); 
  ezassert( slab.start[0] == 0 );
  ezassert( slab.start[1] == 9 );
  ezassert( slab.start[2] == 9 );
  ezassert( slab.start[3] == 20 );

  nc.close();
  
  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );
  
  return true;
}
  
bool test_incSeparateRecStatic(ez::ezTestRunner& runner) {
  ez::ezNcSlicer<int> slicer;
  ez::ezNc nc;
  ez::ezNcSlabIterator<int> it;
  ez::ezNcSlab<size_t> slab;
  ez::ezNcVar *var;
  std::vector<std::string> strings;
  std::string str;

  ezassert( createFile1(runner) );
  nc.filename = filename1(runner);
  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );

  // 4d record var with mixed slicing per dim:
  // force rec slice,index,range,stride
  strings.clear();
  strings.push_back("3");
  strings.push_back("4");
  strings.push_back("9");
  slicer.clear();
  str = "d1";
  slicer.setDimIndices(str, strings, 0, 2);
  strings.clear();
  strings.push_back("9");
  strings.push_back("10");
  strings.push_back("11");
  str = "d2";
  slicer.setDimIndices(str, strings, 0, 2);
  strings.clear();
  strings.push_back("9");
  strings.push_back("19");
  strings.push_back("29");
  str = "d3";
  slicer.setDimIndices(str, strings, 0, 2);
  slicer.build(&nc);
  var = nc.getVar("v7");
  ezassert( var );
  it.setVar(var, &slicer, true); // Force rec slice.
  it.getSlab(&slab); 
  ezassert( slab.start.size() == 4 );
  ezassert( slab.start[0] == 0 );
  ezassert( slab.start[1] == 3 );
  ezassert( slab.start[2] == 9 );
  ezassert( slab.start[3] == 9 );
  ezassert( slab.count.size() == 4 );
  ezassert( slab.count[0] == 1 ); // Should be one because forced slice.
  ezassert( slab.count[1] == 1 );
  ezassert( slab.count[2] == 3 );
  ezassert( slab.count[3] == 1 );
  ezassert( it.isDone() == false );
  while(! it.isDoneStatic() ) it.incStatic();
  it.getSlab(&slab); 
  ezassert( it.isDoneRec() == false ); 
  ezassert( it.isDone() == false );
  ezassert( slab.start[0] == 0 );
  ezassert( slab.start[1] == 9 );
  ezassert( slab.start[2] == 9 );
  ezassert( slab.start[3] == 29 );
  it.incRec();
  it.resetStatic(); // or it.incStatic(); but reset is cheaper. 
  while(! it.isDoneStatic() ) it.incStatic();
  ezassert( it.isDoneRec() == 0);
  it.incRec();
  ezassert( it.isDoneRec() ); 
  ezassert( it.isDone() ); // Should be done since only 2 records.
  it.getSlab(&slab); 
  ezassert( slab.start[0] == 1 );
  ezassert( slab.start[1] == 9 );
  ezassert( slab.start[2] == 9 );
  ezassert( slab.start[3] == 29 );

  nc.close();
  
  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );
  
  return true;
}

int main(int argc, const char* argv[]) {
  ez::ezTestRunner runner;
  runner.getopt(argc, argv);
  #define TEST(N) runner.add(&test_##N, "test_"#N);
  TEST(noslice);
  TEST(incIndexLists);
  TEST(incValueLists);
  TEST(incSeparateRecStatic);
  
  return runner.run();
}