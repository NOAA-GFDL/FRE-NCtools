/*
20111122 rsz Created.
*/

#include "ezTest.hpp"
#include "ezNc.hpp"
#include "ezNcSlabIterator.hpp"
#include "ezNcCopySlabTask.hpp"
#include "ezNcSlicer.hpp"
#include "ezNcCompare.hpp"

#define DEBUGLINE() printf("%s:%d\n", __FILE__,__LINE__);

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
  ezassert( nc.enddef() == 0 ); 

  ez::ezNcBuffers buf;
  ezassert( nc.init(buf) == 0 );
  buf.allocate();

  // Create dummy data.
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

bool test_copy(ez::ezTestRunner& runner) {
  ez::ezNc nc1, nc2;
  ez::ezNcSlabIterator<size_t> it1, it2;
  //ez::ezNcSlab<int> slab1, slab2;
  ez::ezNcCopySlabTask<size_t> task;
  ez::ezNcSlicer<size_t> slicer;
  ez::ezNcVar *var1, *var2;
  ez::ezNcBuffers buf;
  std::vector<std::string> varnames;
  
  ezassert( createFile1(runner) );
  nc1.filename = filename1(runner);
  ezassert( nc1.openReadOnly() == 0 );
  ezassert( nc1.loadAllMetadata() == 0 );
  
  std::string str1, str2, str3;
  str1 = "dr"; str2 = "0"; str3 = "1";
  slicer.setDimStride(str1,str2,str3,str3);
  slicer.build(&nc1);
  nc1.getAllVarNames(varnames);
  nc1.init(buf, varnames);
  buf.allocate();

  nc2.copySettings(nc1);
  nc2.filename = filename2(runner);
  ezassert( nc2.createClobber() == 0 );
  ezassert( nc2.copyDimDefs(nc1) == 0 );
  ezassert( nc2.copyVarDefs(nc1) == 0 );
  ezassert( nc2.enddef() == 0 );
  
  int i;
  for(i=0; i < varnames.size(); ++i) {
    var1 = nc1.getVar(varnames[i]);
    ezassert( var1 );
    var2 = nc2.getVar(varnames[i]);
    ezassert( var2 );
    // Ensure rec dim is sliced into individual records for interleave i/o.
    it1.setVar(var1, &slicer, true);
    it2 = it1;
    it2.var = var2;
    it2.normalize();
    while ( !it1.isDone() ) {
      ///it1.odo.print(','); std::cout << std::endl;
      it1.getSlab(&task.slabSrc);
      it2.getSlab(&task.slabDst);
      ezassert( task.run(buf) == 0 );
      it1.inc();
      it2.inc();
    }
  }
  
  ezassert( nc1.close() == 0 );
  ezassert( nc2.close() == 0 );

  // Check copied results.
  ez::ezNcCompare cmp;
  std::string fn1 = filename1(runner), fn2 = filename2(runner);
  cmp.init(fn1, fn2);
  for(i=0; i < varnames.size(); ++i) {
    ezassert( cmp.varDataEquals(varnames[i], varnames[i]) );
  }
  
  cmp.clear();
  
  if (runner.clean) {
    ezassert( remove(nc1.filename.c_str()) == 0 );
    ezassert( remove(nc2.filename.c_str()) == 0 );
  }

  return true;
}

int main(int argc, const char* argv[]) {
  ez::ezTestRunner runner;
  runner.getopt(argc, argv);
  #define TEST(N) runner.add(&test_##N, "test_"#N);
  TEST(copy);
  
  return runner.run();
}