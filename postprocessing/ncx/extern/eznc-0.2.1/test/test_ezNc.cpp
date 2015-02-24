/*
20111111 rsz Created.
*/

#include "ezTest.hpp"
#include "ezNc.hpp"
#include <cstdio>

#define DEBUGLINE() printf("%s:%d\n", __FILE__, __LINE__); fflush(stdout);

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
  ezassert( nc.reserve(buf, v7) == 0 );
  buf.allocate();
  buf.zeros();

  // Write last record to force other records to be written to file.
  size_t start[4], count[4];
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

bool createFile2(ez::ezTestRunner& runner) {
  ez::ezNc nc;
  
  nc.filename = filename2(runner);
  ezassert( nc.createClobber() == 0 );

  int n1=10;
  int n2=20;
  int n3=30;
  int n4=40;
  int d1 = nc.createDim("x", n1);
  ezassert(d1 != -1);
  int d2 = nc.createDim("y", n2);
  ezassert(d2 != -1);
  int d3 = nc.createDim("z", n3);
  ezassert(d3 != -1);
  int dr = nc.createDim("time", 2, true);
  ezassert(dr != -1);

  int dimidsx[1] = {d1};
  int dimidsy[1] = {d2};
  int dimidsz[1] = {d3};
  int dimids2d[2] = {d2,d1};
  int dimids3d[3] = {d3,d2,d1};
  int dimids0dt[1] = {dr};
  int dimids2dt[3] = {dr,d2,d1};
  int dimids3dt[4] = {dr,d3,d2,d1};
  int v1 = nc.createVar("x", NC_FLOAT, 1, dimidsx);
  ezassert(v1 != -1);
  int v2 = nc.createVar("y", NC_FLOAT, 1, dimidsy);
  ezassert(v2 != -1);
  int v3 = nc.createVar("z", NC_FLOAT, 1, dimidsz);
  ezassert(v3 != -1);
  int v4 = nc.createVar("ps", NC_FLOAT, 2, dimids2d);
  ezassert(v4 != -1);
  int v5 = nc.createVar("rh", NC_FLOAT, 3, dimids3d);
  ezassert(v5 != -1);
  int v6 = nc.createVar("temp", NC_DOUBLE, 3, dimids2dt);
  ezassert(v6 != -1);
  int v7 = nc.createVar("precip", NC_DOUBLE, 4, dimids3dt);
  ezassert(v7 != -1);
  int v8 = nc.createVar("time", NC_DOUBLE, 1, dimids0dt);
  ezassert(v8 != -1);

  ezassert( nc.enddef() == 0 ); 

  ez::ezNcBuffers buf;
  ezassert( nc.reserve(buf, v8) == 0 );
  buf.allocate();
  buf.zeros();

  // Write last record to force other records to be written to file.
  size_t start[4], count[4];
  start[0] = 1;
  start[1] = start[2] = start[3] = 0;
  count[0] = 1;
  count[1] = n1;
  count[2] = n2;
  count[3] = n3;
  ezassert( nc.write(v8, &buf, start, count) == 0 );
  ezassert( nc.close() == 0 );

  return true;
}

bool test_close(ez::ezTestRunner& runner) {
  ez::ezNc nc;
  
  nc.filename = runner.tmpdir + "./tmp.nc";
  ezassert( nc.createClobber() == 0 );
  ezassert( nc.close() == 0 );
  
  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );
    
  return true;
}

bool test_computeSize(ez::ezTestRunner& runner) {
  ez::ezNc nc;
  ez::ezNcVar* v;
  ezassert( createFile1(runner) );
  nc.filename = filename1(runner);

  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );

  for(int i=0; i < 7; ++i)
    nc.computeSize(i);

    v = nc.getVar((unsigned int)0);
  ezassert( v );
  ezassert( v->recsize == 10 );
  ezassert( v->recbytes == 10 );
  ezassert( v->size == 10 );
  ezassert( v->bytes == 10 );

  v = nc.getVar(1);
  ezassert( v );
  ezassert( v->recsize == 10*20 );
  ezassert( v->recbytes == 10*20 );
  ezassert( v->size == 10*20 );
  ezassert( v->bytes == 10*20 );

  v = nc.getVar(2);
  ezassert( v );
  ezassert( v->recsize == 10*20*30 );
  ezassert( v->recbytes == 2*10*20*30 );
  ezassert( v->size == 10*20*30 );
  ezassert( v->bytes == 2*10*20*30 );
 
  v = nc.getVar(3);
  ezassert( v );
  ezassert( v->recsize == 1 );
  ezassert( v->recbytes == 1*4 );
  ezassert( v->size == 2*1 );
  ezassert( v->bytes == 2*1*4 );

  v = nc.getVar(4);
  ezassert( v );
  ezassert( v->recsize == 10 );
  ezassert( v->recbytes == 10*4 );
  ezassert( v->size == 2*10 );
  ezassert( v->bytes == 2*10*4 );
  
  v = nc.getVar(5);
  ezassert( v );
  ezassert( v->recsize == 10*20 );
  ezassert( v->recbytes == 10*20*8 );
  ezassert( v->size == 2*10*20 );
  ezassert( v->bytes == 2*10*20*8 );

  v = nc.getVar(6);
  ezassert( v );
  ezassert( v->recsize == 10*20*30 );
  ezassert( v->recbytes == 10*20*30*8 );
  ezassert( v->size == 2*10*20*30 );
  ezassert( v->bytes == 2*10*20*30*8 );

  ezassert( nc.close() == 0);

  //-----------------------------------
  ezassert( createFile1(runner) );
  nc.clear();
  nc.filename = filename1(runner);

  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );

  std::string str = "v1";
  nc.computeSize(str);
  str = "v2";
  nc.computeSize(str);
  str = "v3";
  nc.computeSize(str);
  str = "v4";
  nc.computeSize(str);
  str = "v5";
  nc.computeSize(str);
  str = "v6";
  nc.computeSize(str);
  str = "v7";
  nc.computeSize(str);

  v = nc.getVar((unsigned int)0);
  ezassert( v );
  ezassert( v->recsize == 10 );
  ezassert( v->recbytes == 10 );
  ezassert( v->size == 10 );
  ezassert( v->bytes == 10 );

  v = nc.getVar(1);
  ezassert( v );
  ezassert( v->recsize == 10*20 );
  ezassert( v->recbytes == 10*20 );
  ezassert( v->size == 10*20 );
  ezassert( v->bytes == 10*20 );

  v = nc.getVar(2);
  ezassert( v );
  ezassert( v->recsize == 10*20*30 );
  ezassert( v->recbytes == 2*10*20*30 );
  ezassert( v->size == 10*20*30 );
  ezassert( v->bytes == 2*10*20*30 );
 
  v = nc.getVar(3);
  ezassert( v );
  ezassert( v->recsize == 1 );
  ezassert( v->recbytes == 1*4 );
  ezassert( v->size == 2*1 );
  ezassert( v->bytes == 2*1*4 );

  v = nc.getVar(4);
  ezassert( v );
  ezassert( v->recsize == 10 );
  ezassert( v->recbytes == 10*4 );
  ezassert( v->size == 2*10 );
  ezassert( v->bytes == 2*10*4 );
  
  v = nc.getVar(5);
  ezassert( v );
  ezassert( v->recsize == 10*20 );
  ezassert( v->recbytes == 10*20*8 );
  ezassert( v->size == 2*10*20 );
  ezassert( v->bytes == 2*10*20*8 );

  v = nc.getVar(6);
  ezassert( v );
  ezassert( v->recsize == 10*20*30 );
  ezassert( v->recbytes == 10*20*30*8 );
  ezassert( v->size == 2*10*20*30 );
  ezassert( v->bytes == 2*10*20*30*8 );

  ezassert( nc.close() == 0);

  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );
  
  return true;
}

bool test_copyAtt(ez::ezTestRunner& runner) {
  ez::ezNc nc1, nc2;
  
  ezassert( createFile1(runner) );
  nc1.filename = filename1(runner);
  ezassert( nc1.openReadOnly() == 0 );
  ezassert( nc1.loadAllMetadata() == 0 );

  nc2.filename = runner.tmpdir + "./tmp2.nc";
  ezassert( nc2.createClobber() == 0 );
  ezassert( nc2.copyAtts(nc1, 0, 0) == 0 );
  ezassert( nc2.enddef() == 0 ); 

  ez::ezNcBuffers atts;
  ezassert( nc2.getAtt(0, "Author", atts) == 0 );
  ezassert( strncmp("Quick Fox", atts.c, 9) == 0 );
  ezassert( nc2.getAtt(0, "age", atts) == 0 );
  ezassert( *atts.i == 2123000 );
  ezassert( nc2.getAtt(0, "microsecond", atts) == 0 );
  ezassert( *atts.f == 1e-6f );
  ezassert( nc2.getAtt(0, "nanosecond", atts) == 0 );
  ezassert( *atts.d == 1e-9 );
  
  ez::ezNc nc3;
  nc3.filename = runner.tmpdir + "./tmp3.nc";
  ezassert( nc3.createClobber() == 0 );
  int d1 = nc3.createDim("d1", 10);
  ezassert(d1 != -1);
  int dimids1d[1] = {d1};
  int v1 = nc3.createVar("v1", NC_BYTE, 1, dimids1d);
  ezassert(v1 != -1);
  int v2 = nc3.createVar("v2", NC_INT, 1, dimids1d);
  ezassert(v2 != -1);

  ezassert( nc3.copyAtts(nc1, "v3", "v1") == 0 );
  ezassert( nc3.copyAtts(nc1, "v4", "v2") == 0 );
  ezassert( nc3.enddef() == 0 ); 

  ezassert( nc3.getAtt("v1", "offset", atts) == 0 );
  ezassert( *atts.s == 32000 );
  ezassert( nc3.getAtt("v2", "units", atts) == 0 );
  ezassert( strncmp("days since 1980-01-01 12:30:45", atts.c, 30) == 0 );
  ezassert( nc3.getAtt("v2", "missing_value", atts) == 0 );
  ezassert( *atts.i == -2123000 );
  ezassert( nc3.getAtt("v2", "pi", atts) == 0 );
  ezassert( *atts.f == 3.14f );
  ezassert( nc3.getAtt("v2", "double_array", atts) == 0 );
  ezassert( atts.d[0] == 1.1 );
  ezassert( atts.d[1] == 2.2 );
  ezassert( atts.d[2] == 3.3 );

  ezassert( nc1.close() == 0);
  ezassert( nc2.close() == 0);
  ezassert( nc3.close() == 0);

  if (runner.clean) {
    ezassert( remove(nc1.filename.c_str()) == 0 );
    ezassert( remove(nc2.filename.c_str()) == 0 );
    ezassert( remove(nc3.filename.c_str()) == 0 );
  }
  
  return true;
}

bool test_copyDimDef(ez::ezTestRunner& runner) {
  ez::ezNc nc1, nc2;
  ez::ezNcDim * dim;
  
  ezassert( createFile1(runner) );
  nc1.filename = filename1(runner);
  ezassert( nc1.openReadOnly() == 0 );
  ezassert( nc1.loadAllMetadata() == 0 );

  nc2.filename = runner.tmpdir + "./tmp2.nc";
  ezassert( nc2.createClobber() == 0 );
  ezassert( nc2.copyDimDef(nc1, "d1") == 0 );
  ezassert( nc2.copyDimDef(nc1, "d2") == 0 );
  ezassert( nc2.copyDimDef(nc1, "d3") == 0 );
  ezassert( nc2.copyDimDef(nc1, 3) == 0 );
  ezassert( nc2.copyDimDef(nc1, 4) == 0 );
  ezassert( nc2.enddef() == 0 ); 

  dim = nc2.getDim((unsigned int)0);
  ezassert( dim && (dim->size == 10) );
  dim = nc2.getDim(1);
  ezassert( dim && (dim->size == 20) );
  dim = nc2.getDim(2);
  ezassert( dim && (dim->size == 30) );
  dim = nc2.getDim(3);
  ezassert( dim && (dim->size == 40) );
  dim = nc2.getDim(4);
  ezassert( dim && (dim->size == 2) );
 
  ezassert( nc1.close() == 0);
  ezassert( nc2.close() == 0);

  if (runner.clean) {
    ezassert( remove(nc1.filename.c_str()) == 0 );
    ezassert( remove(nc2.filename.c_str()) == 0 );
  }
  
  return true;
}

bool test_copyDimDefs(ez::ezTestRunner& runner) {
  ez::ezNc nc1, nc2;
  ez::ezNcDim * dim;
  
  ezassert( createFile1(runner) );
  nc1.filename = filename1(runner);
  ezassert( nc1.openReadOnly() == 0 );
  ezassert( nc1.loadAllMetadata() == 0 );

  nc2.filename = runner.tmpdir + "./tmp2.nc";
  ezassert( nc2.createClobber() == 0 );
  ezassert( nc2.copyDimDefs(nc1) == 0 );
  ezassert( nc2.enddef() == 0 ); 
  dim = nc2.getDim((unsigned int)0);
  ezassert( dim && (dim->size == 10) );
  dim = nc2.getDim(1);
  ezassert( dim && (dim->size == 20) );
  dim = nc2.getDim(2);
  ezassert( dim && (dim->size == 30) );
  dim = nc2.getDim(3);
  ezassert( dim && (dim->size == 40) );
  dim = nc2.getDim(4);

  ezassert( dim && (dim->size == 2) );
  ezassert( nc1.close() == 0);
  ezassert( nc2.close() == 0);

  if (runner.clean) {
    ezassert( remove(nc1.filename.c_str()) == 0 );
    ezassert( remove(nc2.filename.c_str()) == 0 );
  }
  
  return true;
}

bool test_copySettings(ez::ezTestRunner& runner) {
  ez::ezNc nc1, nc2;
  ez::ezNcDim * dim;
  
  ezassert( createFile1(runner) );
  nc1.filename = filename1(runner);
  ezassert( nc1.openReadOnly() == 0 );
  ezassert( nc1.loadAllMetadata() == 0 );

  std::string filename2 = runner.tmpdir + "./tmp2.nc";
  nc2.filename = filename2;
  ezassert( nc2.createClobber() == 0 );
  int oldncid = nc2.ncid;
  nc2.copySettings(nc1);
  ezassert( nc2.filename.compare(nc1.filename) == 0 );
  nc2.filename = filename2;
  #define SEQUALS(V) ezassert( nc2.V == nc1.V );
  SEQUALS(ncid);
  SEQUALS(format);
  SEQUALS(recid);
  SEQUALS(nrec);
  SEQUALS(initialsz);
  SEQUALS(bufrsizehint);
  SEQUALS(h_minfree);
  SEQUALS(v_align);
  SEQUALS(v_minfree);
  SEQUALS(r_align);
  SEQUALS(fill);
  SEQUALS(chunksize);
  SEQUALS(chunknelems);
  SEQUALS(chunkpreemption);
  nc2.ncid = oldncid;
  ezassert( nc2.enddef() == 0 ); 

  ezassert( nc1.close() == 0);
  ezassert( nc2.close() == 0);

  if (runner.clean) {
    ezassert( remove(nc1.filename.c_str()) == 0 );
    ezassert( remove(nc2.filename.c_str()) == 0 );
  }
  
  return true;
}

bool test_copyVarDef(ez::ezTestRunner& runner) {
  ez::ezNc nc1, nc2;
  ez::ezNcVar *v1,*v2;
  
  ezassert( createFile1(runner) );
  nc1.filename = filename1(runner);
  ezassert( nc1.openReadOnly() == 0 );
  ezassert( nc1.loadAllMetadata() == 0 );

  nc2.filename = runner.tmpdir + "./tmp2.nc";
  nc2.nrec = nc1.nrec;
  ezassert( nc2.createClobber() == 0 );
  ezassert( nc2.copyDimDefs(nc1) == 0 );
  
  ezassert( nc2.copyVarDef(nc1, "v1") != -1 );
  ezassert( nc2.copyVarDef(nc1, "v2") != -1 );
  ezassert( nc2.copyVarDef(nc1, "v3") != -1 );
  ezassert( nc2.copyVarDef(nc1, "v4") != -1 );
  ezassert( nc2.copyVarDef(nc1, "v5") != -1 );
  ezassert( nc2.copyVarDef(nc1, "v6") != -1 );
  ezassert( nc2.copyVarDef(nc1, "v7") != -1 );
  ezassert( nc2.enddef() == 0 ); 
  
  v1 = nc1.getVar("v4");
  v2 = nc2.getVar("v4");
  ezassert( v1 );
  ezassert( v2 );
  ezassert( v1->name.compare(v2->name) == 0 );
  ezassert( v1->getNumDims() == v2->getNumDims() );
  ezassert( v1->bytes == v2->bytes );
  ezassert( v1->type == v2->type );
  ezassert( v1->hasrec == v2->hasrec );
  ezassert( v1->hasrec == v2->hasrec );
  ezassert( std::equal( v1->attnames.begin(), v1->attnames.end(), v2->attnames.begin() ) );
  ezassert( v1->chunkstorage == v2->chunkstorage );
  ezassert( std::equal( v1->chunksizes.begin(), v1->chunksizes.end(), v2->chunksizes.begin() ) );
  ezassert( v1->chunksize == v2->chunksize );
  ezassert( v1->chunknelems == v2->chunknelems );
  ezassert( v1->chunkpreemption == v2->chunkpreemption );
  ezassert( v1->no_fill == v2->no_fill );
  ezassert( v1->lfill == v2->lfill );
  ezassert( v1->shuffle == v2->shuffle );
  ezassert( v1->deflate == v2->deflate );
  ezassert( v1->deflate_level == v2->deflate_level );
  ezassert( v1->checksum == v2->checksum );
  ezassert( v1->endian == v2->endian );
  
  ezassert( nc1.close() == 0);
  ezassert( nc2.close() == 0);

  if (runner.clean) {
    ezassert( remove(nc1.filename.c_str()) == 0 );
    ezassert( remove(nc2.filename.c_str()) == 0 );
  }
  
  return true;
}

bool test_copyVarDefs(ez::ezTestRunner& runner) {
  ez::ezNc nc1, nc2;
  ez::ezNcVar *v1,*v2;
  
  ezassert( createFile1(runner) );
  nc1.filename = filename1(runner);
  ezassert( nc1.openReadOnly() == 0 );
  ezassert( nc1.loadAllMetadata() == 0 );

  nc2.filename = runner.tmpdir + "./tmp2.nc";
  nc2.nrec = nc1.nrec;
  ezassert( nc2.createClobber() == 0 );
  ezassert( nc2.copyDimDefs(nc1) == 0 );
  
  ezassert( nc2.copyVarDefs(nc1) == 0 );
  ezassert( nc2.enddef() == 0 ); 
  
  v1 = nc1.getVar("v4");
  v2 = nc2.getVar("v4");
  ezassert( v1 );
  ezassert( v2 );
  ezassert( v1->name.compare(v2->name) == 0 );
  ezassert( v1->getNumDims() == v2->getNumDims() );
  ezassert( v1->bytes == v2->bytes );
  ezassert( v1->type == v2->type );
  ezassert( v1->hasrec == v2->hasrec );
  ezassert( v1->hasrec == v2->hasrec );
  ezassert( std::equal( v1->attnames.begin(), v1->attnames.end(), v2->attnames.begin() ) );
  ezassert( v1->chunkstorage == v2->chunkstorage );
  ezassert( std::equal( v1->chunksizes.begin(), v1->chunksizes.end(), v2->chunksizes.begin() ) );
  ezassert( v1->chunksize == v2->chunksize );
  ezassert( v1->chunknelems == v2->chunknelems );
  ezassert( v1->chunkpreemption == v2->chunkpreemption );
  ezassert( v1->no_fill == v2->no_fill );
  ezassert( v1->lfill == v2->lfill );
  ezassert( v1->shuffle == v2->shuffle );
  ezassert( v1->deflate == v2->deflate );
  ezassert( v1->deflate_level == v2->deflate_level );
  ezassert( v1->checksum == v2->checksum );
  ezassert( v1->endian == v2->endian );
  
  ezassert( nc1.close() == 0);
  ezassert( nc2.close() == 0);

  if (runner.clean) {
    ezassert( remove(nc1.filename.c_str()) == 0 );
    ezassert( remove(nc2.filename.c_str()) == 0 );
  }
  
  return true;
}

bool test_create64(ez::ezTestRunner& runner) {
  ez::ezNc nc1,nc2;
  
  nc1.filename = runner.tmpdir + "./tmp1.nc";
  ezassert( nc1.createClobber() == 0 );
  ezassert( nc1.enddef() == 0 ); 
  ezassert( nc1.close() == 0);

  nc1.clear();
  nc1.filename = runner.tmpdir + "./tmp1.nc";
  ezassert( nc1.create64() != 0 ); // Should fail to overwrite.

  nc2.filename = runner.tmpdir + "./tmp2.nc";
  ezassert( nc2.create64() == 0);
  ezassert( nc2.enddef() == 0 ); 
  ezassert( nc2.close() == 0);

  if (runner.clean) {
    ezassert( remove(nc1.filename.c_str()) == 0 );
    ezassert( remove(nc2.filename.c_str()) == 0 );
  }
  
  return true;
}

bool test_create64Clobber(ez::ezTestRunner& runner) {
  ez::ezNc nc1,nc2;
  
  nc1.filename = runner.tmpdir + "./tmp1.nc";
  ezassert( nc1.createClobber() == 0 );
  ezassert( nc1.enddef() == 0 ); 
  ezassert( nc1.close() == 0);

  nc1.clear();
  nc1.filename = runner.tmpdir + "./tmp1.nc";
  ezassert( nc1.create64Clobber() == 0 );
  ezassert( nc1.enddef() == 0 ); 
  ezassert( nc1.close() == 0);

  if (runner.clean) {
    ezassert( remove(nc1.filename.c_str()) == 0 );
  }
  
  return true;
}

bool test_createHDF(ez::ezTestRunner& runner) {
  ez::ezNc nc1,nc2;
  
  nc1.filename = runner.tmpdir + "./tmp1.nc";
  ezassert( nc1.createClobber() == 0 );
  ezassert( nc1.enddef() == 0 ); 
  ezassert( nc1.close() == 0);

  nc1.clear();
  nc1.filename = runner.tmpdir + "./tmp1.nc";
  ezassert( nc1.createHDF() != 0 ); // Should fail to overwrite.

  nc2.filename = runner.tmpdir + "./tmp2.nc";
  ezassert( nc2.createHDF() == 0);
  ezassert( nc2.enddef() == 0 ); 
  ezassert( nc2.close() == 0);

  if (runner.clean) {
    ezassert( remove(nc1.filename.c_str()) == 0 );
    ezassert( remove(nc2.filename.c_str()) == 0 );
  }
  
  return true;
}

bool test_createHDFClobber(ez::ezTestRunner& runner) {
  ez::ezNc nc1,nc2;
  
  nc1.filename = runner.tmpdir + "./tmp1.nc";
  ezassert( nc1.createClobber() == 0 );
  ezassert( nc1.enddef() == 0 ); 
  ezassert( nc1.close() == 0);

  nc1.clear();
  nc1.filename = runner.tmpdir + "./tmp1.nc";
  ezassert( nc1.createHDFClobber() == 0 );
  ezassert( nc1.enddef() == 0 ); 
  ezassert( nc1.close() == 0);

  if (runner.clean) {
    ezassert( remove(nc1.filename.c_str()) == 0 );
  }
  
  return true;
}

bool test_createNcDim(ez::ezTestRunner& runner) {
  // Non-null creation already done in other unit tests.
  ez::ezNc nc;
  ez::ezNcDim* d = 0;

  ezassert( nc.create(d) == -1 ); // Fails with null nc.
  
  nc.filename = filename1(runner);
  ezassert( nc.createClobber() == 0 );
  
  ezassert( nc.create(d) == -1 ); // Fails with null dim.
  ezassert( nc.enddef() == 0 ); 
  ezassert( nc.close() == 0);

  if (runner.clean) {
    ezassert( remove(nc.filename.c_str()) == 0 );
  }
  
  return true;
}

bool test_createDim(ez::ezTestRunner& runner) {
  // Non-null creation already done in other unit tests.
  ez::ezNc nc;

  ezassert( nc.createDim(0,0) == -1 ); // Fails with null nc.
  
  nc.filename = filename1(runner);
  ezassert( nc.createClobber() == 0 );
  
  ezassert( nc.createDim(0,0) == -1 ); // Fails with null dim.
  ezassert( nc.enddef() == 0 ); 
  ezassert( nc.close() == 0);

  if (runner.clean) {
    ezassert( remove(nc.filename.c_str()) == 0 );
  }
  
  return true;
}

bool test_createNcVar(ez::ezTestRunner& runner) {
  // Non-null creation already done in other unit tests.
  ez::ezNc nc;
  ez::ezNcVar* v = 0;

  ezassert( nc.create(v) == -1 ); // Fails with null nc.
  
  nc.filename = filename1(runner);
  ezassert( nc.createClobber() == 0 );
  
  ezassert( nc.create(v) == -1 ); // Fails with null.
  ezassert( nc.enddef() == 0 ); 
  ezassert( nc.close() == 0);

  if (runner.clean) {
    ezassert( remove(nc.filename.c_str()) == 0 );
  }
  
  return true;
}

bool test_createVar(ez::ezTestRunner& runner) {
  // Non-null creation already done in other unit tests.
  ez::ezNc nc;

  ezassert( nc.createVar(0,0,0,0) == -1 ); // Fails with null nc.
  
  nc.filename = filename1(runner);
  ezassert( nc.createClobber() == 0 );
  
  ezassert( nc.createVar(0,0,0,0) == -1 ); // Fails with null var.
  ezassert( nc.enddef() == 0 ); 
  ezassert( nc.close() == 0);

  if (runner.clean) {
    ezassert( remove(nc.filename.c_str()) == 0 );
  }
  
  return true;
}

bool test_formatToString(ez::ezTestRunner& runner) {
  std::string s;
  
  ez::formatToString(NC_FORMAT_CLASSIC, s);
  ezassert( s.compare("NC_FORMAT_CLASSIC") == 0 );
  
  ez::formatToString(NC_FORMAT_64BIT, s);
  ezassert( s.compare("NC_FORMAT_64BIT") == 0 );
  
  ez::formatToString(NC_FORMAT_NETCDF4, s);
  ezassert( s.compare("NC_FORMAT_NETCDF4") == 0 );
  
  ez::formatToString(NC_FORMAT_NETCDF4_CLASSIC, s);
  ezassert( s.compare("NC_FORMAT_NETCDF4_CLASSIC") == 0 );

  return true;
}

bool test_getAllVarBytes(ez::ezTestRunner& runner) {
  ez::ezNc nc;
  size_t nb;
  
  nb = nc.getAllVarBytes();
  ezassert( nb == 0 );
  
  ezassert( createFile1(runner) );
  nc.filename = filename1(runner);

  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );
  
  nb = nc.getAllVarBytes();
  ezassert( nb == (10 + 10*20 + 2*10*20*30 + 2*1*4 + 2*10*4 + 2*10*20*8 + 2*10*20*30*8) );
  ezassert( nc.close() == 0);

  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );
  
  return true;  
}

bool test_getAllVarIds(ez::ezTestRunner& runner) {
  ez::ezNc nc;
  ezassert( createFile1(runner) );
  nc.filename = filename1(runner);

  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );
  
  std::vector<int> ids;
  nc.getAllVarIds(ids);
  ezassert( ids.size() == 7 );
  ezassert( ids[0] == 0 );
  ezassert( ids[6] == 6 );
  
  ezassert( nc.close() == 0);

  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );
  
  return true;    
}

bool test_getAllVarNames(ez::ezTestRunner& runner) {
  std::vector<std::string> names;
  ez::ezNc nc;
  ezassert( createFile1(runner) );
  nc.filename = filename1(runner);

  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );
  
  nc.getAllVarNames(names);
  ezassert( names.size() == 7 );
  ezassert( names[0].compare("v1") == 0 );
  ezassert( names[6].compare("v7") == 0 );
  
  ezassert( nc.close() == 0);

  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );
  
  return true;  
}

bool test_getDim(ez::ezTestRunner& runner) {
  ez::ezNc nc;
  ezassert( createFile1(runner) );
  nc.filename = filename1(runner);

  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );
  
  ez::ezNcDim* d;
  d = nc.getDim((unsigned int)0);
  ezassert( d );
  
  d = nc.getDim((char*)0);
  ezassert( d == 0 );

  d = nc.getDim("");
  ezassert( d == 0 );

  std::string s;
  d = nc.getDim(s);
  ezassert( d == 0 );

  d = nc.getDim(10);
  ezassert( d == 0);

  d = nc.getDim("d4");
  ezassert( d );

  s = "dr";
  ezassert( d );
  
  ezassert( nc.close() == 0);

  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );
  
  return true;
}

bool test_getDimIds(ez::ezTestRunner& runner) {
  std::vector<std::string> names;
  std::vector<int> ids;
  ez::ezNc nc;
  ezassert( createFile1(runner) );
  nc.filename = filename1(runner);

  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );
  
  names.push_back("d2");
  names.push_back("d1");
  names.push_back("d3");
  names.push_back("dr");
  names.push_back("dmissing");
  names.push_back("d4");
  nc.getDimIds(names, ids);
  
  ezassert( ids.size() == 6 );
  ezassert( ids[0] == 1 );
  ezassert( ids[1] == 0 );
  ezassert( ids[2] == 2 );
  ezassert( ids[3] == 4 );
  ezassert( ids[4] < 0 );
  ezassert( ids[5] == 3 );

  ezassert( nc.close() == 0);

  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );
  
  return true;
}

bool test_getDimSizesId(ez::ezTestRunner& runner) {
  std::vector<std::string> names;
  std::vector<int> ids;
  ez::ezNc nc;
  ezassert( createFile1(runner) );
  nc.filename = filename1(runner);

  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );
  
  std::vector<size_t> sizes;
  
  nc.getDimSizes(999, sizes);
  ezassert( sizes.empty() );
  
  nc.getDimSizes(6, sizes);
  ezassert( sizes.size() == 4 );
  ezassert( sizes[0] == 2 );
  ezassert( sizes[1] == 10 );
  ezassert( sizes[2] == 20 );
  ezassert( sizes[3] == 30 );
  
  ezassert( nc.close() == 0);

  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );
  
  return true;  
}

bool test_getDimSizesString(ez::ezTestRunner& runner) {
  ez::ezNc nc;
  ezassert( createFile1(runner) );
  nc.filename = filename1(runner);

  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );
  
  std::vector<size_t> sizes;
  
  std::string s;
  nc.getDimSizes(s, sizes);
  ezassert( sizes.empty() );

  s = "missing";
  nc.getDimSizes(s, sizes);
  ezassert( sizes.empty() );
  
  s = "v7";
  nc.getDimSizes(s, sizes);
  ezassert( sizes.size() == 4 );
  ezassert( sizes[0] == 2 );
  ezassert( sizes[1] == 10 );
  ezassert( sizes[2] == 20 );
  ezassert( sizes[3] == 30 );
  
  ezassert( nc.close() == 0);

  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );
  
  return true;  
}

bool test_getDimNamesForVars(ez::ezTestRunner& runner) {
  std::vector<std::string> vars;
  std::set<std::string> dims;
  ez::ezNc nc;
  ezassert( createFile1(runner) );
  nc.filename = filename1(runner);

  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );
   
  nc.getDimNamesForVars(vars, dims);
  ezassert( dims.empty() );

  vars.push_back("v1");
  vars.push_back("v2");
  vars.push_back("v3");
  vars.push_back("v4");
  vars.push_back("v5");
  vars.push_back("v6");
  vars.push_back("v7");
  nc.getDimNamesForVars(vars, dims);
  ezassert( dims.size() == 4 );
  ezassert( dims.count("d1") );
  ezassert( dims.count("d2") );
  ezassert( dims.count("d3") );
  ezassert( dims.count("dr") );
  
  ezassert( nc.close() == 0);

  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );
  
  return true;  
}

bool test_getMaxNumDims(ez::ezTestRunner& runner) {
  std::vector<std::string> vars;
  ez::ezNc nc;
  ezassert( createFile1(runner) );
  nc.filename = filename1(runner);

  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );
   
  ezassert( nc.getMaxNumDims(vars) == 0 );
  vars.push_back("v1");
  ezassert( nc.getMaxNumDims(vars) == 1 );
  vars.push_back("v2");
  ezassert( nc.getMaxNumDims(vars) == 2 );
  vars.push_back("v3");
  ezassert( nc.getMaxNumDims(vars) == 3 );
  vars.push_back("v4");
  ezassert( nc.getMaxNumDims(vars) == 3 );
  vars.push_back("v5");
  ezassert( nc.getMaxNumDims(vars) == 3 );
  vars.push_back("v6");
  ezassert( nc.getMaxNumDims(vars) == 3 );
  vars.push_back("v7");
  ezassert( nc.getMaxNumDims(vars) == 4 );
  
  ezassert( nc.close() == 0);

  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );
  
  return true;  
}

bool test_getNumVars(ez::ezTestRunner& runner) {
  ez::ezNc nc;
  ezassert( createFile1(runner) );
  nc.filename = filename1(runner);

  ezassert( nc.getNumVars() == 0 );
  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.getNumVars() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );
  ezassert( nc.getNumVars() == 7 );
  
  ezassert( nc.close() == 0);

  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );
  
  return true;     
}

bool test_getNumStaticVars(ez::ezTestRunner& runner) {
  ez::ezNc nc;
  ezassert( createFile1(runner) );
  nc.filename = filename1(runner);

  ezassert( nc.getNumStaticVars() == 0 );
  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.getNumStaticVars() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );
  ezassert( nc.getNumStaticVars() == 3 );
  
  ezassert( nc.close() == 0);

  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );
  
  return true;     
}

bool test_getNumRecVars(ez::ezTestRunner& runner) {
  ez::ezNc nc;
  ezassert( createFile1(runner) );
  nc.filename = filename1(runner);

  ezassert( nc.getNumRecVars() == 0 );
  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.getNumRecVars() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );
  ezassert( nc.getNumRecVars() == 4 );
  
  ezassert( nc.close() == 0);

  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );
  
  return true;     
}

bool test_getRecVarBytes(ez::ezTestRunner& runner) {
  ez::ezNc nc;
  size_t nb;
  
  nb = nc.getAllVarBytes();
  ezassert( nb == 0 );
  
  ezassert( createFile1(runner) );
  nc.filename = filename1(runner);

  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );
  
  nb = nc.getRecVarBytes();
  ezassert( nb == (2*1*4 + 2*10*4 + 2*10*20*8 + 2*10*20*30*8) );
  ezassert( nc.close() == 0);

  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );
  
  return true;  
}

bool test_getVar(ez::ezTestRunner& runner) {
  ez::ezNc nc;
  char* c = 0;
  std::string s;
  unsigned int i = 0;
  
  ezassert( nc.getVar(i) == 0 );
  ezassert( nc.getVar(c) == 0 );
  ezassert( nc.getVar("") == 0 );
  ezassert( nc.getVar(s) == 0 );
  
  ezassert( createFile1(runner) );
  nc.filename = filename1(runner);

  ezassert( nc.openReadOnly() == 0 );

  ezassert( nc.getVar(i) == 0 );
  ezassert( nc.getVar(c) == 0 );
  ezassert( nc.getVar("") == 0 );
  ezassert( nc.getVar(s) == 0 );

  ezassert( nc.loadAllMetadata() == 0 );
  ezassert( nc.getVar(999) == 0 );
  ezassert( nc.getVar(c) == 0 );
  ezassert( nc.getVar("") == 0 );
  ezassert( nc.getVar(s) == 0 );

  i = 2;
  ezassert( nc.getVar(i) );
  ezassert( nc.getVar(5) );
  ezassert( nc.getVar("v7") );
  s = "v4";
  ezassert( nc.getVar(s) );

  ezassert( nc.close() == 0);

  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );
  
  return true;
}

bool test_getRecVarNames(ez::ezTestRunner& runner) {
  ez::ezNc nc;
  std::vector<std::string> strings;

  nc.getRecVarNames(strings);
  ezassert( strings.empty() );

  ezassert( createFile1(runner) );
  nc.filename = filename1(runner);

  ezassert( nc.openReadOnly() == 0 );
  nc.getRecVarNames(strings);
  ezassert( strings.empty() );
  
  ezassert( nc.loadAllMetadata() == 0 );

  nc.getRecVarNames(strings);
  ezassert( strings.size() == 4 );
  ezassert( strings[0].compare("v4") == 0 );
  ezassert( strings[1].compare("v5") == 0 );
  ezassert( strings[2].compare("v6") == 0 );
  ezassert( strings[3].compare("v7") == 0 );

  ezassert( nc.close() == 0);

  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );
  
  return true;
}

bool test_getRecVarIds(ez::ezTestRunner& runner) {
  ez::ezNc nc;
  std::vector<int> ids;

  nc.getRecVarIds(ids);
  ezassert( ids.empty() );

  ezassert( createFile1(runner) );
  nc.filename = filename1(runner);

  ezassert( nc.openReadOnly() == 0 );
  nc.getRecVarIds(ids);
  ezassert( ids.empty() );
  
  ezassert( nc.loadAllMetadata() == 0 );

  nc.getRecVarIds(ids);
  ezassert( ids.size() == 4 );
  ezassert( ids[0] == 3 );
  ezassert( ids[1] == 4 );
  ezassert( ids[2] == 5 );
  ezassert( ids[3] == 6 );

  ezassert( nc.close() == 0);

  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );
  
  return true;
}

bool test_getStaticVarBytes(ez::ezTestRunner& runner) {
  ez::ezNc nc;
  size_t nb;
  
  nb = nc.getStaticVarBytes();
  ezassert( nb == 0 );
  
  ezassert( createFile1(runner) );
  nc.filename = filename1(runner);

  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );
  
  nb = nc.getStaticVarBytes();
  ezassert( nb == (10 + 10*20 + 2*10*20*30) );
  ezassert( nc.close() == 0);

  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );
  
  return true;  
}

bool test_getStaticVarNames(ez::ezTestRunner& runner) {
  ez::ezNc nc;
  std::vector<std::string> strings;

  nc.getStaticVarNames(strings);
  ezassert( strings.empty() );

  ezassert( createFile1(runner) );
  nc.filename = filename1(runner);

  ezassert( nc.openReadOnly() == 0 );
  nc.getStaticVarNames(strings);
  ezassert( strings.empty() );
  
  ezassert( nc.loadAllMetadata() == 0 );

  nc.getStaticVarNames(strings);
  ezassert( strings.size() == 3 );
  ezassert( strings[0].compare("v1") == 0 );
  ezassert( strings[1].compare("v2") == 0 );
  ezassert( strings[2].compare("v3") == 0 );

  ezassert( nc.close() == 0);

  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );
  
  return true;
}

bool test_getStaticVarIds(ez::ezTestRunner& runner) {
  ez::ezNc nc;
  std::vector<int> ids;

  nc.getStaticVarIds(ids);
  ezassert( ids.empty() );

  ezassert( createFile1(runner) );
  nc.filename = filename1(runner);

  ezassert( nc.openReadOnly() == 0 );
  nc.getStaticVarIds(ids);
  ezassert( ids.empty() );
  
  ezassert( nc.loadAllMetadata() == 0 );

  nc.getStaticVarIds(ids);
  ezassert( ids.size() == 3 );
  ezassert( ids[0] == 0 );
  ezassert( ids[1] == 1 );
  ezassert( ids[2] == 2 );

  ezassert( nc.close() == 0);

  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );
  
  return true;
}

bool test_getVarSizeTZYX(ez::ezTestRunner& runner) {
  ez::ezNc nc;
  size_t nt, nz, ny, nx;
  
  nc.getVarSizeTZYX(0,nt,nz,ny,nx);
  ezassert( nt == 0 );
  ezassert( nz == 0 );
  ezassert( ny == 0 );
  ezassert( nx == 0 );
  
  ezassert( createFile1(runner) );
  nc.filename = filename1(runner);

  ezassert( nc.openReadOnly() == 0 );
  nc.getVarSizeTZYX(0,nt,nz,ny,nx);
  ezassert( nt == 0 );
  ezassert( nz == 0 );
  ezassert( ny == 0 );
  ezassert( nx == 0 );
  
  ezassert( nc.loadAllMetadata() == 0 );
  
  nc.getVarSizeTZYX(0,nt,nz,ny,nx);
  ezassert( nt == 0 );
  ezassert( nz == 0 );
  ezassert( ny == 0 );
  ezassert( nx == 10 );

  nc.getVarSizeTZYX(1,nt,nz,ny,nx);
  ezassert( nt == 0 );
  ezassert( nz == 0 );
  ezassert( ny == 10 );
  ezassert( nx == 20 );
  
  std::string s = "v3";
  nc.getVarSizeTZYX(s,nt,nz,ny,nx);
  ezassert( nt == 0 );
  ezassert( nz == 10 );
  ezassert( ny == 20 );
  ezassert( nx == 30 );

  s = "v4";
  nc.getVarSizeTZYX(s,nt,nz,ny,nx);
  ezassert( nt == 2 );
  ezassert( nz == 0 );
  ezassert( ny == 0 );
  ezassert( nx == 0 );

  s = "v5";
  nc.getVarSizeTZYX(s,nt,nz,ny,nx);
  ezassert( nt == 2 );
  ezassert( nz == 0 );
  ezassert( ny == 0 );
  ezassert( nx == 10 );
  
  s = "v6";
  nc.getVarSizeTZYX(s,nt,nz,ny,nx);
  ezassert( nt == 2 );
  ezassert( nz == 0 );
  ezassert( ny == 10 );
  ezassert( nx == 20 );

  nc.getVarSizeTZYX(6,nt,nz,ny,nx);
  ezassert( nt == 2 );
  ezassert( nz == 10 );
  ezassert( ny == 20 );
  ezassert( nx == 30 );
  
  ezassert( nc.close() == 0);

  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );
  
  return true;
}

bool test_getVars(ez::ezTestRunner& runner) {
  ez::ezNc nc;
  std::vector<ez::ezNcVar*> vars;
  std::vector<std::string> patterns;
  bool regex, regex_nocase, exclude, coordvars;
  
  //-------------------------------------
  // Base case.
  regex = 0;
  regex_nocase = 0;
  exclude = 0;
  coordvars = 0;
  nc.getVars(vars, patterns, regex, regex_nocase, coordvars, exclude);
  ezassert( vars.empty() );
  
  nc.filename = filename2(runner);
  ezassert( createFile2(runner) );
  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );

  //-------------------------------------
  vars.clear();
  patterns.push_back("temp");
  patterns.push_back("rh");
  nc.getVars(vars, patterns, regex, regex_nocase, coordvars, exclude);
  ezassert( vars.size() == 2 );
  ezassert( vars[0]->name.compare("rh") == 0 );
  ezassert( vars[1]->name.compare("temp") == 0 );
  
  //-------------------------------------
  vars.clear();
  patterns.clear();
  patterns.push_back("^[pr]");
  regex = 1;
  nc.getVars(vars, patterns, regex, regex_nocase, coordvars, exclude);
  ezassert( vars.size() == 3 );
  ezassert( vars[0]->name.compare("ps") == 0 );
  ezassert( vars[1]->name.compare("rh") == 0 );
  ezassert( vars[2]->name.compare("precip") == 0 );

  //-------------------------------------
  vars.clear();
  patterns.clear();
  patterns.push_back("T");
  regex = 1;
  regex_nocase = 1;
  nc.getVars(vars, patterns, regex, regex_nocase, coordvars, exclude);
  ezassert( vars.size() == 2 );
  ezassert( vars[0]->name.compare("temp") == 0 );
  ezassert( vars[1]->name.compare("time") == 0 );

  //-------------------------------------
  vars.clear();
  patterns.clear();
  patterns.push_back("T");
  regex = 1;
  regex_nocase = 1;
  coordvars = 1;
  nc.getVars(vars, patterns, regex, regex_nocase, coordvars, exclude);
  ezassert( vars.size() == 4 );
  ezassert( vars[0]->name.compare("x") == 0 );
  ezassert( vars[1]->name.compare("y") == 0 );
  ezassert( vars[2]->name.compare("temp") == 0 );
  ezassert( vars[3]->name.compare("time") == 0 );

  //-------------------------------------
  vars.clear();
  patterns.clear();
  patterns.push_back("T");
  regex = 1;
  regex_nocase = 1;
  coordvars = 0;
  exclude = 1;
  nc.getVars(vars, patterns, regex, regex_nocase, coordvars, exclude);
  ezassert( vars.size() == 6 );
  ezassert( vars[0]->name.compare("x") == 0 );
  ezassert( vars[1]->name.compare("y") == 0 );
  ezassert( vars[2]->name.compare("z") == 0 );
  ezassert( vars[3]->name.compare("ps") == 0 );
  ezassert( vars[4]->name.compare("rh") == 0 );
  ezassert( vars[5]->name.compare("precip") == 0 );

  //-------------------------------------
  vars.clear();
  patterns.clear();
  patterns.push_back("T");
  regex = 1;
  regex_nocase = 1;
  coordvars = 1;
  exclude = 1;
  nc.getVars(vars, patterns, regex, regex_nocase, coordvars, exclude);
  ezassert( vars.size() == 7 );
  ezassert( vars[0]->name.compare("x") == 0 );
  ezassert( vars[1]->name.compare("y") == 0 );
  ezassert( vars[2]->name.compare("z") == 0 );
  ezassert( vars[3]->name.compare("ps") == 0 );
  ezassert( vars[4]->name.compare("rh") == 0 );
  ezassert( vars[5]->name.compare("precip") == 0 );
  ezassert( vars[6]->name.compare("time") == 0 );

  ezassert( nc.close() == 0);

  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );
  
  return true;
}

bool test_hasrec(ez::ezTestRunner& runner) {
  ez::ezNc nc;

  ezassert( nc.hasrec() == 0 );

  ezassert( createFile1(runner) );
  nc.filename = filename1(runner);

  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.hasrec() );
  
  ezassert( nc.loadAllMetadata() == 0 );

  ezassert( nc.hasrec() );

  ezassert( nc.close() == 0);

  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );
  
  return true;
}

bool test_init(ez::ezTestRunner& runner) {
  ez::ezNc nc;
  ez::ezNcBuffers buf;
  
  ezassert( nc.reserve(buf, 0) != 0 );

  ezassert( createFile1(runner) );
  nc.filename = filename1(runner);
  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );

  ezassert( nc.reserve(buf, 0) == 0 );

  std::vector<int> ids;
  ezassert( nc.reserve(buf, ids) == 0 );
  ids.push_back(3);
  ids.push_back(4);
  ezassert( nc.reserve(buf, ids) == 0 );

  std::vector<std::string> names;
  ezassert( nc.reserve(buf, names) == 0 );
  names.push_back("v7");
  names.push_back("v2");
  ezassert( nc.reserve(buf, names) == 0 );
  
  ezassert( nc.close() == 0);

  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );
  
  return true;
}

bool test_isHDF(ez::ezTestRunner& runner) {
  ez::ezNc nc1,nc2;
  ezassert( nc1.isHDF() == 0 );
  
  nc1.filename = runner.tmpdir + "./tmp1.nc";
  ezassert( nc1.createClobber() == 0 );
  ezassert( nc1.enddef() == 0 ); 
  ezassert( nc1.isHDF() == 0 );
  ezassert( nc1.close() == 0);

  nc2.filename = runner.tmpdir + "./tmp2.nc";
  ezassert( nc2.createHDFClobber() == 0);
  ezassert( nc2.enddef() == 0 );
  ezassert( nc2.isHDF() );  
  ezassert( nc2.close() == 0);

  if (runner.clean) {
    ezassert( remove(nc1.filename.c_str()) == 0 );
    ezassert( remove(nc2.filename.c_str()) == 0 );
  }
  
  return true;
}

bool test_loadDimsMetadata(ez::ezTestRunner& runner) {
  ez::ezNc nc;
  
  ezassert( createFile1(runner) );
  nc.filename = filename1(runner);

  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadDimsMetadata() == 0 );

  ezassert( nc.close() == 0);

  if (runner.clean) {
    ezassert( remove(nc.filename.c_str()) == 0 );
  }
  
  return true;  
}
  
bool test_loadAllMetadata(ez::ezTestRunner& runner) {
  ez::ezNc nc;
  
  ezassert( createFile1(runner) );
  nc.filename = filename1(runner);

  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );

  ezassert( nc.close() == 0);

  if (runner.clean) {
    ezassert( remove(nc.filename.c_str()) == 0 );
  }
  
  return true;  
}

bool test_openWritable(ez::ezTestRunner& runner) {
  ez::ezNc nc1;
  
  ezassert( createFile1(runner) );
  nc1.filename = filename1(runner);
  ezassert( nc1.openWritable() == 0 );
  ezassert( nc1.loadAllMetadata() == 0 );
  
  std::vector<std::string> vals;
  vals.push_back("1.234");
  vals.push_back("5.678");
  ezassert( nc1.putAtt(0, "New Att", NC_FLOAT, vals, 0, 1) == 0 );
  ezassert( nc1.enddef() == 0 ); 
  ezassert( nc1.close() == 0);

  nc1.clear();
  nc1.filename = filename1(runner);
  ezassert( nc1.openReadOnly() == 0 );
  ezassert( nc1.loadAllMetadata() == 0 );
  
  ez::ezNcBuffers buf;
  ezassert( nc1.getAtt(0, "New Att", buf) == 0);
  ezassert( buf.nf == 2 );
  ezassert( buf.f[0] == 1.234f );
  ezassert( buf.f[1] == 5.678f );

  ezassert( nc1.close() == 0);

  if (runner.clean) {
    ezassert( remove(nc1.filename.c_str()) == 0 );
  }
  
  return true;  
}

bool test_putAtt(ez::ezTestRunner& runner) {
  ez::ezNc nc;
  nc.filename = filename1(runner);
  ezassert( nc.createClobber() == 0 );
  
  std::vector<std::string> vals;
  ezassert( nc.putAtt(0, "should be skipped", NC_FLOAT, vals, 0, -1) == 0 );

  vals.push_back("1.234");
  ezassert( nc.putAtt(0, "string", NC_CHAR, vals, 0, 0) == 0 );

  vals[0] = "1234";
  ezassert( nc.putAtt(0, "short", NC_SHORT, vals, 0, 0) == 0 );

  ezassert( nc.putAtt(0, "int", NC_INT, vals, 0, 0) == 0 );
  
  vals[0] = "1.234";
  vals.push_back("5.678");
  ezassert( nc.putAtt(0, "two floats", NC_FLOAT, vals, 0, 1) == 0 );
  ezassert( nc.putAtt(0, "two doubles", NC_DOUBLE, vals, 0, 1) == 0 );

  ezassert( nc.enddef() == 0 ); 
  ezassert( nc.close() == 0);

  nc.clear();
  nc.filename = filename1(runner);
  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );

  ez::ezNcBuffers atts;
  ezassert( nc.getAtt(0, "string", atts) == 0 );
  ezassert( strncmp("1.234", atts.c, 5) == 0 );
  ezassert( nc.getAtt(0, "short", atts) == 0 );
  ezassert( *atts.s == 1234 );
  ezassert( nc.getAtt(0, "int", atts) == 0 );
  ezassert( *atts.i == 1234 );
  ezassert( nc.getAtt(0, "two floats", atts) == 0 );
  ezassert( atts.nf == 2 );
  ezassert( atts.f[0] == 1.234f );
  ezassert( atts.f[1] == 5.678f );
  ezassert( nc.getAtt(0, "two doubles", atts) == 0 );
  ezassert( atts.nd == 2 );
  ezassert( atts.d[0] == 1.234 );
  ezassert( atts.d[1] == 5.678 );
  
  ezassert( nc.close() == 0);

  if (runner.clean) {
    ezassert( remove(nc.filename.c_str()) == 0 );
  }
  
  return true;    
}

bool test_read(ez::ezTestRunner& runner) {
  ez::ezNc nc;
  ez::ezNcBuffers buf;
  size_t start[4], count[4];
  int n1=10, n2=20, n3=30;
  
  ezassert( createFile1(runner) );
  nc.filename = filename1(runner);
  ezassert( nc.read("v7", &buf, start, count) != 0 );
  ezassert( nc.read(nc.getVarId("v7"), &buf, start, count) != 0 );
  ezassert( nc.read((char*)0, &buf, start, count) != 0 );
  ezassert( nc.read((int)0, &buf, start, count) != 0 );

  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );

  int varid = nc.getVarId("v3");
  ezassert( nc.reserve(buf, varid) == 0 );
  buf.allocate();
  buf.fill(78);
  start[0] = 0;
  start[1] = start[2] = start[3] = 0;
  count[0] = n1;
  count[1] = n2;
  count[2] = n3;
  ezassert( nc.read(varid, &buf, start, count) == 0 );
  ezassert( buf.s[0] != 78 );
  ezassert( buf.s[n1*n2*n3-1] != 78 );
  
  ezassert( nc.reserve(buf, "v7") == 0 );
  buf.allocate();
  buf.fill(78);
  start[0] = 1;
  start[1] = start[2] = start[3] = 0;
  count[0] = 1;
  count[1] = n1;
  count[2] = n2;
  count[3] = n3;
  ezassert( nc.read("v7", &buf, start, count) == 0 );
  ezassert( buf.d[0] != 78 );
  ezassert( buf.d[n1*n2*n3-1] != 78 );
  
  ezassert( nc.close() == 0);

  if (runner.clean) {
    ezassert( remove(nc.filename.c_str()) == 0 );
  }
  
  return true;
}

bool test_clear(ez::ezTestRunner& runner) {
  ez::ezNc nc;
  
  ezassert( createFile1(runner) );
  nc.filename = filename1(runner);
  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );
  nc.close();
  nc.clear();
  ezassert( nc.filename.empty() );
  ezassert( nc.ncid < 0 );
  ezassert( nc.recid < 0 );
  ezassert( nc.nrec == 0 );
  ezassert( nc.varsNameMap.empty() );
  ezassert( nc.dimsNameMap.empty() );
  ezassert( nc.varsIdVector.empty() );
  ezassert( nc.dimsIdVector.empty() );
  ezassert( nc.fill == false );
  ezassert( nc.initialsz == 0 );
  
  nc.filename = filename1(runner);

  if (runner.clean) {
    ezassert( remove(nc.filename.c_str()) == 0 );
  }
  
  return true;
}

bool test_sizeOf(ez::ezTestRunner& runner) {
  ezassert( ez::sizeOf(NC_CHAR) == 1);
  ezassert( ez::sizeOf(NC_BYTE) == 1);
  ezassert( ez::sizeOf(NC_UBYTE) == 1);
  ezassert( ez::sizeOf(NC_SHORT) == 2);
  ezassert( ez::sizeOf(NC_USHORT) == 2);
  ezassert( ez::sizeOf(NC_INT) == 4);
  ezassert( ez::sizeOf(NC_UINT) == 4);
  ezassert( ez::sizeOf(NC_INT64) == 8);
  ezassert( ez::sizeOf(NC_UINT64) == 8);
  ezassert( ez::sizeOf(NC_FLOAT) == 4);
  ezassert( ez::sizeOf(NC_DOUBLE) == 8);
  return true;
}

bool test_typeToString(ez::ezTestRunner& runner) {
  std::string s;
  
  ez::ezNc::typeToString(NC_CHAR, s);
  ezassert( s.compare("char") == 0 );
  ez::ezNc::typeToString(NC_BYTE, s);
  ezassert( s.compare("unsigned char") == 0 );
  ez::ezNc::typeToString(NC_UBYTE, s);
  ezassert( s.compare("unsigned char") == 0 );
  ez::ezNc::typeToString(NC_SHORT, s);
  ezassert( s.compare("short") == 0 );
  ez::ezNc::typeToString(NC_USHORT, s);
  ezassert( s.compare("unsigned short") == 0 );
  ez::ezNc::typeToString(NC_INT, s);
  ezassert( s.compare("int") == 0 );
  ez::ezNc::typeToString(NC_UINT, s);
  ezassert( s.compare("unsigned int") == 0 );
  ez::ezNc::typeToString(NC_INT64, s);
  ezassert( s.compare("long long") == 0 );
  ez::ezNc::typeToString(NC_UINT64, s);
  ezassert( s.compare("unsigned long long") == 0 );
  ez::ezNc::typeToString(NC_FLOAT, s);
  ezassert( s.compare("float") == 0 );
  ez::ezNc::typeToString(NC_DOUBLE, s);
  ezassert( s.compare("double") == 0 );

  return true;
}

bool test_write(ez::ezTestRunner& runner) {
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
  size_t start[4], count[4];

  ezassert( nc.reserve(buf, v1) == 0 );
  buf.allocate();
  buf.fill(1);
  start[0] = 0; start[1] = start[2] = start[3] = 0;
  count[0] = n1;
  ezassert( nc.write(v1, &buf, start, count) == 0 );

  ezassert( nc.reserve(buf, v2) == 0 );
  buf.allocate();
  buf.fill(2);
  count[0] = n1; count[1] = n2;
  ezassert( nc.write(v2, &buf, start, count) == 0 );

  ezassert( nc.reserve(buf, v3) == 0 );
  buf.allocate();
  buf.fill(3);
  count[0] = n1; count[1] = n2; count[2] = n3;
  ezassert( nc.write(v3, &buf, start, count) == 0 );

  ezassert( nc.reserve(buf, v4) == 0 );
  buf.allocate();
  buf.fill(4);
  count[0] = 1;
  ezassert( nc.write(v4, &buf, start, count) == 0 );

  ezassert( nc.reserve(buf, v5) == 0 );
  buf.allocate();
  buf.fill(5);
  count[0] = 1; count[1] = n1;
  ezassert( nc.write(v5, &buf, start, count) == 0 );

  ezassert( nc.reserve(buf, v6) == 0 );
  buf.allocate();
  buf.fill(6);
  count[0] = 1; count[1] = n1; count[2] = n2; count[3] = n3;
  ezassert( nc.write(v6, &buf, start, count) == 0 );

  ezassert( nc.reserve(buf, v7) == 0 );
  buf.allocate();
  buf.fill(7);
  start[0] = 1; 
  count[0] = 1; count[1] = n1; count[2] = n2; count[3] = n3;
  ezassert( nc.write(v7, &buf, start, count) == 0 );
  
  ezassert( nc.close() == 0 );

  nc.clear();
  nc.filename = filename1(runner);
  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );
  
  ezassert( nc.reserve(buf, "v1") == 0 );
  buf.allocate();
  buf.fill(78);
  start[0] = 0;
  count[0] = n1; 
  ezassert( nc.read("v1", &buf, start, count) == 0 );
  ezassert( buf.uc[0] == 1 );
  ezassert( buf.uc[n1-1] == 1 );

  ezassert( nc.reserve(buf, "v2") == 0 );
  buf.allocate();
  buf.fill(78);
  count[0] = n1; count[1] = n2; 
  ezassert( nc.read("v2", &buf, start, count) == 0 );
  ezassert( buf.c[0] == 2 );
  ezassert( buf.c[n1*n2-1] == 2 );

  ezassert( nc.reserve(buf, "v3") == 0 );
  buf.allocate();
  buf.fill(78);
  count[0] = n1; count[1] = n2; count[2] = n3; 
  ezassert( nc.read("v3", &buf, start, count) == 0 );
  ezassert( buf.s[0] == 3 );
  ezassert( buf.s[n1*n2*n3-1] == 3 );

  ezassert( nc.reserve(buf, "v4") == 0 );
  buf.allocate();
  buf.fill(78);
  start[0] = 0;
  count[0] = 1;
  ezassert( nc.read("v4", &buf, start, count) == 0 );
  ezassert( buf.i[0] == 4 );

  ezassert( nc.reserve(buf, "v5") == 0 );
  buf.allocate();
  buf.fill(78);
  start[0] = 0;
  count[0] = 1; count[1] = n1;
  ezassert( nc.read("v5", &buf, start, count) == 0 );
  ezassert( buf.f[0] == 5 );
  ezassert( buf.f[n1-1] == 5 );

  ezassert( nc.reserve(buf, "v6") == 0 );
  buf.allocate();
  buf.fill(78);
  start[0] = 0;
  count[0] = 1; count[1] = n1; count[2] = n2;
  ezassert( nc.read("v6", &buf, start, count) == 0 );
  ezassert( buf.d[0] == 6 );
  ezassert( buf.d[n1*n2-1] == 6 );

  ezassert( nc.reserve(buf, "v7") == 0 );
  buf.allocate();
  buf.fill(78);
  start[0] = 1;
  count[0] = 1; count[1] = n1; count[2] = n2; count[3] = n3;
  ezassert( nc.read("v7", &buf, start, count) == 0 );
  ezassert( buf.d[0] == 7 );
  ezassert( buf.d[n1*n2*n3-1] == 7 );

  nc.close();
  
  if (runner.clean) {
    ezassert( remove(nc.filename.c_str()) == 0 );
  }

  return true;
}
  
int main(int argc, const char* argv[]) {
  ez::ezTestRunner runner;
  runner.getopt(argc, argv);
  #define TEST(N) runner.add(&test_##N, "test_"#N);
  TEST(close);
  TEST(computeSize);
  TEST(copyAtt);
  TEST(copyDimDef);
  TEST(copyDimDefs);
  TEST(copySettings);
  TEST(copyVarDef);
  TEST(copyVarDefs);
  TEST(create64);
  TEST(create64Clobber);
  TEST(createHDF);
  TEST(createHDFClobber);
  TEST(createNcDim);
  TEST(createDim);
  TEST(createNcVar);
  TEST(createVar);
  TEST(formatToString);
  TEST(getAllVarBytes);
  TEST(getAllVarIds);
  TEST(getAllVarNames);
  TEST(getDim);
  TEST(getDimIds);
  TEST(getDimSizesId);
  TEST(getDimSizesString);
  TEST(getDimNamesForVars);
  TEST(getMaxNumDims);
  TEST(getNumVars);
  TEST(getNumStaticVars);
  TEST(getNumRecVars);
  TEST(getRecVarBytes);
  TEST(getVar);
  TEST(getRecVarNames);
  TEST(getRecVarIds);  
  TEST(getStaticVarBytes);
  TEST(getStaticVarNames);
  TEST(getStaticVarIds);
  TEST(getVarSizeTZYX);
  TEST(getVars);
  TEST(hasrec);
  TEST(init);
  TEST(isHDF);
  TEST(loadDimsMetadata);
  TEST(loadAllMetadata);
  TEST(openWritable);
  TEST(putAtt);
  TEST(read);
  TEST(clear);  
  TEST(sizeOf);  
  TEST(typeToString);  
  TEST(write);
  
  return runner.run();
}