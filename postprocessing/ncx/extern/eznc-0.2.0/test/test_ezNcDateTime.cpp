/*
20111204 rsz Created.
*/

#include "ezTest.hpp"
#include "ezNcDateTime.hpp"

std::string filename(ez::ezTestRunner& runner, const char * i) {
  std::string f = runner.tmpdir + "./tmp";
  f.append(i);
  f.append(".nc");
  return f;
}

bool createFile1dStaticTime(ez::ezTestRunner& runner, const char* i, const char* tname) {
  ez::ezNc nc;
  
  nc.clear();
  nc.filename = filename(runner, i);
  ezassert( nc.createClobber() == 0 );
  
  int n1=10;
  int d1 = nc.createDim("d1", n1);
  ezassert(d1 != -1);
  int d2 = nc.createDim("d2", n1);
  ezassert(d2 != -1);
  
  int dimids1d1[1] = {d1};
  int dimids1d2[1] = {d2};
  int dimids2d[2] = {d1,d2};
  
  ezassert( nc.createVar("v1", NC_FLOAT, 1, dimids1d1) != -1 );
  ezassert( nc.createVar("v2", NC_FLOAT, 1, dimids1d2) != -1 );
  ezassert( nc.createVar(tname, NC_FLOAT, 1, dimids1d2) != -1 );
  ezassert( nc.createVar("v3", NC_FLOAT, 1, dimids2d) != -1 );

  ezassert( nc.close() == 0 );
  
  return true;
}

bool createFile1dRecTime(ez::ezTestRunner& runner, const char* i, const char* tname) {
  ez::ezNc nc;
  
  nc.clear();
  nc.filename = filename(runner, i);
  ezassert( nc.createClobber() == 0 );
  
  int n1=10;
  int d1 = nc.createDim("d1", n1);
  ezassert(d1 != -1);
  int d2 = nc.createDim("d2", n1);
  ezassert(d2 != -1);
  int dr = nc.createDim("dr", NC_UNLIMITED);
  ezassert(dr != -1);
  
  int dimids1d1[1] = {d1};
  int dimids1d2[1] = {d2};
  int dimids1dr[1] = {dr};
  int dimids2d[2] = {d1,d2};
  
  ezassert( nc.createVar("v1", NC_FLOAT, 1, dimids1d1) != -1 );
  ezassert( nc.createVar("v2", NC_FLOAT, 1, dimids1d2) != -1 );
  ezassert( nc.createVar("v3", NC_FLOAT, 1, dimids1dr) != -1 );
  ezassert( nc.createVar(tname, NC_FLOAT, 1, dimids1dr) != -1 );
  ezassert( nc.createVar("v4", NC_FLOAT, 1, dimids1dr) != -1 );
  ezassert( nc.createVar("v5", NC_FLOAT, 1, dimids2d) != -1 );

  ezassert( nc.close() == 0 );
  
  return true;
}

bool createFile1dRecTimeAtt(ez::ezTestRunner& runner, const char* i, const char* tname, const char* att, const char* attval) {
  ez::ezNc nc;
  
  nc.clear();
  nc.filename = filename(runner, i);
  ezassert( nc.createClobber() == 0 );
  
  int n1=10;
  int d1 = nc.createDim("d1", n1);
  ezassert(d1 != -1);
  int d2 = nc.createDim("d2", n1);
  ezassert(d2 != -1);
  int dr = nc.createDim("dr", NC_UNLIMITED);
  ezassert(dr != -1);
  
  int dimids1d1[1] = {d1};
  int dimids1d2[1] = {d2};
  int dimids1dr[1] = {dr};
  int dimids2d[2] = {d1,d2};
  
  ezassert( nc.createVar("v1", NC_FLOAT, 1, dimids1d1) != -1 );
  ezassert( nc.createVar("v2", NC_FLOAT, 1, dimids1d2) != -1 );
  ezassert( nc.createVar("v3", NC_FLOAT, 1, dimids1dr) != -1 );
  ezassert( nc.createVar(tname, NC_FLOAT, 1, dimids1dr) != -1 );
  ezassert( nc.createVar("v4", NC_FLOAT, 1, dimids1dr) != -1 );
  ezassert( nc.createVar("v5", NC_FLOAT, 1, dimids2d) != -1 );

  std::vector<std::string> strings;
  strings.resize(1);
  strings[0] = attval;
  ezassert( nc.putAtt(tname, att, NC_CHAR, strings, 0, 0) == 0);

  ezassert( nc.close() == 0 );
  
  return true;
}

bool createFile1dStaticTimeAtt(ez::ezTestRunner& runner, const char* i, const char* tname, const char* att, const char* attval) {
  ez::ezNc nc;
  
  nc.clear();
  nc.filename = filename(runner, i);
  ezassert( nc.createClobber() == 0 );
  
  int n1=10;
  int d1 = nc.createDim("d1", n1);
  ezassert(d1 != -1);
  int d2 = nc.createDim("d2", n1);
  ezassert(d2 != -1);
  int dr = nc.createDim("dr", NC_UNLIMITED);
  ezassert(dr != -1);
  
  int dimids1d1[1] = {d1};
  int dimids1d2[1] = {d2};
  int dimids1dr[1] = {dr};
  int dimids2d[2] = {d1,d2};
  
  ezassert( nc.createVar("v1", NC_FLOAT, 1, dimids1d1) != -1 );
  ezassert( nc.createVar(tname, NC_FLOAT, 1,dimids1d1) != -1 );
  ezassert( nc.createVar("v2", NC_FLOAT, 1, dimids1d2) != -1 );
  ezassert( nc.createVar("v3", NC_FLOAT, 1, dimids1dr) != -1 );
  ezassert( nc.createVar("v4", NC_FLOAT, 1, dimids1dr) != -1 );
  ezassert( nc.createVar("v5", NC_FLOAT, 1, dimids2d) != -1 );

  std::vector<std::string> strings;
  strings.resize(1);
  strings[0] = attval;
  ezassert( nc.putAtt(tname, att, NC_CHAR, strings, 0, 0) == 0);

  ezassert( nc.close() == 0 );
  
  return true;
}

bool test_autodetect(ez::ezTestRunner& runner) {
  ez::ezNcDateTime dt;
  ez::ezNcVar* v;
  ez::ezNc nc;

  //-----------------------------------
  // Test simple static 1d vars with a "time" name.
  // Note, these are detected as last resort if all CF conventions fail.
  ezassert( createFile1dStaticTime(runner, "1", "time") );
  nc.filename = filename(runner, "1");

  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );
  
  v = dt.autodetect(&nc);
  ezassert( v );
  ezassert( v->name.compare("time") == 0 );
  
  ezassert( nc.close() == 0 );
  
  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );

  //-------------------------------------
  nc.clear();
  ezassert( createFile1dStaticTime(runner, "2", "rec") );
  nc.filename = filename(runner, "2");

  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );
  
  v = dt.autodetect(&nc);
  ezassert( v );
  ezassert( v->name.compare("rec") == 0 );
  
  ezassert( nc.close() == 0 );
  
  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );

  //-------------------------------------
  nc.clear();
  ezassert( createFile1dStaticTime(runner, "3", "record") );
  nc.filename = filename(runner, "3");

  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );
  
  v = dt.autodetect(&nc);
  ezassert( v );
  ezassert( v->name.compare("record") == 0 );
  
  ezassert( nc.close() == 0 );
  
  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );

  //-------------------------------------
  nc.clear();
  ezassert( createFile1dStaticTime(runner, "4", "t") );
  nc.filename = filename(runner, "4");

  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );
  
  v = dt.autodetect(&nc);
  ezassert( v );
  ezassert( v->name.compare("t") == 0 );
  
  ezassert( nc.close() == 0 );
  
  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );

  //--------------------------------------
  // Simple non-CF case, but with rec vars.
  ezassert( createFile1dRecTime(runner, "11", "time") );
  nc.filename = filename(runner, "11");

  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );
  
  v = dt.autodetect(&nc);
  ezassert( v );
  ezassert( v->name.compare("time") == 0 );
  
  ezassert( nc.close() == 0 );
  
  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );

  //--------------------------------------
  ezassert( createFile1dRecTime(runner, "12", "t") );
  nc.filename = filename(runner, "12");

  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );
  
  v = dt.autodetect(&nc);
  ezassert( v );
  ezassert( v->name.compare("t") == 0 );
  
  ezassert( nc.close() == 0 );
  
  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );

  //--------------------------------------
  ezassert( createFile1dRecTime(runner, "13", "record") );
  nc.filename = filename(runner, "13");

  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );
  
  v = dt.autodetect(&nc);
  ezassert( v );
  ezassert( v->name.compare("record") == 0 );
  
  ezassert( nc.close() == 0 );
  
  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );

  //--------------------------------------
  ezassert( createFile1dRecTime(runner, "14", "rec") );
  nc.filename = filename(runner, "14");

  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );
  
  v = dt.autodetect(&nc);
  ezassert( v );
  ezassert( v->name.compare("rec") == 0 );
  
  ezassert( nc.close() == 0 );
  
  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );

  //--------------------------------------
  // CF cases with rec vars.
  ezassert( createFile1dRecTimeAtt(runner, "111", "vtime", "standard_name", "time") );
  nc.filename = filename(runner, "111");

  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );
  
  v = dt.autodetect(&nc);
  ezassert( v );
  ezassert( v->name.compare("vtime") == 0 );
  
  ezassert( nc.close() == 0 );
  
  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );

  //--------------------------------------
  // CF cases with rec vars.
  ezassert( createFile1dRecTimeAtt(runner, "111", "vtime", "standard_name", "time") );
  nc.filename = filename(runner, "111");

  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );
  
  v = dt.autodetect(&nc);
  ezassert( v );
  ezassert( v->name.compare("vtime") == 0 );
  
  ezassert( nc.close() == 0 );
  
  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );

  //--------------------------------------
  ezassert( createFile1dRecTimeAtt(runner, "112", "vtime", "long_name", "time") );
  nc.filename = filename(runner, "112");

  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );
  
  v = dt.autodetect(&nc);
  ezassert( v );
  ezassert( v->name.compare("vtime") == 0 );
  
  ezassert( nc.close() == 0 );
  
  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );

  //--------------------------------------
  ezassert( createFile1dRecTimeAtt(runner, "113", "vtime", "axis", "T") );
  nc.filename = filename(runner, "113");

  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );
  
  v = dt.autodetect(&nc);
  ezassert( v );
  ezassert( v->name.compare("vtime") == 0 );
  
  ezassert( nc.close() == 0 );
  
  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );

  //--------------------------------------
  ezassert( createFile1dRecTimeAtt(runner, "114", "vtime", "units", "hours since 1980-0-0") );
  nc.filename = filename(runner, "114");

  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );
  
  v = dt.autodetect(&nc);
  ezassert( v );
  ezassert( v->name.compare("vtime") == 0 );
  
  ezassert( nc.close() == 0 );
  
  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );

  //--------------------------------------
  ezassert( createFile1dRecTimeAtt(runner, "115", "vtime", "units", "days since 1980-0-0") );
  nc.filename = filename(runner, "115");

  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );
  
  v = dt.autodetect(&nc);
  ezassert( v );
  ezassert( v->name.compare("vtime") == 0 );
  
  ezassert( nc.close() == 0 );
  
  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );

  //--------------------------------------
  ezassert( createFile1dRecTimeAtt(runner, "116", "vtime", "units", "years since 1980-0-0") );
  nc.filename = filename(runner, "116");

  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );
  
  v = dt.autodetect(&nc);
  ezassert( v );
  ezassert( v->name.compare("vtime") == 0 );
  
  ezassert( nc.close() == 0 );
  
  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );

  //--------------------------------------
  ezassert( createFile1dRecTimeAtt(runner, "117", "vtime", "calendar", "julian") );
  nc.filename = filename(runner, "117");

  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );
  
  v = dt.autodetect(&nc);
  ezassert( v );
  ezassert( v->name.compare("vtime") == 0 );
  
  ezassert( nc.close() == 0 );
  
  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );

  //--------------------------------------
  // CF cases with static vars.
  ezassert( createFile1dStaticTimeAtt(runner, "1111", "vtime", "standard_name", "time") );
  nc.filename = filename(runner, "1111");

  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );
  
  v = dt.autodetect(&nc);
  ezassert( v );
  ezassert( v->name.compare("vtime") == 0 );
  
  ezassert( nc.close() == 0 );
  
  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );

  //--------------------------------------
  ezassert( createFile1dStaticTimeAtt(runner, "1112", "vtime", "axis", "T") );
  nc.filename = filename(runner, "1112");

  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );
  
  v = dt.autodetect(&nc);
  ezassert( v );
  ezassert( v->name.compare("vtime") == 0 );
  
  ezassert( nc.close() == 0 );
  
  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );

  //--------------------------------------
  ezassert( createFile1dStaticTimeAtt(runner, "1113", "vtime", "units", "days since 2000-0-0") );
  nc.filename = filename(runner, "1113");

  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );
  
  v = dt.autodetect(&nc);
  ezassert( v );
  ezassert( v->name.compare("vtime") == 0 );
  
  ezassert( nc.close() == 0 );
  
  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );

  //--------------------------------------
  ezassert( createFile1dStaticTimeAtt(runner, "1114", "vtime", "calendar", "noleap") );
  nc.filename = filename(runner, "1114");

  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );
  
  v = dt.autodetect(&nc);
  ezassert( v );
  ezassert( v->name.compare("vtime") == 0 );
  
  ezassert( nc.close() == 0 );
  
  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );

  return true;
}


bool test_setUnits(ez::ezTestRunner& runner) {
  ez::ezNcDateTime d;
  std::string u;
  
  u = "";
  d.setUnits(u);
  ezassert( d.baseline.value == 0.0 );
  
  return true;
}

int main(int argc, const char* argv[]) {
  ez::ezTestRunner runner;
  runner.getopt(argc, argv);
  #define TEST(N) runner.add(&test_##N, "test_"#N);
  TEST(autodetect);
  TEST(setUnits);
  
  return runner.run();
}