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
  int dr = nc.createDim("dr", 1, true);
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
  int dr = nc.createDim("dr", 1, true);
  ezassert(dr != -1);
  
  int dimids1d1[1] = {d1};
  int dimids1d2[1] = {d2};
  int dimids1dr[1] = {dr};
  int dimids2d[2] = {d1,d2};
  int vt;
  
  ezassert( nc.createVar("v1", NC_FLOAT, 1, dimids1d1) != -1 );
  ezassert( nc.createVar("v2", NC_FLOAT, 1, dimids1d2) != -1 );
  ezassert( nc.createVar("v3", NC_FLOAT, 1, dimids1dr) != -1 );
  ezassert( (vt = nc.createVar(tname, NC_FLOAT, 1, dimids1dr)) != -1 );
  ezassert( nc.createVar("v4", NC_FLOAT, 1, dimids1dr) != -1 );
  ezassert( nc.createVar("v5", NC_FLOAT, 1, dimids2d) != -1 );

  std::vector<std::string> strings;
  strings.resize(1);
  strings[0] = attval;
  ezassert( nc.putAtt(tname, att, NC_CHAR, strings, 0, 0) == 0);

  ezassert( nc.close() == 0 );
  
  return true;
}

bool createFile1dRecTimeAttUnits(ez::ezTestRunner& runner, const char* filenum, const char* tname, const char* unitsString, const char* calendar, float toffset=0.0) {
  ez::ezNc nc;
  
  nc.clear();
  nc.filename = filename(runner, filenum);
  ezassert( nc.createClobber() == 0 );
  
  int n1=10; 
  int nt=10;
  int d1 = nc.createDim("d1", n1);
  ezassert(d1 != -1);
  int d2 = nc.createDim("d2", n1);
  ezassert(d2 != -1);
  int dr = nc.createDim("dr", nt, true);
  ezassert(dr != -1);
  
  int dimids1d1[1] = {d1};
  int dimids1d2[1] = {d2};
  int dimids1dr[1] = {dr};
  int dimids2d[2] = {d1,d2};
  int vt;
  
  ezassert( nc.createVar("v1", NC_FLOAT, 1, dimids1d1) != -1 );
  ezassert( nc.createVar("v2", NC_FLOAT, 1, dimids1d2) != -1 );
  ezassert( nc.createVar("v3", NC_FLOAT, 1, dimids1dr) != -1 );
  ezassert( (vt = nc.createVar(tname, NC_FLOAT, 1, dimids1dr)) != -1 );
  ezassert( nc.createVar("v4", NC_FLOAT, 1, dimids1dr) != -1 );
  ezassert( nc.createVar("v5", NC_FLOAT, 1, dimids2d) != -1 );

  std::vector<std::string> strings;
  strings.resize(1);
  strings[0] = "time";
  ezassert( nc.putAtt(tname, "standard_name", NC_CHAR, strings, 0, 0) == 0);

  strings[0] = unitsString;
  if (strings[0].size())
    ezassert( nc.putAtt(tname, "units", NC_CHAR, strings, 0, 0) == 0);

  strings[0] = calendar;
  if (strings[0].size())
    ezassert( nc.putAtt(tname, "calendar", NC_CHAR, strings, 0, 0) == 0);
  
  ezassert( nc.enddef() == 0 );
  ez::ezNcBuffers buf;
  size_t start[1], count[1];
  
  // Write rec dim data.
  ezassert( nc.reserve(buf, vt) == 0 );
  buf.allocate();

  int i;
  count[0] = 1;
  for(i=0; i < nt; ++i) {
    buf.f[0] = i + toffset;
    start[0] = i;
    ezassert( nc.write(vt, &buf, start, count) == 0 );
  }
  
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
  int dr = nc.createDim("dr", 1, true);
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
  ezassert( d.baselineValue == 0.0 );

  u = "years since 0000-01-01 00:00";
  d.setUnits(u);
  ezassert( d.baselineValue == 1.0 );
  ezassert( d.unit == ez::ezDateTime::YEARS );
  ezassert( d.tzHours == 0 );
  ezassert( d.tzMinutes == 0 );

  u = "days since 1980-12-31 18:12:34 +3:30";
  d.setUnits(u);
  ezassert( d.baselineValue == 1980365.6128935185 );
  ezassert( d.unit == ez::ezDateTime::DAYS );
  ezassert( d.tzHours == 3 );
  ezassert( d.tzMinutes == 30 );

//  printf("baselineValue = %g\n", d.baselineValue);
  
  return true;
}

bool test_fromUnits(ez::ezTestRunner& runner) {
  ez::ezNcDateTime d;
  std::string u;
  double t[] = {0,1,2};
  int nt = 3;
  std::vector<double> v;
  
  // Convert days "t" array to double values relative to "u" baseline.
  u = "years since 1980-12-31 18:12:34 +3:30";
  d.setUnits(u);
  d.fromUnits<double>(t, nt, v);
  ezassert( d.unit == ez::ezDateTime::YEARS );
  ezassert( ez::round_frac(v[0],9) == 1980365.612893518 );
  ezassert( ez::round_frac(v[1],9) == 1981365.612893518 );
  ezassert( ez::round_frac(v[2],9) == 1982365.612893518 );

  u = "months since 1980-1-1 18:12:34 +3:30";
  d.setUnits(u);
  d.fromUnits<double>(t, nt, v);
  ezassert( d.unit == ez::ezDateTime::MONTHS );
  ezassert( ez::round_frac(v[0],9) == 1980001.612893518 );
  ezassert( ez::round_frac(v[1],9) == 1980032.612893518 );
  ezassert( ez::round_frac(v[2],9) == 1980060.612893518 );
  
  u = "days since 1980-12-31 18:12:34 +3:30";
  d.setUnits(u);
  d.fromUnits<double>(t, nt, v);
  ezassert( d.unit == ez::ezDateTime::DAYS );
  ezassert( ez::round_frac(v[0],9) == 1980365.612893518 );
  ezassert( ez::round_frac(v[1],9) == 1981001.612893518 );
  ezassert( ez::round_frac(v[2],9) == 1981002.612893518 );

  d.clear();
  u = "hours since 1980-12-31 18:12:34 +3:30";
  d.setUnits(u);
  d.fromUnits<double>(t, nt, v);
  ezassert( d.unit == ez::ezDateTime::HOURS );
  ezassert( ez::round_frac(v[0],10) == 1980365.6128935183 );
  ezassert( ez::round_frac(v[1],10) == 1980365.6545601853 );
  ezassert( ez::round_frac(v[2],10) == 1980365.6962268515 );

  d.clear();
  u = "minutes since 1980-12-31 18:12:34 +3:30";
  d.setUnits(u);
  d.fromUnits<double>(t, nt, v);
  ezassert( d.unit == ez::ezDateTime::MINUTES );
  ezassert( ez::round_frac(v[0],10) == 1980365.6128935183 );
  ezassert( ez::round_frac(v[1],10) == 1980365.6135879627 );
  ezassert( ez::round_frac(v[2],10) == 1980365.6142824071 );

  d.clear();
  u = "seconds since 1980-12-31 18:12:34 +3:30";
  d.setUnits(u);
  d.fromUnits<double>(t, nt, v);
  ezassert( d.unit == ez::ezDateTime::SECONDS );
  ezassert( ez::round_frac(v[0],10) == 1980365.6128935183 );
  ezassert( ez::round_frac(v[1],10) == 1980365.6129050923 );
  ezassert( ez::round_frac(v[2],10) == 1980365.6129166663 );

//  printf("v = %.10f %.10f %.10f\n", v[0], v[1], v[2]);
//  printf("r = %.10f %.10f %.10f\n", ez::round_frac(v[0],10), ez::round_frac(v[1],10), ez::round_frac(v[2],10));
  
  return true;
}

bool test_build(ez::ezTestRunner& runner) {
  ez::ezNcDateTime d;
  std::string u;
  ez::ezNc nc;
  
  ezassert( createFile1dRecTimeAttUnits(runner, "111", "vtime", "days since 1980-12-31 18:12:34 +3:30", "") );
  nc.filename = filename(runner, "111");

  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );
  
  ezassert( d.build(&nc) == 0 );
  ezassert( d.baselineValue == 1980365.6128935185 );
  ezassert( d.calendar == ez::ezDateTime::NOLEAP );
  ezassert( d.tzHours == 3 );
  ezassert( d.tzMinutes == 30 );
  ezassert( d.unit == ez::ezDateTime::DAYS );
  
  ezassert( nc.close() == 0 );
  
  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );

  return true;
}

bool test_contains(ez::ezTestRunner& runner) {
  ez::ezNcDateTime d;
  std::string u;

  ez::ezNc nc;
  
  ezassert( createFile1dRecTimeAttUnits(runner, "111", "vtime", "days since 1980-12-31 18:12:34 +3:30", "") );
  nc.filename = filename(runner, "111");

  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );
  
  ezassert( d.build(&nc) == 0 );
  
  ezassert( d.contains(0) == 0);
  ezassert( d.contains(1980000) == 0);
  ezassert( d.contains(1980365) == 0);
  ezassert( d.contains(1981001) == 1);
  ezassert( d.contains(1981002) == 1);
  ezassert( d.contains(1981011) == 0);
  ezassert( d.contains(2000001) == 0);
  
  ezassert( nc.close() == 0 );
  
  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );

  return true;
}

bool test_getRange(ez::ezTestRunner& runner) {
  ez::ezNcDateTime d;
  int ibegin, iend;
  ez::ezNc nc;
  
  ezassert( createFile1dRecTimeAttUnits(runner, "111", "vtime", "days since 1980-12-31 18:12:34 +3:30", "") );
  nc.filename = filename(runner, "111");

  ezassert( nc.openReadOnly() == 0 );
  ezassert( nc.loadAllMetadata() == 0 );
  
  ezassert( d.build(&nc) == 0 );

  d.getRange(0, 1, ibegin, iend);
  ezassert( ibegin == -1 );
  ezassert( iend == -1 );

  d.getRange(1981001, 1981009, ibegin, iend);
  ezassert( ibegin == 1 );
  //printf("%s:%d iend=%d\nvalues[8]=%.17f\nvalues[9]=%.17f\n", __FILE__, __LINE__, iend, d.values[8], d.values[9]);
  ezassert( iend == 8 );
  
  d.getRange(1981002, 1981005, ibegin, iend);
  ezassert( ibegin == 2 );
  ezassert( iend == 4 );

  d.getRange(0, 1981005, ibegin, iend);
  ezassert( ibegin == 0 );
  ezassert( iend == 4 );

  d.getRange(1981001, 2000123, ibegin, iend);
  ezassert( ibegin == 1 );
  ezassert( iend == 9 );

  ezassert( nc.close() == 0 );
  
  if (runner.clean)
    ezassert( remove(nc.filename.c_str()) == 0 );

  return true;
}

bool test_equals(ez::ezTestRunner& runner) {
  ez::ezNcDateTime d1,d2;
  ez::ezNc nc1, nc2;
  
  ezassert( createFile1dRecTimeAttUnits(runner, "1", "vtime", "days since 1980-12-31 18:12:34 +3:30", "") );
  nc1.filename = filename(runner, "1");
  nc2.filename = filename(runner, "1");

  ezassert( nc1.openReadOnly() == 0 );
  ezassert( nc1.loadAllMetadata() == 0 );
  ezassert( d1.build(&nc1) == 0 );

  ezassert( nc2.openReadOnly() == 0 );
  ezassert( nc2.loadAllMetadata() == 0 );
  ezassert( d2.build(&nc2) == 0 );
  
  ezassert( d1.equals(d2) == 1 );
  
  ezassert( nc1.close() == 0 );
  ezassert( nc2.close() == 0 );

  //--------------------------------------------
  // Different calendar.
  d1.clear(); 
  d2.clear();
  ezassert( createFile1dRecTimeAttUnits(runner, "2", "vtime", "days since 1980-12-31 18:12:34 +3:30", "ALL_LEAP") );
  
  nc1.filename = filename(runner, "1");
  nc2.filename = filename(runner, "2");

  ezassert( nc1.openReadOnly() == 0 );
  ezassert( nc1.loadAllMetadata() == 0 );
  ezassert( d1.build(&nc1) == 0 );

  ezassert( nc2.openReadOnly() == 0 );
  ezassert( nc2.loadAllMetadata() == 0 );
  ezassert( d2.build(&nc2) == 0 );
  
  ezassert( d1.equals(d2) == 0 );
  ezassert( d1.valuesEqual(d2) == 0 );
  
  ezassert( nc1.close() == 0 );
  ezassert( nc2.close() == 0 );

  //--------------------------------------------
  // Different unit.
  d1.clear(); 
  d2.clear();
  ezassert( createFile1dRecTimeAttUnits(runner, "3", "vtime", "hours since 1980-12-31 18:12:34 +3:30", "ALL_LEAP") );
  
  nc1.filename = filename(runner, "1");
  nc2.filename = filename(runner, "3");

  ezassert( nc1.openReadOnly() == 0 );
  ezassert( nc1.loadAllMetadata() == 0 );
  ezassert( d1.build(&nc1) == 0 );

  ezassert( nc2.openReadOnly() == 0 );
  ezassert( nc2.loadAllMetadata() == 0 );
  ezassert( d2.build(&nc2) == 0 );
  
  ezassert( d1.equals(d2) == 0 );
  
  ezassert( nc1.close() == 0 );
  ezassert( nc2.close() == 0 );

  //--------------------------------------------
  // Different timezone hours, but same UTC time.
  d1.clear(); 
  d2.clear();
  ezassert( createFile1dRecTimeAttUnits(runner, "4", "vtime", "days since 1980-12-31 17:12:34 +2:30", "") );
  
  nc1.filename = filename(runner, "1");
  nc2.filename = filename(runner, "4");

  ezassert( nc1.openReadOnly() == 0 );
  ezassert( nc1.loadAllMetadata() == 0 );
  ezassert( d1.build(&nc1) == 0 );

  ezassert( nc2.openReadOnly() == 0 );
  ezassert( nc2.loadAllMetadata() == 0 );
  ezassert( d2.build(&nc2) == 0 );
  
  ezassert( d1.equals(d2) == 0 );
  ezassert( d1.valuesEqual(d2) == 1 );
  
  ezassert( nc1.close() == 0 );
  ezassert( nc2.close() == 0 );

  //--------------------------------------------
  // Different timezone minutes, but same UTC time.
  d1.clear(); 
  d2.clear();
  ezassert( createFile1dRecTimeAttUnits(runner, "5", "vtime", "days since 1980-12-31 18:02:34 +3:20", "") );
  
  nc1.filename = filename(runner, "1");
  nc2.filename = filename(runner, "5");

  ezassert( nc1.openReadOnly() == 0 );
  ezassert( nc1.loadAllMetadata() == 0 );
  ezassert( d1.build(&nc1) == 0 );

  ezassert( nc2.openReadOnly() == 0 );
  ezassert( nc2.loadAllMetadata() == 0 );
  ezassert( d2.build(&nc2) == 0 );
  
  ezassert( d1.equals(d2) == 0 );
  ezassert( d1.valuesEqual(d2) == 1 );
  
  ezassert( nc1.close() == 0 );
  ezassert( nc2.close() == 0 );

  //--------------------------------------------
  // Different values.
  d1.clear(); 
  d2.clear();
  ezassert( createFile1dRecTimeAttUnits(runner, "6", "vtime", "days since 1980-12-31 18:12:34 +3:30", "", 1.0) );
  
  nc1.filename = filename(runner, "1");
  nc2.filename = filename(runner, "6");

  ezassert( nc1.openReadOnly() == 0 );
  ezassert( nc1.loadAllMetadata() == 0 );
  ezassert( d1.build(&nc1) == 0 );

  ezassert( nc2.openReadOnly() == 0 );
  ezassert( nc2.loadAllMetadata() == 0 );
  ezassert( d2.build(&nc2) == 0 );
  
  ezassert( d1.equals(d2) == 0 );
  ezassert( d1.valuesEqual(d2) == 0 );
  
  ezassert( nc1.close() == 0 );
  ezassert( nc2.close() == 0 );

  if (runner.clean) {
    ezassert( remove(filename(runner, "1").c_str()) == 0 );
    ezassert( remove(filename(runner, "2").c_str()) == 0 );
    ezassert( remove(filename(runner, "3").c_str()) == 0 );
    ezassert( remove(filename(runner, "4").c_str()) == 0 );
    ezassert( remove(filename(runner, "5").c_str()) == 0 );
    ezassert( remove(filename(runner, "6").c_str()) == 0 );
  }
  
  return true;
}

bool test_compare(ez::ezTestRunner& runner) {
  ez::ezNcDateTime d1,d2;
  ez::ezNc nc1, nc2;
 
  ezassert( createFile1dRecTimeAttUnits(runner, "1", "vtime", "days since 1980-12-31 18:12:34 +3:30", "NOLEAP", 0.0) );
  ezassert( createFile1dRecTimeAttUnits(runner, "2", "vtime", "days since 1980-12-31 18:12:34 +3:30", "NOLEAP", 5.0) );

  nc1.filename = filename(runner, "1");
  nc2.filename = filename(runner, "2");

  ezassert( nc1.openReadOnly() == 0 );
  ezassert( nc1.loadAllMetadata() == 0 );
  ezassert( d1.build(&nc1) == 0 );

  ezassert( nc2.openReadOnly() == 0 );
  ezassert( nc2.loadAllMetadata() == 0 );
  ezassert( d2.build(&nc2) == 0 );
  
  ezassert( ez::ezNcDateTime::compare(d2,d1) == 0 );
  ezassert( ez::ezNcDateTime::compare(d1,d2) == 1 );
 
  ezassert( nc1.close() == 0 );
  ezassert( nc2.close() == 0 );

  if (runner.clean) {
    ezassert( remove(filename(runner, "1").c_str()) == 0 );
    ezassert( remove(filename(runner, "2").c_str()) == 0 );
  }
  
  return true;
}

int main(int argc, const char* argv[]) {
  ez::ezTestRunner runner;
  runner.getopt(argc, argv);
  #define TEST(N) runner.add(&test_##N, "test_"#N);
  TEST(autodetect);
  TEST(setUnits);
  TEST(fromUnits);
  TEST(build);
  TEST(contains);
  TEST(getRange);
  TEST(equals);
  TEST(compare);
  
  return runner.run();
}