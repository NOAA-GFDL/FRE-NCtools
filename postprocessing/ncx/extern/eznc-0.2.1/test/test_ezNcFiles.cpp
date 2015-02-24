/*
20120101 rsz Created.
*/

#include "ezTest.hpp"
#include "ezNcFiles.hpp"

std::string createFilename(ez::ezTestRunner& runner, const char* id) {
  std::string fname = runner.tmpdir;
  fname.append(id);
  fname.append(".nc");
  return fname;
}

bool createFile(ez::ezTestRunner& runner, std::string& fname, int toffset) {
  ez::ezNc nc;
  nc.filename = fname;  
  ezassert( nc.createClobber() == 0 );

  int n = 10;
  int dr = nc.createDim("time", n, true);
  ezassert(dr != -1);

  int dimids[1] = {dr};
  int v = nc.createVar("time", NC_DOUBLE, 1, dimids);
  ezassert(v != -1);
  
  std::vector<std::string> strings;
  strings.resize(1);
  strings[0] = "time";
  ezassert( nc.putAtt("time", "standard_name", NC_CHAR, strings, 0, 0) == 0);
  strings[0] = "years since 1980-01-01";
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

bool test_getSortedByName(ez::ezTestRunner& runner) {
  ez::ezNcFiles files;
  std::vector<std::string> names;
  std::vector<ez::ezNc*> sorted;
  
  files.getSortedByName(sorted);
  ezassert( sorted.empty() );
  
  names.resize(3);
  names[0] = createFilename(runner, "first");
  names[1] = createFilename(runner, "2nd");
  names[2] = createFilename(runner, "andlast");
  
  createFile(runner, names[0], 0);
  createFile(runner, names[1], 2);
  createFile(runner, names[2], 5);
  files.setFileNames(names);

  sorted.clear();
  files.getSortedByName(sorted);
  ezassert( sorted.size() == 3 );
  ezassert( sorted[0]->filename.compare(names[1]) == 0 );
  ezassert( sorted[1]->filename.compare(names[2]) == 0 );
  ezassert( sorted[2]->filename.compare(names[0]) == 0 );
  
  files.clear();
  
  if (runner.clean) {
    ezassert( remove(names[0].c_str()) == 0 );
    ezassert( remove(names[1].c_str()) == 0 );
    ezassert( remove(names[2].c_str()) == 0 );
  }
  
  return true;
}

bool test_getSortedByDateTime(ez::ezTestRunner& runner) {
  ez::ezNcFiles files;
  std::vector<std::string> names;
  std::vector<ez::ezNc*> sorted;
  
  names.resize(3);
  names[0] = createFilename(runner, "2nd");
  names[1] = createFilename(runner, "andlast");
  names[2] = createFilename(runner, "first");
  
  createFile(runner, names[0], 2);
  createFile(runner, names[1], 5);
  createFile(runner, names[2], 0);

  files.setFileNames(names);
  files.load();
  files.getSortedByDatetime(sorted);
  
  ezassert( sorted.size() == 3 );
  ezassert( sorted[0]->filename.compare(names[2]) == 0 );
  ezassert( sorted[1]->filename.compare(names[0]) == 0 );
  ezassert( sorted[2]->filename.compare(names[1]) == 0 );
  
  if (runner.clean) {
    ezassert( remove(names[0].c_str()) == 0 );
    ezassert( remove(names[1].c_str()) == 0 );
    ezassert( remove(names[2].c_str()) == 0 );
  }

  return true;
}

int main(int argc, const char* argv[]) {
  ez::ezTestRunner runner;
  runner.getopt(argc, argv);
  #define TEST(N) runner.add(&test_##N, "test_"#N);
  TEST(getSortedByName);
  TEST(getSortedByDateTime);
  
  return runner.run();
}