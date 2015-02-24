/*
20111109 rsz Created. All tests pass.
20111124 rsz Added regex. All pass.
*/

#include "ezTest.hpp"
#include "ezNcUtil.hpp"

bool test_DirName(ez::ezTestRunner& runner) {
  std::string path, dir;
  
  path = "";
  ez::DirName(path, dir);
  ezassert( dir.empty() );
  
  path = "/";
  ez::DirName(path, dir);
  ezassert( dir.compare("/") == 0 );
  
  path = "./";
  ez::DirName(path, dir);
  ezassert( dir.compare("./") == 0 );
  
  path = "./.";
  ez::DirName(path, dir);
  ezassert( dir.compare("./") == 0 );
  
  path = "./a";
  ez::DirName(path, dir);
  ezassert( dir.compare("./") == 0 );

  path = "./asdf";
  ez::DirName(path, dir);
  ezassert( dir.compare("./") == 0 );

  path = "/asdf";
  ez::DirName(path, dir);
  ezassert( dir.compare("/") == 0 );

  path = "/asdf/jkl";
  ez::DirName(path, dir);
  ezassert( dir.compare("/asdf/") == 0 );

  path = "c:\\";
  ez::DirName(path, dir);
  ezassert( dir.compare("c:\\") == 0 );

  path = "c:\\some\\path\\file";
  ez::DirName(path, dir);
  ezassert( dir.compare("c:\\some\\path\\") == 0 );

  return true;
}

bool test_FileExists(ez::ezTestRunner& runner) {
  ezassert( ez::FileExists("test_ezNcUtil") );
  ezassert( ! ez::FileExists("XYZ") );
  return true;
}

bool test_isRevSeq(ez::ezTestRunner& runner) {
  std::vector<int> list;

  list.clear();
  list.push_back(0);

  ezassert( ez::isRevSeq<int>(list) );
    
  list.clear();
  list.push_back(1);
  list.push_back(2);
  list.push_back(3);
  
  ezassert( ! ez::isRevSeq<int>(list) );
  
  list.clear();
  list.push_back(1);
  list.push_back(3);
  
  ezassert( ! ez::isRevSeq<int>(list) );

  list.clear();
  list.push_back(2);
  list.push_back(1);
  list.push_back(0);

  ezassert( ez::isRevSeq<int>(list) );

  list.clear();
  list.push_back(-2);
  list.push_back(-1);
  list.push_back(0);

  ezassert( ! ez::isRevSeq<int>(list) );
  
  list.clear();
  list.push_back(-2);
  list.push_back(-1);
  list.push_back(1);

  ezassert( ! ez::isRevSeq<int>(list) );

  return true;
}

bool test_isSeq(ez::ezTestRunner& runner) {
  std::vector<int> list;

  list.clear();
  list.push_back(0);

  ezassert( ez::isSeq<int>(list) );
    
  list.clear();
  list.push_back(1);
  list.push_back(2);
  list.push_back(3);
  
  ezassert( ez::isSeq<int>(list) );
  
  list.clear();
  list.push_back(1);
  list.push_back(3);
  
  ezassert( ! ez::isSeq<int>(list) );

  list.clear();
  list.push_back(2);
  list.push_back(1);
  list.push_back(0);

  ezassert( ! ez::isSeq<int>(list) );

  list.clear();
  list.push_back(-2);
  list.push_back(-1);
  list.push_back(0);

  ezassert( ez::isSeq<int>(list) );
  
  list.clear();
  list.push_back(-2);
  list.push_back(-1);
  list.push_back(1);

  ezassert( ! ez::isSeq<int>(list) );

  return true;
}

bool test_isStrideSeq(ez::ezTestRunner& runner) {
  std::vector<int> list;

  list.clear();
  list.push_back(0);

  ezassert( ez::isStrideSeq<int>(list) );
    
  list.clear();
  list.push_back(1);
  list.push_back(2);
  list.push_back(3);
  
  ezassert( ez::isStrideSeq<int>(list) );
  
  list.clear();
  list.push_back(1);
  list.push_back(3);
  list.push_back(5);
  
  ezassert( ez::isStrideSeq<int>(list) );

  list.clear();
  list.push_back(22);
  list.push_back(12);
  list.push_back(2);

  ezassert( ez::isStrideSeq<int>(list) );

  list.clear();
  list.push_back(-20);
  list.push_back(-10);
  list.push_back(0);

  ezassert( ez::isStrideSeq<int>(list) );
  
  list.clear();
  list.push_back(-2);
  list.push_back(-1);
  list.push_back(1);

  ezassert( ! ez::isStrideSeq<int>(list) );

  list.clear();
  list.push_back(20);
  list.push_back(10);
  list.push_back(1);

  ezassert( ! ez::isStrideSeq<int>(list) );

  std::vector<double> dlist;
  dlist.push_back(100.5);
  dlist.push_back(200.5);
  dlist.push_back(300.5);
  
  ezassert( ez::isStrideSeq<double>(dlist) );

  dlist[2] = 300.1;
  ezassert( ! ez::isStrideSeq<double>(dlist) );

  return true;
}

bool test_process_mem_usage(ez::ezTestRunner& runner) {
  double vm, rss;
  ez::process_mem_usage(vm, rss);
  ezassert( vm > 0 );
  ezassert( rss > 0 );
  
  return true;
}

bool test_osQueryPerformance(ez::ezTestRunner& runner) {
  long long tic, toc;
  tic = ez::osQueryPerformance();
  
  // Make fake loop to cause time delay.
  std::vector<double> res;
  for(int i=0; i < 1000000; ++i) res.push_back(i/(i+1));
  res.resize(0);
  
  toc = ez::osQueryPerformance();
  ezassert( (toc - tic) > 0 );
  
  return true;
}

bool test_seq(ez::ezTestRunner& runner) {
  std::vector<int> list;
  ez::seq<int>(list, 0);
  ezassert( list.size() == 0 );

  ez::seq<int>(list, 1);
  ezassert( list.size() == 1 );
  ezassert( list[0] == 0 );

  list.clear();
  ez::seq<int>(list, 10);
  ezassert( list.size() == 10 );
  ezassert( list[9] == 9 );

  return true;
}

bool test_seqPOD(ez::ezTestRunner& runner) {
  int* null=0;
  int list[10];
  int garbage = list[0];
  
  ez::seq<int>(null, 0); // Should not segfault.
  ez::seq<int>(null, 10); // Should not segfault.

  ez::seq<int>(list, 0);
  ezassert( garbage == list[0] );
  
  garbage = list[1];
  ez::seq<int>(list, 1);
  ezassert( list[0] == 0 );
  ezassert( garbage == list[1] );

  ez::seq<int>(list, 10);
  ezassert( list[0] == 0 );
  ezassert( list[5] == 5 );
  ezassert( list[9] == 9 );

  return true;
}

bool test_seqMinMax(ez::ezTestRunner& runner) {
  std::vector<int> list;
  ez::seq<int>(list, 0, 0);
  ezassert( list.size() == 1 );
  ezassert( list[0] == 0 );

  list.clear();
  ez::seq<int>(list, 0, 1);
  ezassert( list.size() == 2 );
  ezassert( list[0] == 0 );
  ezassert( list[1] == 1 );

  list.clear();
  ez::seq<int>(list, 0, -1);
  ezassert( list.size() == 2 );
  ezassert( list[0] == 0 );
  ezassert( list[1] == -1 );

  list.clear();
  ez::seq<int>(list, 0, 10);
  ezassert( list.size() == 11 );
  ezassert( list[0] == 0 );
  ezassert( list[10] == 10 );

  list.clear();
  ez::seq<int>(list, -10, 0);
  ezassert( list.size() == 11 );
  ezassert( list[0] == -10 );
  ezassert( list[10] == 0 );

  return true;
}

bool test_seqMinMaxPOD(ez::ezTestRunner& runner) {
  int list[20];
  int *null=0;
  int garbage;
  
  ez::seq<int>(null, 0, 0);
  ez::seq<int>(null, 0, 10);
  ez::seq<int>(null, 10, 0);
  
  garbage = list[1];
  ez::seq<int>(list, 0, 0);
  ezassert( list[0] == 0 );
  ezassert( list[1] == garbage );

  garbage = list[2];
  ez::seq<int>(list, 0, 1);
  ezassert( list[0] == 0 );
  ezassert( list[1] == 1 );
  ezassert( list[2] == garbage );

  ez::seq<int>(list, 0, -1);
  ezassert( list[0] == 0 );
  ezassert( list[1] == -1 );

  ez::seq<int>(list, 0, 10);
  ezassert( list[0] == 0 );
  ezassert( list[10] == 10 );

  ez::seq<int>(list, -10, 0);
  ezassert( list[0] == -10 );
  ezassert( list[10] == 0 );

  return true;
}

bool test_StringToNcType(ez::ezTestRunner& runner) {
  ezassert( ez::StringToNcType("byte") == NC_BYTE );
  ezassert( ez::StringToNcType("ubyte") == NC_UBYTE );
  ezassert( ez::StringToNcType("char") == NC_CHAR );
  ezassert( ez::StringToNcType("uchar") == NC_UBYTE );
  ezassert( ez::StringToNcType("short") == NC_SHORT );
  ezassert( ez::StringToNcType("ushort") == NC_USHORT );
  ezassert( ez::StringToNcType("int") == NC_INT );
  ezassert( ez::StringToNcType("uint") == NC_UINT );
  ezassert( ez::StringToNcType("int64") == NC_INT64 );
  ezassert( ez::StringToNcType("uint64") == NC_UINT64 );
  ezassert( ez::StringToNcType("long") == NC_INT );
  ezassert( ez::StringToNcType("float") == NC_FLOAT );
  ezassert( ez::StringToNcType("double") == NC_DOUBLE );
  ezassert( ez::StringToNcType("string") == NC_STRING );
  ezassert( ez::StringToNcType("nibble") == NC_NAT );
  
  return true;
}

bool test_StringsToValues(ez::ezTestRunner& runner) {
  std::vector<std::string> strings;
  std::vector<int> values;
  
  ez::StringsToValues<int>(strings, 0, -1, values);
  ezassert( values.empty() );
  
  strings.clear();
  strings.push_back("0");
  strings.push_back("1");
  ez::StringsToValues<int>(strings, 0, 1, values);
  ezassert( values.size() == 2 );
  ezassert( values[0] == 0 );
  ezassert( values[1] == 1 );

  strings.clear();
  strings.push_back("0");
  strings.push_back("1");
  strings.push_back("11");
  ez::StringsToValues<int>(strings, 0, 2, values);
  ezassert( values.size() == 3 );
  ezassert( values[0] == 0 );
  ezassert( values[2] == 11 );

  strings.clear();
  strings.push_back("-1234");
  strings.push_back("-11");
  strings.push_back("-22");
  ez::StringsToValues<int>(strings, 0, 0, values);
  ezassert( values.size() == 1 );
  ezassert( values[0] == -1234 );

  strings.clear();
  strings.push_back("0");
  strings.push_back("-11");
  strings.push_back("-22");
  ez::StringsToValues<int>(strings, 1, 2, values);
  ezassert( values.size() == 2 );
  ezassert( values[0] == -11 );
  ezassert( values[1] == -22 );

  strings.clear();
  strings.push_back("0");
  strings.push_back("111");
  strings.push_back("-222");
  ez::StringsToValues<int>(strings, 1, 1, values);
  ezassert( values.size() == 1 );
  ezassert( values[0] == 111 );
  
  return true;
}

bool test_Regex(ez::ezTestRunner& runner) {
  std::vector<std::string> strings, matches;
  std::string pattern;
  
  // Empty case.
  ezassert( ez::regex(pattern, strings, matches, false, false) == 0 );
  ezassert( matches.empty() );

  strings.push_back("fdsa");
  strings.push_back("ASDF");
  strings.push_back("a,./;");
  strings.push_back("asdf");

  pattern = "a";
  matches.clear();
  ezassert( ez::regex(pattern, strings, matches, false, false) == 0 );
  ezassert( matches.size() == 3 );

  pattern = "asdf";
  matches.clear();
  ezassert( ez::regex(pattern, strings, matches, false, false) == 0 );
  ezassert( matches.size() == 1 );
  ezassert( matches[0].compare("asdf") == 0 );
  
  pattern = "^a[[:alnum:]]";
  matches.clear();
  ezassert( ez::regex(pattern, strings, matches, false, false) == 0 );
  ezassert( matches.size() == 1 );
  ezassert( matches[0].compare("asdf") == 0 );

  // Ignore case.
  pattern = "^a[[:alnum:]]";
  matches.clear();
  ezassert( ez::regex(pattern, strings, matches, false, true) == 0 );
  ezassert( matches.size() == 2 );
  ezassert( matches[0].compare("ASDF") == 0 );
  ezassert( matches[1].compare("asdf") == 0 );
  
  // Fail union without extended.
  pattern = "asdf|fdsa";
  matches.clear();
  ezassert( ez::regex(pattern, strings, matches, false, false) == 0 );
  ezassert( matches.empty() );

  // Extended supports union.
  pattern = "asdf|ASDF";
  matches.clear();
  ezassert( ez::regex(pattern, strings, matches, true, false) == 0 );
  ezassert( matches.size() == 2 );
  ezassert( matches[0].compare("ASDF") == 0 );
  ezassert( matches[1].compare("asdf") == 0 );

  // Extended, no case.
  pattern = "asdf|fdsa";
  matches.clear();
  ezassert( ez::regex(pattern, strings, matches, true, true) == 0 );
  ezassert( matches.size() == 3 );
  ezassert( matches[0].compare("fdsa") == 0 );
  ezassert( matches[1].compare("ASDF") == 0 );
  ezassert( matches[2].compare("asdf") == 0 );
  
  return true;
}

bool test_TmpFilename(ez::ezTestRunner& runner) {
  std::string string;

  ez::TmpFilename(string, 0);
  ezassert( string.size() > 0 );
  ezassert( string.compare(0,4,"tmp.") == 0 );

  string.clear();
  ez::TmpFilename(string, ".ext");
  ezassert( string.size() > 0 );
  ezassert( string.compare(10,4,".ext") == 0 );
  
  string = "/tmp";
  ez::TmpFilename(string, "");
  ezassert( string.size() == (4+11) );
  ezassert( string.compare(0,9,"/tmp/tmp.") == 0 );

  string = "/tmp";
  ez::TmpFilename(string, ".ext");
  ezassert( string.size() == (5+10+4) );
  ezassert( string.compare(15,9,".ext") == 0 );
  
  return true;
}

int main(int argc, const char* argv[]) {
  ez::ezTestRunner runner;
  runner.getopt(argc, argv);
  runner.add(&test_DirName, "test_DirName");
  runner.add(&test_FileExists, "test_FileExists");
  runner.add(&test_isRevSeq, "test_isRevSeq");
  runner.add(&test_isSeq, "test_isSeq");
  runner.add(&test_isStrideSeq, "test_isStrideSeq");
  runner.add(&test_process_mem_usage, "test_process_mem_usage");
  runner.add(&test_osQueryPerformance, "test_osQueryPerformance");
  runner.add(&test_Regex, "test_Regex");
  runner.add(&test_seq, "test_seq");
  runner.add(&test_seqPOD, "test_seqPOD");
  runner.add(&test_seqMinMax, "test_seqMinMax");
  runner.add(&test_seqMinMaxPOD, "test_seqMinMaxPOD");
  runner.add(&test_StringToNcType, "test_StringToNcType");
  runner.add(&test_StringsToValues, "test_StringsToValues");
  runner.add(&test_TmpFilename, "test_TmpFilename");
  
  return runner.run();
}