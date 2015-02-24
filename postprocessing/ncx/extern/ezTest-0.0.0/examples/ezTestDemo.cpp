/*
Demo of using ezTestRunner class.
./ezTestDemo --verbose --color

20111107 rsz Created.
*/

#include "ezTest.hpp"

bool test_ezassert(ez::ezTestRunner& runner) {
  ezassert(runner.clean);
  return true;
}
bool test_clean(ez::ezTestRunner& runner) {
  return runner.clean;
}
bool test_tmpdir(ez::ezTestRunner& runner) {
  return runner.tmpdir.compare("/tmp")==0;
}
bool test_verbose(ez::ezTestRunner& runner) {
  return runner.verbose;
}

int main(int argc, const char* argv[]) {
  ez::ezTestRunner t;
  t.getopt(argc, argv);
  t.add(&test_ezassert, "test_EzAssert. See that clean is on by using ezassert test.");
  t.add(&test_clean, "test_Clean. See if --clean is on.");
  t.add(&test_tmpdir, "test_TmpDir. See if --tmpdir is set to \"/tmp\".");
  t.add(&test_verbose, "test_Verbose. See if --verbose was used.");

  // Do setup here and/or in functions.
  
  int status = t.run();
  
  // Do tear down here and/or in functions.
  
  return status;
}