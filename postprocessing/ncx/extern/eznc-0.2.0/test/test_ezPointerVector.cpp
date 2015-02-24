/*
20111111 rsz Created.
*/

#include "ezTest.hpp"
#include "ezPointerVector.hpp"

bool test_reset(ez::ezTestRunner& runner) {
  ez::ezPointerVector< std::vector<int> > v;
  
  for(int i=0; i < 10; ++i) {
    v.vector.push_back( new std::vector<int> );
    v.vector[i]->resize(100);
  }
  
  v.reset();
  ezassert( v.vector.empty() );
  
  return true;
}

int main(int argc, const char* argv[]) {
  ez::ezTestRunner runner;
  runner.getopt(argc, argv);
  runner.add(&test_reset, "test_reset");
  
  return runner.run();
}