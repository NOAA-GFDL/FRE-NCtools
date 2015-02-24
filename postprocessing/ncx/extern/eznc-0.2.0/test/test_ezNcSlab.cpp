/*
20111111 rsz Created.
*/

#include "ezTest.hpp"
#include "ezNcSlab.hpp"

bool test_countEqual(ez::ezTestRunner& runner) {
  ez::ezNcSlab<int> s1, s2;
  
  ezassert( !s1.countEqual(0) );
  ezassert( s1.countEqual(&s2) );
  
  s1.count.push_back(1);
  s1.count.push_back(2);
  s1.count.push_back(3);
  s2.count.push_back(1);
  s2.count.push_back(2);

  ezassert( !s1.countEqual(&s2) );

  s2.count.push_back(3);
  ezassert( s1.countEqual(&s2) );

  s2.count[2] = 13;
  ezassert( !s1.countEqual(&s2) );

  return true;
}

bool test_clear(ez::ezTestRunner& runner) {
  ez::ezNcSlab<int> s1;
  
  s1.count.push_back(1);
  s1.start.push_back(1);
  s1.clear();
  
  ezassert( s1.count.empty() );
  ezassert( s1.start.empty() );
  ezassert( s1.var == 0 );
  
  return true;
}

bool test_setVar(ez::ezTestRunner& runner) {
  ez::ezNcSlab<int> s1;
  ez::ezNcVar v;
  
  s1.setVar(&v);
  
  ezassert( s1.count.empty() );
  ezassert( s1.start.empty() );
  
  v.dimids.push_back(1);
  v.dimids.push_back(2);

  s1.setVar(&v);
  
  ezassert( s1.count.size() == 2 );
  ezassert( s1.start.size() == 2 );
  ezassert( s1.count[1] == 0 );
  ezassert( s1.start[1] == 0 );
  
  return true;
}

int main(int argc, const char* argv[]) {
  ez::ezTestRunner runner;
  runner.getopt(argc, argv);
  #define TEST(N) runner.add(&test_##N, "test_"#N);
  TEST(countEqual);
  TEST(clear);
  TEST(setVar);
  
  return runner.run();
}