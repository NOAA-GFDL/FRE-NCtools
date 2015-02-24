/*
20111111 rsz Created.
*/

#include "ezTest.hpp"
#include "ezNcVariant.hpp"

bool test_sizes(ez::ezTestRunner& runner) {
  ezassert( sizeof(char) == 1 );
  ezassert( sizeof(unsigned char) == 1 );
  ezassert( sizeof(short) == 2 );
  ezassert( sizeof(unsigned short) == 2 );
  ezassert( sizeof(int) == 4 );
  ezassert( sizeof(unsigned int) == 4 );
  ezassert( sizeof(long long) == 8 );
  ezassert( sizeof(unsigned long long) == 8 );
  ezassert( sizeof(float) == 4 );
  ezassert( sizeof(double) == 8 );
  
  return true;
}

bool test_get(ez::ezTestRunner& runner) {
  ez::ezNcVariant v;
  
  ezassert( v.get<unsigned char>() == 0 );
  ezassert( v.get<char>() == 0 );
  ezassert( v.get<short>() == 0 );
  ezassert( v.get<unsigned short>() == 0 );
  ezassert( v.get<int>() == 0 );
  ezassert( v.get<unsigned int>() == 0 );
  ezassert( v.get<long long>() == 0 );
  ezassert( v.get<unsigned long long>() == 0 );
  ezassert( v.get<float>() == 0.0f );
  ezassert( v.get<double>() == 0.0 );

  v.set<int>(100);
  ezassert( v.get<unsigned char>() == 100 );
  ezassert( v.get<char>() == 100 );
  ezassert( v.get<short>() == 100 );
  ezassert( v.get<unsigned short>() == 100 );
  ezassert( v.get<int>() == 100 );
  ezassert( v.get<unsigned int>() == 100 );
  ezassert( v.get<long long>() == 100 );
  ezassert( v.get<unsigned long long>() == 100 );
  ezassert( v.get<float>() == 100.0f );
  ezassert( v.get<double>() == 100.0 );
  
  return true;
}

bool test_set(ez::ezTestRunner& runner) {
  ez::ezNcVariant v;
  
  v.type = NC_UBYTE;
  v.set<char>(101);
  ezassert( v.uc == 101 );

  v.type = NC_BYTE;
  v.set<char>(102);
  ezassert( v.uc == 102 );
  
  v.type = NC_SHORT;
  v.set<char>(103);
  ezassert( v.s == 103 );

  v.type = NC_USHORT;
  v.set<char>(104);
  ezassert( v.us == 104 );

  v.type = NC_INT;
  v.set<char>(105);
  ezassert( v.i == 105 );

  v.type = NC_UINT;
  v.set<char>(106);
  ezassert( v.ui == 106 );

  v.type = NC_INT64;
  v.set<char>(107);
  ezassert( v.l == 107 );

  v.type = NC_UINT64;
  v.set<char>(108);
  ezassert( v.ul == 108 );

  v.type = NC_FLOAT;
  v.set<char>(109);
  ezassert( v.f == 109 );

  v.type = NC_DOUBLE;
  v.set<char>(110);
  ezassert( v.d == 110 );

  v.type = NC_UBYTE;
  v.set<double>(101.);
  ezassert( v.uc == 101 );

  v.type = NC_BYTE;
  v.set<double>(102.);
  ezassert( v.uc == 102 );
  
  v.type = NC_SHORT;
  v.set<double>(103.);
  ezassert( v.s == 103 );

  v.type = NC_USHORT;
  v.set<double>(104.);
  ezassert( v.us == 104 );

  v.type = NC_INT;
  v.set<double>(105.);
  ezassert( v.i == 105 );

  v.type = NC_UINT;
  v.set<double>(106.);
  ezassert( v.ui == 106 );

  v.type = NC_INT64;
  v.set<double>(107.);
  ezassert( v.l == 107 );

  v.type = NC_UINT64;
  v.set<double>(108.);
  ezassert( v.ul == 108 );

  v.type = NC_FLOAT;
  v.set<double>(109.);
  ezassert( v.f == 109 );

  v.type = NC_DOUBLE;
  v.set<double>(110.);
  ezassert( v.d == 110 );

  return true;
}

int main(int argc, const char* argv[]) {
  ez::ezTestRunner runner;
  runner.getopt(argc, argv);
  runner.add(&test_sizes, "test_sizes");
  runner.add(&test_get, "test_get");
  runner.add(&test_set, "test_set");
  
  return runner.run();
}