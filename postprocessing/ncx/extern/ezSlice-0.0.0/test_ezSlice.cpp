/*
20111031 rsz Created. Tests pass.
*/
#include "ezSlice.hpp"
#include "ezTest.hpp"

void print(int mn, int mx, int st) {
  std::cout << "[ ";
  
  if (st > 0) {
    for(; mn <= mx; mn += st)
      std::cout << mn << " ";
  } else {
    for(; mn >= mx; mn += st)
      std::cout << mn << " ";  
  }
  
  std::cout << "]\n";
}
void print(std::vector<int> & list) {
  std::cout << "[ ";
  for(int i=0; i < list.size(); ++i)
    std::cout << list[i] << " ";
    
  std::cout << "]\n";
}

bool test_setCharPtrs(ez::ezTestRunner& runner) {
  ez::ezSlice<int> sli;
  int min,max,stride,n;
  std::vector<int> idx;
  n = 10;
  
  if (runner.verbose) std::cout << "in: 2 0 0\n";
  sli.set("2",0,0);
  sli.get(min,max,stride,n);
  ezassert( min == 2 );
  ezassert( max == 9 );
  ezassert( stride == 1 );
  ezassert( sli.size(n) == 8 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);
  
  if (runner.verbose) std::cout << "\nin: 2 5 0\n";
  sli.set("2","5",0);
  sli.get(min,max,stride,n);
  ezassert( min == 2 );
  ezassert( max == 5 );
  ezassert( stride == 1 );
  ezassert( sli.size(n) == 4 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);

  if (runner.verbose) std::cout << "\nin: 2 7 2\n";
  sli.set("2","7","2");
  sli.get(min,max,stride,n);
  ezassert( min == 2 );
  ezassert( max == 7 );
  ezassert( stride == 2 );
  ezassert( sli.size(n) == 3 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);

  if (runner.verbose) std::cout << "\nin: 0 5 0\n";
  sli.set(0,"5",0);
  sli.get(min,max,stride,n);
  ezassert( min == 0 );
  ezassert( max == 5 );
  ezassert( stride == 1 );
  ezassert( sli.size(n) == 6 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);

  if (runner.verbose) std::cout << "\nin: 0 5 3\n";
  sli.set(0,"5","3");
  sli.get(min,max,stride,n);
  ezassert( min == 0 );
  ezassert( max == 5 );
  ezassert( stride == 3 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);
  ezassert( sli.size(n) == 2 );

  if (runner.verbose) std::cout << "\nin: 0 5 -1\n";
  sli.set(0,"5","-1");
  sli.get(min,max,stride,n);
  ezassert( min == 9 );
  ezassert( max == 5 );
  ezassert( stride == -1 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);
  ezassert( sli.size(n) == 5 );

  if (runner.verbose) std::cout << "\nin: -9 0 0\n";
  sli.set("-9",0,0);
  sli.get(min,max,stride,n);
  ezassert( min == 1 );
  ezassert( max == 9 );
  ezassert( stride == 1 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);
  ezassert( sli.size(n) == 9 );

  if (runner.verbose) std::cout << "\nin: -9 5 0\n";
  sli.set("-9","5",0);
  sli.get(min,max,stride,n);
  ezassert( min == 1 );
  ezassert( max == 5 );
  ezassert( stride == 1 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);
  ezassert( sli.size(n) == 5 );

  if (runner.verbose) std::cout << "\nin: -9 -6 0\n";
  sli.set("-9","-6",0);
  sli.get(min,max,stride,n);
  ezassert( min == 1 );
  ezassert( max == 4 );
  ezassert( stride == 1 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);
  ezassert( sli.size(n) == 4 );

  if (runner.verbose) std::cout << "\nin: -9 -6 2\n";
  sli.set("-9","-6","2");
  sli.get(min,max,stride,n);
  ezassert( min == 1 );
  ezassert( max == 4 );
  ezassert( stride == 2 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);
  ezassert( sli.size(n) == 2 );

  if (runner.verbose) std::cout << "\nin: -1 -6 -2\n";
  sli.set("-1","-6","-2");
  sli.get(min,max,stride,n);
  ezassert( min == 9 );
  ezassert( max == 4 );
  ezassert( stride == -2 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);
  ezassert( sli.size(n) == 3 );

  if (runner.verbose) std::cout << "\nin: 0 4 2\n";
  sli.set("0","4","2");
  sli.get(min,max,stride,n);
  ezassert( min == 0 );
  ezassert( max == 4 );
  ezassert( stride == 2 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);
  ezassert( sli.size(n) == 3 );

  return true;
}

bool test_setString(ez::ezTestRunner& runner) {
  ez::ezSlice<int> sli;
  int min,max,stride,n;
  std::vector<int> idx;
  n = 10;
  std::string str;

  str = "2:";
  if (runner.verbose) std::cout << "\nin: " << str << "\n";
  sli.set(str,':');
  sli.get(min,max,stride,n);
  ezassert( min == 2 );
  ezassert( max == 9 );
  ezassert( stride == 1 );
  sli.get(idx,n);
  ezassert( idx.size() == 8);
  ezassert( idx[0] == 2 );
  ezassert( idx[7] == 9 );
  if (runner.verbose) print(idx);

  str = "2:5";
  if (runner.verbose) std::cout << "\nin: " << str << "\n";
  sli.set(str,':');
  sli.get(min,max,stride,n);
  ezassert( min == 2 );
  ezassert( max == 5 );
  ezassert( stride == 1 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);
  ezassert( idx.size() == 4);
  ezassert( idx[0] == 2 );
  ezassert( idx[3] == 5 );

  str = "2:7:2";
  if (runner.verbose) std::cout << "\nin: " << str << "\n";
  sli.set(str,':');
  sli.get(min,max,stride,n);
  ezassert( min == 2 );
  ezassert( max == 7 );
  ezassert( stride == 2 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);
  ezassert( idx.size() == 3);
  ezassert( idx[0] == 2 );
  ezassert( idx[2] == 6 );

  str = ":5";
  if (runner.verbose) std::cout << "\nin: " << str << "\n";
  sli.set(str,':');
  sli.get(min,max,stride,n);
  ezassert( min == 0 );
  ezassert( max == 5 );
  ezassert( stride == 1 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);
  ezassert( idx.size() == 6);
  ezassert( idx[0] == 0 );
  ezassert( idx[5] == 5 );

  str = ":5:3";
  if (runner.verbose) std::cout << "\nin: " << str << "\n";
  sli.set(str,':');
  sli.get(min,max,stride,n);
  ezassert( min == 0 );
  ezassert( max == 5 );
  ezassert( stride == 3 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);
  ezassert( idx.size() == 2);
  ezassert( idx[0] == 0 );
  ezassert( idx[1] == 3 );

  str = ":5:-1";
  if (runner.verbose) std::cout << "\nin: " << str << "\n";
  sli.set(str,':');
  sli.get(min,max,stride,n);
  ezassert( min == 9 );
  ezassert( max == 5 );
  ezassert( stride == -1 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);
  ezassert( idx.size() == 5);
  ezassert( idx[0] == 9 );
  ezassert( idx[4] == 5 );

  str = "-9:";
  if (runner.verbose) std::cout << "\nin: " << str << "\n";
  sli.set(str,':');
  sli.get(min,max,stride,n);
  ezassert( min == 1 );
  ezassert( max == 9 );
  ezassert( stride == 1 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);
  ezassert( idx.size() == 9);
  ezassert( idx[0] == 1 );
  ezassert( idx[8] == 9 );

  str = "-9:5";
  if (runner.verbose) std::cout << "\nin: " << str << "\n";
  sli.set(str,':');
  sli.get(min,max,stride,n);
  ezassert( min == 1 );
  ezassert( max == 5 );
  ezassert( stride == 1 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);
  ezassert( idx.size() == 5);
  ezassert( idx[0] == 1 );
  ezassert( idx[4] == 5 );

  str = "-9:-6";
  if (runner.verbose) std::cout << "\nin: " << str << "\n";
  sli.set(str,':');
  sli.get(min,max,stride,n);
  ezassert( min == 1 );
  ezassert( max == 4 );
  ezassert( stride == 1 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);
  ezassert( idx.size() == 4);
  ezassert( idx[0] == 1 );
  ezassert( idx[3] == 4 );

  str = "-9:-6:2";
  if (runner.verbose) std::cout << "\nin: " << str << "\n";
  sli.set(str,':');
  sli.get(min,max,stride,n);
  ezassert( min == 1 );
  ezassert( max == 4 );
  ezassert( stride == 2 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);
  ezassert( idx.size() == 2);
  ezassert( idx[0] == 1 );
  ezassert( idx[1] == 3 );

  str = "-1:-6:-2";
  if (runner.verbose) std::cout << "\nin: " << str << "\n";
  sli.set(str,':');
  sli.get(min,max,stride,n);
  ezassert( min == 9 );
  ezassert( max == 4 );
  ezassert( stride == -2 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);
  ezassert( idx.size() == 3);
  ezassert( idx[0] == 9 );
  ezassert( idx[2] == 5 );

  str = "::2";
  if (runner.verbose) std::cout << "\nin: " << str << "\n";
  sli.set(str,':');
  sli.get(min,max,stride,n);
  ezassert( min == 0 );
  ezassert( max == 9 );
  ezassert( stride == 2 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);
  ezassert( idx.size() == 5);
  ezassert( idx[0] == 0 );
  ezassert( idx[4] == 8 );

  str = "::-3";
  if (runner.verbose) std::cout << "\nin: " << str << "\n";
  sli.set(str,':');
  sli.get(min,max,stride,n);
  ezassert( min == 9 );
  ezassert( max == 0 );
  ezassert( stride == -3 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);
  ezassert( idx.size() == 4);
  ezassert( idx[0] == 9 );
  ezassert( idx[3] == 0 );

  str = "...";
  if (runner.verbose) std::cout << "\nin: " << str << "\n";
  sli.set(str,':');
  sli.get(min,max,stride,n);
  ezassert( min == 0 );
  ezassert( max == 9 );
  ezassert( stride == 1 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);
  ezassert( idx.size() == 10);
  ezassert( idx[0] == 0 );
  ezassert( idx[9] == 9 );

  str = ":";
  if (runner.verbose) std::cout << "\nin: " << str << "\n";
  sli.set(str,':');
  sli.get(min,max,stride,n);
  ezassert( min == 0 );
  ezassert( max == 9 );
  ezassert( stride == 1 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);
  ezassert( idx.size() == 10);
  ezassert( idx[0] == 0 );
  ezassert( idx[9] == 9 );

  str = "::";
  if (runner.verbose) std::cout << "\nin: " << str << "\n";
  sli.set(str,':');
  sli.get(min,max,stride,n);
  ezassert( min == 0 );
  ezassert( max == 9 );
  ezassert( stride == 1 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);
  ezassert( idx.size() == 10);
  ezassert( idx[0] == 0 );
  ezassert( idx[9] == 9 );

  str = "3";
  if (runner.verbose) std::cout << "\nin: " << str << "\n";
  sli.set(str,':');
  sli.get(min,max,stride,n);
  ezassert( min == 3 );
  ezassert( max == 3 );
  ezassert( stride == 1 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);
  ezassert( idx.size() == 1);
  ezassert( idx[0] == 3 );

  return true;
}

bool test_offset(ez::ezTestRunner& runner) {
  ez::ezSlice<int> sli;
  int min,max,stride,n;
  std::vector<int> idx;
  n = 10;
  std::string str;

  str = "2:";
  if (runner.verbose) std::cout << "\nin: " << str << "\n";
  sli.set(str,':',-1);
  sli.get(min,max,stride,n);
  ezassert( min == 1 );
  ezassert( max == 9 );
  ezassert( stride == 1 );
  sli.get(idx,n);
  ezassert( idx.size() == 9);
  ezassert( idx[0] == 1 );
  ezassert( idx[8] == 9 );
  if (runner.verbose) print(idx);

  str = "2:5";
  if (runner.verbose) std::cout << "\nin: " << str << "\n";
  sli.set(str,':',-1);
  sli.get(min,max,stride,n);
  ezassert( min == 1 );
  ezassert( max == 4 );
  ezassert( stride == 1 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);
  ezassert( idx.size() == 4);
  ezassert( idx[0] == 1 );
  ezassert( idx[3] == 4 );

  str = "2:7:2";
  if (runner.verbose) std::cout << "\nin: " << str << "\n";
  sli.set(str,':',-1);
  sli.get(min,max,stride,n);
  ezassert( min == 1 );
  ezassert( max == 6 );
  ezassert( stride == 2 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);
  ezassert( idx.size() == 3);
  ezassert( idx[0] == 1 );
  ezassert( idx[2] == 5 );
  
  str = ":5";
  if (runner.verbose) std::cout << "\nin: " << str << "\n";
  sli.set(str,':',-1);
  sli.get(min,max,stride,n);
  ezassert( min == 0 );
  ezassert( max == 4 );
  ezassert( stride == 1 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);
  ezassert( idx.size() == 5);
  ezassert( idx[0] == 0 );
  ezassert( idx[4] == 4 );
  
  str = ":5:3";
  if (runner.verbose) std::cout << "\nin: " << str << "\n";
  sli.set(str,':',-1);
  sli.get(min,max,stride,n);
  ezassert( min == 0 );
  ezassert( max == 4 );
  ezassert( stride == 3 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);
  ezassert( idx.size() == 2);
  ezassert( idx[0] == 0 );
  ezassert( idx[1] == 3 );

  str = ":5:-1";
  if (runner.verbose) std::cout << "\nin: " << str << "\n";
  sli.set(str,':',-1);
  sli.get(min,max,stride,n);
  ezassert( min == 9 );
  ezassert( max == 4 );
  ezassert( stride == -1 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);
  ezassert( idx.size() == 6);
  ezassert( idx[0] == 9 );
  ezassert( idx[5] == 4 );

  str = "-9:";
  if (runner.verbose) std::cout << "\nin: " << str << "\n";
  sli.set(str,':',-1);
  sli.get(min,max,stride,n);
  ezassert( min == 1 );
  ezassert( max == 9 );
  ezassert( stride == 1 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);
  ezassert( idx.size() == 9);
  ezassert( idx[0] == 1 );
  ezassert( idx[8] == 9 );

  str = "-9:5";
  if (runner.verbose) std::cout << "\nin: " << str << "\n";
  sli.set(str,':',-1);
  sli.get(min,max,stride,n);
  ezassert( min == 1 );
  ezassert( max == 4 );
  ezassert( stride == 1 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);
  ezassert( idx.size() == 4);
  ezassert( idx[0] == 1 );
  ezassert( idx[3] == 4 );

  str = "-9:-6";
  if (runner.verbose) std::cout << "\nin: " << str << "\n";
  sli.set(str,':',-1);
  sli.get(min,max,stride,n);
  ezassert( min == 1 );
  ezassert( max == 4 );
  ezassert( stride == 1 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);
  ezassert( idx.size() == 4);
  ezassert( idx[0] == 1 );
  ezassert( idx[3] == 4 );

  str = "-9:-6:2";
  if (runner.verbose) std::cout << "\nin: " << str << "\n";
  sli.set(str,':',-1);
  sli.get(min,max,stride,n);
  ezassert( min == 1 );
  ezassert( max == 4 );
  ezassert( stride == 2 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);
  ezassert( idx.size() == 2);
  ezassert( idx[0] == 1 );
  ezassert( idx[1] == 3 );

  str = "-1:-6:-2";
  if (runner.verbose) std::cout << "\nin: " << str << "\n";
  sli.set(str,':',-1);
  sli.get(min,max,stride,n);
  ezassert( min == 9 );
  ezassert( max == 4 );
  ezassert( stride == -2 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);
  ezassert( idx.size() == 3);
  ezassert( idx[0] == 9 );
  ezassert( idx[2] == 5 );

  str = "::2";
  if (runner.verbose) std::cout << "\nin: " << str << "\n";
  sli.set(str,':',-1);
  sli.get(min,max,stride,n);
  ezassert( min == 0 );
  ezassert( max == 9 );
  ezassert( stride == 2 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);
  ezassert( idx.size() == 5);
  ezassert( idx[0] == 0 );
  ezassert( idx[4] == 8 );

  str = "::-3";
  if (runner.verbose) std::cout << "\nin: " << str << "\n";
  sli.set(str,':',-1);
  sli.get(min,max,stride,n);
  ezassert( min == 9 );
  ezassert( max == 0 );
  ezassert( stride == -3 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);
  ezassert( idx.size() == 4);
  ezassert( idx[0] == 9 );
  ezassert( idx[3] == 0 );

  str = "...";
  if (runner.verbose) std::cout << "\nin: " << str << "\n";
  sli.set(str,':',-1);
  sli.get(min,max,stride,n);
  ezassert( min == 0 );
  ezassert( max == 9 );
  ezassert( stride == 1 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);
  ezassert( idx.size() == 10);
  ezassert( idx[0] == 0 );
  ezassert( idx[9] == 9 );

  str = ":";
  if (runner.verbose) std::cout << "\nin: " << str << "\n";
  sli.set(str,':',-1);
  sli.get(min,max,stride,n);
  ezassert( min == 0 );
  ezassert( max == 9 );
  ezassert( stride == 1 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);
  ezassert( idx.size() == 10);
  ezassert( idx[0] == 0 );
  ezassert( idx[9] == 9 );

  str = "::";
  if (runner.verbose) std::cout << "\nin: " << str << "\n";
  sli.set(str,':',-1);
  sli.get(min,max,stride,n);
  ezassert( min == 0 );
  ezassert( max == 9 );
  ezassert( stride == 1 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);
  ezassert( idx.size() == 10);
  ezassert( idx[0] == 0 );
  ezassert( idx[9] == 9 );

  str = "3";
  if (runner.verbose) std::cout << "\nin: " << str << "\n";
  sli.set(str,':',-1);
  sli.get(min,max,stride,n);
  ezassert( min == 2 );
  ezassert( max == 2 );
  ezassert( stride == 1 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);
  ezassert( idx.size() == 1);
  ezassert( idx[0] == 2 );

  str = "1:5:2";
  if (runner.verbose) std::cout << "\nin: " << str << "\n";
  sli.set(str,':',-1);
  sli.get(min,max,stride,n);
  ezassert( min == 0 );
  ezassert( max == 4 );
  ezassert( stride == 2 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);
  ezassert( sli.size(n) == 3 );

  str = "3:7:2";
  if (runner.verbose) std::cout << "\nin: " << str << "\n";
  sli.set(str,':',-1);
  sli.get(min,max,stride,n);
  ezassert( min == 2 );
  ezassert( max == 6 );
  ezassert( stride == 2 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);
  ezassert( sli.size(n) == 3 );

  str = "3:7:2";
  if (runner.verbose) std::cout << "\nin: " << str << "\n";
  sli.set("3", "7", "2", -1);
  sli.get(min,max,stride,n);
  ezassert( min == 2 );
  ezassert( max == 6 );
  ezassert( stride == 2 );
  sli.get(idx,n);
  if (runner.verbose) print(idx);
  ezassert( sli.size(n) == 3 );

  return true;
}
  
int main(int argc, const char* argv[]) {
  ez::ezTestRunner t;
  t.getopt(argc, argv);
  t.add(&test_setCharPtrs, "test_setCharPtrs.");
  t.add(&test_setString, "test_setString.");
  t.add(&test_offset, "test_offset.");

  // Do setup here and/or in functions.
  
  int status = t.run();
  
  // Do tear down here and/or in functions.
  
  return status;
}