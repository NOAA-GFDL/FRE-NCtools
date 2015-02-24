/*
20111119 rsz Created.
*/
#include "ezStringUtil.hpp"
#include "ezTest.hpp"

#define DEBUGLINE() printf("%s:%d\n", __FILE__, __LINE__);

bool test_compare_nocase(ez::ezTestRunner& runner) {
	std::string s1;
	std::string s2;
  
  s2 = "0";
  ezassert( ez::compare_nocase(s1,s2) );

  s1 = "1";
  s2 = "0";
  ezassert( !ez::compare_nocase(s1,s2) );

  s1 = "0";
  s2 = "1";
  ezassert( ez::compare_nocase(s1,s2) );

  s1 = "0";
  s2 = "0";
  ezassert( !ez::compare_nocase(s1,s2) );

  s1 = "a";
  s2 = "A";
  ezassert( !ez::compare_nocase(s1,s2) );

  s1 = "abc";
  s2 = "ABC";
  ezassert( !ez::compare_nocase(s1,s2) );
  
  s1 = "abC";
  s2 = "Abc";
  ezassert( !ez::compare_nocase(s1,s2) );

  s1 = "abC";
  s2 = "Abcd";
  ezassert( ez::compare_nocase(s1,s2) );

  s1 = "abcd";
  s2 = "abc";
  ezassert( !ez::compare_nocase(s1,s2) );

  s1 = "abcd";
  s2 = "aBC";
  ezassert( !ez::compare_nocase(s1,s2) );

  s1 = "a b c";
  s2 = "a B C";
  ezassert( !ez::compare_nocase(s1,s2) );
  
  s1 = "a b c";
  s2 = "a B C ";
  ezassert( ez::compare_nocase(s1,s2) );

  return true;
}

bool test_find_first(ez::ezTestRunner& runner) {
  std::vector<std::string> v;
  std::string s;
  
  ezassert( ez::find_first(v, s) < 0 );
  
  v.push_back("asdf");
  v.push_back("fdsa");
  v.push_back("ASDF");

  ezassert( ez::find_first(v, s) < 0 );
  s = "ASD";
  ezassert( ez::find_first(v, s) < 0 );
  
  s = "ASDFE";
  ezassert( ez::find_first(v, s) < 0 );

  s = "ASDF";
  ezassert( ez::find_first(v, s) == 2 );

  ezassert( ez::find_first(v, "") < 0 );
  ezassert( ez::find_first(v, "fdsa") == 1 );
  ezassert( ez::find_first(v, "ASD") < 0 );
  ezassert( ez::find_first(v, "ASDFE") < 0 );
  ezassert( ez::find_first(v, "ASDF") == 2 );

  return true;
}

bool test_lower(ez::ezTestRunner& runner) {
  std::vector<std::string> v;
  std::string s;

  ez::lower(s);
  ezassert( s.compare("") == 0 );

  s = "asdf";
  ez::lower(s);
  ezassert( s.compare("asdf") == 0 );

  s = "ASDF";
  ez::lower(s);
  ezassert( s.compare("asdf") == 0 );

  s = "asdF_1234_JKL;";
  ez::lower(s);
  ezassert( s.compare("asdf_1234_jkl;") == 0 );

  return true;
}

bool test_not_in_both(ez::ezTestRunner& runner) {
  std::vector<std::string> v1, v2, v3;

  ez::not_in_both(v1, v2, v3);
  ezassert( v3.empty() );
  
  v1.push_back("1");
  v1.push_back("3");
  v1.push_back("5");
  v1.push_back("7");

  ez::not_in_both(v1, v2, v3);
  ezassert( v3.size() == 4 );
  ezassert( v3[0].compare("1") == 0 );
  ezassert( v3[1].compare("3") == 0 );
  ezassert( v3[2].compare("5") == 0 );
  ezassert( v3[3].compare("7") == 0 );

  v2.push_back("1");
  v2.push_back("2");
  v2.push_back("3");
  v2.push_back("4");
  
  ez::not_in_both(v1, v2, v3);
  ezassert( v3.size() == 4 );
  ezassert( v3[0].compare("2") == 0 );
  ezassert( v3[1].compare("4") == 0 );
  ezassert( v3[2].compare("5") == 0 );
  ezassert( v3[3].compare("7") == 0 );
    
  return true;
}

bool test_not_in_second(ez::ezTestRunner& runner) {
  std::vector<std::string> v1, v2, v3;

  ez::not_in_second(v1, v2, v3);
  ezassert( v3.empty() );
  
  v1.push_back("1");
  v1.push_back("3");
  v1.push_back("5");
  v1.push_back("7");

  ez::not_in_second(v1, v2, v3);
  ezassert( v3.size() == 4 );
  ezassert( v3[0].compare("1") == 0 );
  ezassert( v3[1].compare("3") == 0 );
  ezassert( v3[2].compare("5") == 0 );
  ezassert( v3[3].compare("7") == 0 );

  v2.push_back("1");
  v2.push_back("2");
  v2.push_back("3");
  v2.push_back("4");
  
  ez::not_in_second(v1, v2, v3);
  ezassert( v3.size() == 2 );
  ezassert( v3[0].compare("5") == 0 );
  ezassert( v3[1].compare("7") == 0 );

  return true;
}

bool test_intersection(ez::ezTestRunner& runner) {
  std::vector<std::string> v1, v2, v3;

  ez::intersection(v1, v2, v3);
  ezassert( v3.empty() );
  
  v1.push_back("1");
  v1.push_back("3");
  v1.push_back("5");
  v1.push_back("7");

  ez::intersection(v1, v2, v3);
  ezassert( v3.empty() );

  v2.push_back("1");
  v2.push_back("2");
  v2.push_back("3");
  v2.push_back("4");
  
  ez::intersection(v1, v2, v3);
  ezassert( v3.size() == 2 );
  ezassert( v3[0].compare("1") == 0 );
  ezassert( v3[1].compare("3") == 0 );

  return true;
}

bool test_sort_nocase(ez::ezTestRunner& runner) {
  std::vector<std::string> v1;
  
  ez::sort_nocase(v1);
  ezassert( v1.empty() );

  v1.push_back("z");
  v1.push_back("g");
  v1.push_back("a");
  v1.push_back("f");
  
  ez::sort_nocase(v1);
  ezassert( v1.size() == 4 );
  ezassert( v1[0].compare("a") == 0 );
  ezassert( v1[1].compare("f") == 0 );
  ezassert( v1[2].compare("g") == 0 );
  ezassert( v1[3].compare("z") == 0 );
  
  v1.clear();
  v1.push_back("z");
  v1.push_back("G");
  v1.push_back("A");
  v1.push_back("f");

  ez::sort_nocase(v1);
  ezassert( v1.size() == 4 );
  ezassert( v1[0].compare("A") == 0 );
  ezassert( v1[1].compare("f") == 0 );
  ezassert( v1[2].compare("G") == 0 );
  ezassert( v1[3].compare("z") == 0 );

  v1.clear();
  v1.push_back("zz");
  v1.push_back("Ggg");
  v1.push_back("Aaa");
  v1.push_back("f");

  ez::sort_nocase(v1);
  ezassert( v1.size() == 4 );
  ezassert( v1[0].compare("Aaa") == 0 );
  ezassert( v1[1].compare("f") == 0 );
  ezassert( v1[2].compare("Ggg") == 0 );
  ezassert( v1[3].compare("zz") == 0 );

  return true;
}

bool test_split(ez::ezTestRunner& runner) {
  std::string s1;
  std::vector<std::string> v1;

  ez::split(s1, ',', v1);
  ezassert( v1.empty() );
  
  s1 = ",";
  v1.clear();
  ez::split(s1, ',', v1);
  ezassert( v1.empty() );

  s1 = ",,";
  v1.clear();
  ez::split(s1, ',', v1);
  ezassert( v1.empty() );
  
  s1 = " ";
  v1.clear();
  ez::split(s1, ',', v1);
  ezassert( v1.size() == 1 );
  ezassert( v1[0].compare(" ") == 0 );

  s1 = "asdf";
  v1.clear();
  ez::split(s1, ',', v1);
  ezassert( v1.size() == 1 );
  ezassert( v1[0].compare("asdf") == 0 );
  
  s1 = "a,b";
  v1.clear();
  ez::split(s1, ',', v1);
  ezassert( v1.size() == 2 );
  ezassert( v1[0].compare("a") == 0 );
  ezassert( v1[1].compare("b") == 0 );

  s1 = "1234,asdf,;'[]";
  v1.clear();
  ez::split(s1, ',', v1);
  ezassert( v1.size() == 3 );
  ezassert( v1[0].compare("1234") == 0 );
  ezassert( v1[1].compare("asdf") == 0 );
  ezassert( v1[2].compare(";'[]") == 0 );

  s1.clear();
  s1.reserve(49000);
  char tmp[8];
  for(int i=0; i < 10000; ++i) {
    sprintf(tmp,"%d,", i);
    s1.append(tmp);
  }
  v1.clear();
  ez::split(s1, ',', v1);
  ezassert( v1.size() == 10000 );
  ezassert( v1[0].compare("0") == 0 );
  ezassert( v1[9999].compare("9999") == 0 );
  
  return true;
}

bool test_vector_union(ez::ezTestRunner& runner) {
  std::vector<std::string> v1, v2, v3;
  
  ez::vector_union(v1, v2, v3);
  ezassert( v3.empty() );
  
  v1.push_back("1");
  v1.push_back("3");
  v1.push_back("5");
  v1.push_back("7");

  ez::vector_union(v1, v2, v3);
  ezassert( v3.size() == 4 );
  ezassert( v3[0].compare("1") == 0 );
  ezassert( v3[1].compare("3") == 0 );
  ezassert( v3[2].compare("5") == 0 );
  ezassert( v3[3].compare("7") == 0 );

  v2.push_back("1");
  v2.push_back("2");
  v2.push_back("3");
  v2.push_back("4");
  
  ez::vector_union(v1, v2, v3);
  ezassert( v3.size() == 6 );
  ezassert( v3[0].compare("1") == 0 );
  ezassert( v3[1].compare("2") == 0 );
  ezassert( v3[2].compare("3") == 0 );
  ezassert( v3[3].compare("4") == 0 );
  ezassert( v3[4].compare("5") == 0 );
  ezassert( v3[5].compare("7") == 0 );

  return true;
}

int main(int argc, const char* argv[]) {
  ez::ezTestRunner runner;
  runner.getopt(argc, argv);
  #define TEST(N) runner.add(&test_##N, "test_"#N);
  TEST(compare_nocase);
  TEST(find_first);
  TEST(lower);
  TEST(not_in_both);
  TEST(not_in_second);
  TEST(intersection);
  TEST(sort_nocase);
  TEST(split);
  TEST(vector_union);
  
  return runner.run();
}