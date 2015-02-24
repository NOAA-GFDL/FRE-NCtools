/*
Simple test framework. See ezTest.cpp for example.
On exit, returns number of test cases that failed, or 0 if all passed.

20111107 rsz Created.
20111116 rsz Added color.
*/

#ifndef EZTEST_H
#define EZTEST_H

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <string.h>
#include <sys/time.h> // gettimeofday

#define ezassert(x) \
if (x) {} \
else {\
  printf("%s:%u: failed assertion \"%s\"\n", __FILE__, __LINE__, #x); \
  fflush(stdout); \
  return false; \
}

#define EZTESTVERSION "0.0.0"

#define EZRED "\x1b[91m"
#define EZGREEN "\x1b[92m"
#define EZYELLOW "\x1b[93m"
#define EZBLUE "\x1b[94m"
#define EZEND "\x1b[0m"

namespace ez {
//########################################################################
// http://stackoverflow.com/questions/3283804/c-get-milliseconds-since-some-date
long long ezOsQueryPerformance() {
#ifdef WIN32
  LARGE_INTEGER llPerf = {0};
  QueryPerformanceCounter(&llPerf);
  return llPerf.QuadPart * 1000ll / ( g_llFrequency.QuadPart / 1000ll);
#else
  struct timeval stTimeVal;
  gettimeofday(&stTimeVal, NULL);
  return stTimeVal.tv_sec * 1000000ll + stTimeVal.tv_usec;
#endif
};
//########################################################################
class ezTestRunner;
// Function signature for all tests. Return a bool and input an ezTestRunner reference.
typedef bool (*PointerToTest)(ezTestRunner&);

class ezTestRunner {
public:
  // Ctor.
  ezTestRunner() : clean(true), color(false), verbose(false) {};
  // Add test function and name.
  void add(PointerToTest t, const char* n) { tests.push_back(t); names.push_back(n); };
  // Simple option parsing.
  void getopt(int argc, const char* argv[]);
  // Print version info.
  void printVersion();
  // Run the tests. Return 0 on all tests passing, 1 if any tests failed.
  int run();
  // Prints help message.
  void usage(const char* cmd);
  
  // Array of test numbers to run.
  std::vector<int> ids;
  // Test names/descriptions printed to console.
  std::vector<std::string> names;
  // No clean option flag, so tests know if temporary files should be erased or not.
  bool clean;
  // Enable colorized terminal messages.
  bool color;
  // Array of tests.
  std::vector<PointerToTest> tests;
  // Path for temporary files.
  std::string tmpdir;
  // Verbose flag.
  bool verbose;
};
//#####################################################################
void ezTestRunner::getopt(int argc, const char* argv[]) {
  if (argc < 2) return;
  
  // Look for any digits. They are test numbers to run. If none, run all.
  int i=1;
  ids.reserve(argc-1);
  for(; i < argc; ++i) {
    if(isdigit(argv[i][0])) {
      ids.push_back(atoi(argv[i]));
    } else if(strcmp(argv[i], "--clean")==0) {
      ++i;
      clean = atoi(argv[i]);
    } else if(strcmp(argv[i], "--color")==0) {
      color = true;
    } else if(strcmp(argv[i], "--tmpdir")==0) {
      ++i;
      tmpdir = argv[i];
    } else if(strcmp(argv[i], "--help")==0) {
      usage(argv[0]);
      exit(1);
    } else if(strcmp(argv[i], "--verbose")==0) {
      verbose = true;
    } else if(strcmp(argv[i], "--version")==0) {
      printVersion();
      exit(1);
    } else {
      fprintf(stderr, "ERROR: Unrecognized option \"%s\".\n", argv[i]);
    }
  }
}
//#####################################################################
void ezTestRunner::printVersion() {
  printf("ezTestRunner %s Copyright (C) 2011 Remik Ziemlinski\n", EZTESTVERSION);
  printf("This software is provided as-is for free without warranty.\n");
}
//#####################################################################
int ezTestRunner::run() {
  // Setup color constants;
  const char *R="", *G="", *B="", *Y="", *E="";
  if (color) {
    R = EZRED;
    G = EZGREEN;
    Y = EZYELLOW;
    B = EZBLUE;
    E = EZEND;
  }

  int i,n;
  if(ids.empty()) {
    n = tests.size();
    ids.resize(n);
    for(i=0; i < n; ++i)
      ids[i] = i+1;
  }
  
  n = ids.size();
  int nPass = 0;
  int testidx;
  int nValid = tests.size();
  long long tic, toc, totalTic, totalToc;
  bool pass;
  totalTic = ezOsQueryPerformance();
  
  for(i=0; i < n; ++i) {
    testidx = ids[i]-1;
    if ( (testidx >= nValid) || (testidx < 0) ) {
      printf("%s%d. Invalid test.\nFAIL%s\n", R, testidx+1, E);
      fflush(stdout);
    } else {
      printf("%s%d. %s%s\n", B, ids[i], names[testidx].c_str(), E);
      fflush(stdout);
  
      tic = ezOsQueryPerformance();
      pass = (*tests[testidx])(*this);
      toc = ezOsQueryPerformance();
      
      if (pass) {
        if (verbose)
          printf("%sPASS", G);
        ++nPass;
      } else {
        printf("%sFAIL", R);
      }
      if (verbose || !pass)
        printf(" (%.5f sec)%s\n", (toc-tic)/1.0e6, E);
    }
    if (verbose)
      printf("%s-----------------------------------------------------------------%s\n", B, E);
  }
  totalToc = ezOsQueryPerformance();
  
  const char* finalcolor = Y;
  if (nPass == n) 
    finalcolor = G;
  else if (nPass == 0)
    finalcolor = R;
    
  printf("%s%d of %d tests passed in %.5f seconds.%s\n", finalcolor, nPass, n, (totalToc-totalTic)/1.e6, E);
  
  return n - nPass;
}
//#####################################################################
void ezTestRunner::usage(const char* cmd) {
  printf("%s -- Simple unit test framework runner.\n\n", cmd);
  printf("USAGE: %s [OPTIONS] [test numbers]\n\n", cmd);
  printf("OPTIONS:\n\n");
  printf("  --help         Prints this message.\n");
  printf("  --color        Colorize terminal messages.\n");
  printf("  --clean C      Temporary files should be erased? 0 or 1\n");
  printf("  --tmpdir S     Use S as the path for temporary files.\n");
  printf("  --verbose      Hint for verbose output.\n");
  printf("  --version      Print version and exit.\n");
  printf("\n");
  printf("EXAMPLES:\n\n");
  printf("  %s --tmpdir /tmp/ --noclean 1 3 5 10\n", cmd);
  printf("\n");
  printVersion();
}
//#####################################################################
};
#endif // EZTEST_H