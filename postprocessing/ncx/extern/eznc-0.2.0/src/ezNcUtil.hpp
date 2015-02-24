/*
*/
/*
CHANGELOG
v0.0.0 20110514 rsz Created.
*/
#ifndef EZNCUTIL_H_
#define EZNCUTIL_H_

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <ctime>
#include <sys/time.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <limits>
#include <regex.h>
#include "netcdf.h"
#include <math.h>

namespace ez {
//###########################################################################
// Get path prefix to a filename. Has trailing slash. Linux or Windows.
void DirName(std::string& path, std::string& prefix) {
  if (path.empty()) return;
  
  size_t i = path.find_last_of("/\\");
  if (i == std::string::npos) return;
  
  prefix.assign( path.substr(0, i+1) );
}
//###########################################################################
// Test if file exists given path and name.
bool FileExists( const char* FileName ) {
  FILE* fp = fopen( FileName, "rb" );
  if( fp != NULL ) {
    fclose( fp );
    return true;
  }

  return false;
};
//###########################################################################
// Converts format to a human readable version.
void formatToString(int f, std::string & s) {
  s.clear();
  
  switch(f) {
    case NC_FORMAT_CLASSIC: s = "NC_FORMAT_CLASSIC"; break;
    case NC_FORMAT_64BIT: s = "NC_FORMAT_64BIT"; break;
    case NC_FORMAT_NETCDF4: s = "NC_FORMAT_NETCDF4"; break;
    case NC_FORMAT_NETCDF4_CLASSIC: s = "NC_FORMAT_NETCDF4_CLASSIC"; break;
    default: s = "UNKNOWN"; break;
  }
};
//###########################################################################
int HandleNcStatus(int status) {
  if (status != NC_NOERR) {
    fprintf(stderr, "%s\n", nc_strerror(status));
  }
  return status;
};
//###########################################################################
// Check if list of numbers is a whole number increasing or decreasing sequence.
template <typename T>
bool isRevSeq(std::vector<T>& list) {
  int n = list.size();
  if (n == 1) {
    return true;
  }
  
  n -= 1;
  int i = 0;
  T delta;
  for(; i < n; ++i) {
    // If abs(stride), aka delta, between 2 elements is > +/-1, then false.
    delta = list[i+1]-list[i];
    if( (delta != -1) ) return false;
  }
  
  return true;
}
//###########################################################################
// Check if list of numbers is a whole number increasing or decreasing sequence.
template <typename T>
bool isSeq(std::vector<T>& list) {
  int n = list.size();
  if (n == 1) {
    return true;
  }
  
  n -= 1;
  int i = 0;
  T delta;
  for(; i < n; ++i) {
    // If abs(stride), aka delta, between 2 elements is > +/-1, then false.
    delta = list[i+1]-list[i];
    if( (delta != 1) ) return false;
  }
  
  return true;
}
//###########################################################################
// Check if list of numbers is a strided increasing or decreasing sequence.
template <typename T>
bool isStrideSeq(std::vector<T>& list) {
  int n = list.size();
  if (n <= 2) {
    return true;
  }
  
  n -= 1;
  T prevdelta = list[1]-list[0];
  T delta;
  int i = 1;
  for(; i < n; ++i) {
    delta = list[i+1]-list[i];
    if (delta != prevdelta) return false;
    prevdelta = delta;
  }
  
  return true;
}
//###########################################################################
// process_mem_usage(double &, double &) - takes two doubles by reference,
// attempts to read the system-dependent data for a process' virtual memory
// size and resident set size, and return the results in KB.
//
// On failure, returns 0.0, 0.0
void process_mem_usage(double& vm_usage, double& resident_set) {
  using std::ios_base;
  using std::ifstream;
  using std::string;

  vm_usage     = 0.0;
  resident_set = 0.0;

  // 'file' stat seems to give the most reliable results
  //
  ifstream stat_stream("/proc/self/stat",ios_base::in);

  // dummy vars for leading entries in stat that we don't care about
  //
  string pid, comm, state, ppid, pgrp, session, tty_nr;
  string tpgid, flags, minflt, cminflt, majflt, cmajflt;
  string utime, stime, cutime, cstime, priority, nice;
  string O, itrealvalue, starttime;

  // the two fields we want
  //
  unsigned long vsize;
  long rss;

  stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
             >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
             >> utime >> stime >> cutime >> cstime >> priority >> nice
             >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

  stat_stream.close();

  long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
  vm_usage     = vsize / 1024.0;
  resident_set = rss * page_size_kb;
};
//###########################################################################
// http://stackoverflow.com/questions/3283804/c-get-milliseconds-since-some-date
long long osQueryPerformance() {
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
//###########################################################################
// Given a regex pattern, choose each string that has a match and append to output vector.
// Returns 0 if no fatal issues.
int regex(std::string& pattern, std::vector<std::string>& strings, std::vector<std::string>& matches, bool extended=false, bool ignorecase=false) {
  regex_t regex_struct;
  int regflags = (extended ? REG_EXTENDED : 0) | (ignorecase ? REG_ICASE : 0);
  
  int status = regcomp(&regex_struct, pattern.c_str(), regflags);
  if (status) {
    fprintf(stderr, "ERROR: Could not compile regular expression \"%s\".\n", pattern.c_str());
    regfree(&regex_struct);
    return 1;
  }
  
  int n = strings.size();
  int i = 0;
  char errmsg[128];
  for(; i < n; ++i) {
    status = regexec(&regex_struct, strings[i].c_str(), 0, 0, 0);
    switch(status) {
      case 0:
        matches.push_back(strings[i]);
        break;
      case REG_NOMATCH:
        break;
      default:
        regerror(status, &regex_struct, errmsg, sizeof(errmsg));
        fprintf(stderr, "ERROR: Regular expression failed: %s\n", errmsg);
        regfree(&regex_struct);
        return 1;
    }
  }
  
  regfree(&regex_struct);
  return 0;
}
//###########################################################################
double round(double r) {
  return (r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5);
}
//###########################################################################
// Return simple sequence from zero to n-1. 'list' must be pre-allocated.
template <typename T>
void seq(T* list, int n) {
  if (list == 0) return;
  
  for(int i=0; i < n; ++i)
    list[i] = i;
};
//###########################################################################
// Return simple sequence from zero to n-1.
template <typename T>
void seq(std::vector<T>& list, int n) {
  list.clear();
  list.resize(n);
  for(int i=0; i < n; ++i)
    list[i] = i;
};
//###########################################################################
// Return simple sequence from min to max inclusive.  min can be greater than max for reversed sequence.
template <typename T>
void seq(T* list, T mn, T mx) {
  if (list == 0) return;
  int n,i;
  
  if (mn < mx) {
    n = mx-mn+1;
    for(i=0; i < n; ++i)
      list[i] = mn+i;
  } else {
    n = mn-mx+1;
    for(i=0; i < n; ++i)
      list[i] = mn-i;  
  }
};
//###########################################################################
// Return simple sequence from min to max inclusive.  min can be greater than max for reversed sequence.
template <typename T>
void seq(std::vector<T>& list, T mn, T mx) {
  list.clear();
  int n,i;
  
  if (mn < mx) {
    n = mx-mn+1;
    list.resize(n);
    for(i=0; i < n; ++i)
      list[i] = mn+i;
  } else {
    n = mn-mx+1;
    list.resize(n);
    for(i=0; i < n; ++i)
      list[i] = mn-i;  
  }
};
//###########################################################################
int sizeOf(nc_type type) {
  switch(type) {
    case NC_CHAR: return sizeof(char); 
    case NC_BYTE: case NC_UBYTE: return sizeof(unsigned char);
    case NC_SHORT: return sizeof(short);
    case NC_USHORT: return sizeof(unsigned short);
    case NC_INT: return sizeof(int);
    case NC_UINT: return sizeof(unsigned int);
    case NC_INT64: return sizeof(long long);
    case NC_UINT64: return sizeof(unsigned long long);
    case NC_FLOAT: return sizeof(float);
    case NC_DOUBLE: return sizeof(double);
    default: return 0;
  }
};
//###########################################################################
/* 
Convert a vector with "min[,max[,stride]]" to list of numbers.

offset: If using fortran style, then set offset=-1.
max: Maximum allowable number in array. Helps define limit when generating sequence.

Examples:
str=0, offset=0, max=3, v=0,1,2,3
str=1, offset=-1, max=3, v=0,1,2,3
str=0,3, offset=0, max=5, v=0,1,2,3
str=1,3, offset=-1, max=5, v=0,1,2
str=1,10,2, offset=0, max=15, v=1,3,5,7,9
str=1,10,2, offset=-1, max=15, v=0,2,4,6,8 
*/
void MinMaxStrideToArray(std::vector<size_t> & mms, std::vector<size_t> & out, int offset=0, int max=0) {
  int mini=0;
  int maxi=max;
  int stride=1;
  
  switch(mms.size()) {
    case 3:
      stride = mms[2];
    case 2:
      maxi = mms[1]+offset;
      maxi = (max < maxi) ? max : maxi;
    case 1:
      mini = mms[0]+offset;
      break;
    default:
      return;
  }
   
  out.clear();
  for(int i=mini; i <= maxi; i += stride) {
    out.push_back(i);
  }
};
//###########################################################################
// Parses comma-delimited string.
// mms: comma-delimited string, such as 0,9,3.
// out: the substrings parsed and converted to list of numerical values.
// offset: scalar added to parsed values. fortran indexing would use -1.
// max: optional ceiling value if mms has less than 2 values.
void MinMaxStrideToArray(std::string & mms, std::vector<size_t> & out, int offset=0, int max=0) {
  std::stringstream ss(mms);
  std::vector<size_t> vec;
  size_t i;
  while (ss >> i) {
    vec.push_back(i);
    if (ss.peek() == ',') ss.ignore();
  }
  
  MinMaxStrideToArray(vec, out, offset, max);
};
//###########################################################################
// Parses comma-delimited string.
void MinMaxStrideToArray(const char* mms, std::vector<size_t> & out, int offset=0, int max=0) {
  std::stringstream ss(mms);
  std::vector<size_t> vec;
  size_t i;
  while (ss >> i) {
    vec.push_back(i);
    if (ss.peek() == ',') ss.ignore();
  }
  
  MinMaxStrideToArray(vec, out, offset, max);
};
//###########################################################################
// Expected input strings: 
// byte,char,short,long,float,double,string,int,int64,ubyte,uchar,ushort,uint,uint64
nc_type StringToNcType(const char* str) {
  if (str==0) return NC_NAT;
  
  int n = strlen(str);
  if (n < 2) return NC_NAT;
  switch(str[0]){
    case 'u':
      switch(str[1]) {
        case 'b': case 'c': return NC_UBYTE;
        case 's': return NC_USHORT;
        case 'i':
          switch(n) {
            case 4: return NC_UINT;
            default:return NC_UINT64;
          }
        default: return NC_NAT;
      }
    case 'b': return NC_BYTE;
    case 'c': return NC_CHAR;
    case 's': 
      switch(str[1]) {
        case 'h': return NC_SHORT;
        case 't': return NC_STRING;
        default: return NC_NAT;
      }
    case 'l': return NC_INT;
    case 'f': return NC_FLOAT;
    case 'd': return NC_DOUBLE;
    case 'i':
      switch(n) {
        case 3: return NC_INT;
        default:return NC_INT64;
      }
    default: return NC_NAT;
  } 
}
//###########################################################################
// Converts strings to preallocated array. Better to use vector version below.
// Example1: 
//   inputs: <"12","-1","0">, 1, 2, int[2]
//   output: int[2] = [-1,0]
template <typename T>
void StringsToValues(std::vector<std::string> & vec, int first, int last, T* arr) {
  if (arr == 0) return;
  int ctr=0;
  for(; first <= last; ++first) {
    std::stringstream ss(vec[first]);
    ss >> arr[ctr++];
  }
};
//###########################################################################
// Converts strings to vector of values.
template <typename T>
void StringsToValues(std::vector<std::string> & vec, int first, int last, std::vector<T>& arr) {
  arr.clear();
  arr.resize(last-first+1);
  StringsToValues(vec, first, last, &(arr[0]));
};
//###########################################################################
// Appends to output string. so allows user to prefix with a path.
void TmpFilename(std::string& out, const char* suffix) {
  char *tmp = new char[out.size() + 12];
  sprintf(tmp, "%s%stmp.XXXXXX", out.c_str(), (out.size() ? "/" : ""));
  
  int fd = -1;
  fd = mkstemp(tmp);
  if (fd != -1) {
    // mkstemp creates an empty file, so if works, we must delete it.
    remove(tmp);
    out = tmp;
  } else {
    out.clear();
  }
  
  if (out.size() && suffix)
    out.append(suffix);
    
  delete [] tmp;
};
//###########################################################################
}
#endif // EZNCUTIL_H_