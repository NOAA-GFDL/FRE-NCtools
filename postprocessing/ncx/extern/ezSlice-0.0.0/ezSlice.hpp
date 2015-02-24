/*
20111031 rsz Created.
*/
#ifndef EZNCSLICE_H
#define EZNCSLICE_H

#include <stdio.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <math.h> // ceil

namespace ez {

template <typename T>
class ezSlice {
// Parses and describes a dimension's index min,max,stride slicing for 1 axis.
// min,max,stride can be negatives. If min or max < 0, then they are 
// relative to end of a sequence. Negative stride denotes reverse order,
// but only makes sense if min > max.
// Assumes 0-based inclusive indexing where 0:3:1 refers to 4 elements.
// All indices returned by "get" functions are 0-based since this is C++ afterall. Fortran/Matlab 1-based indices can be "set" if offset=-1.
//
// Acceptable formats old comma style:
//  mn,mx,st
//  mn,mx[,]
//  ,mx,st
//  ,mx[,]
//  mn[,][,]
//  ,,st
//  mn,,st
//
// Accepted better colon style:
//  mx:mx:st
//  mx:mx[:]
//  :mx:st
//  :mx[:]
//  mn[:][:]
//  ::st
//  mx::st
//  ...
public:
  // Lowest index inclusive.
  T min;
  // Upper bound index (might not be including due to over-reaching stride).
  T max;
  // For fortran/matlab 1-based indices, use offset=-1 to make min and max 0-based. Does not affect stride or negative offsets since fortran/matlab don't know about them. Applied only within 'get' functions.
  int offset;
  // Index skip rate (default is 1).
  int stride;
  // If not supplied, and open-ended depending on real upper bounds.
  bool minMissing;
  // If not supplied, and open-ended depending on real upper bounds.
  bool maxMissing;
  // If not supplied.
  bool strideMissing;
  // If not a slice and just a simple index.
  bool singleton;
  
  // Resets.
  void clear() { min = 0; max = 0; stride = 1; offset = 0; minMissing = true; maxMissing = true; strideMissing = true; singleton = false; };
  // Ctor.
  ezSlice() { clear(); };
  // Get min,max,stride within context of length n.
  void get(T& _min, T& _max, T& _stride, T n);
  // Get list of indices for current min,max,stride.
  void get(std::vector<T>& indices, T n);
  // Print slice under context of size n for an array of indices 0 to n-1;
  void print(T n);
  // Set specific values.
  void set(T _min, T _max, T _stride, int offset=0);
  // Set with strings containing values.
  void set(const char* mn=0, const char* mx=0, const char* st=0, int offset=0);
  // Parse and set values found.
  void set(std::string& s, char delim=',', int offset=0);
  // Get number of elements in context of length n list.
  T size(T n);
};
//##################################################################
template <typename T>
void ezSlice<T>::get(T& _min, T& _max, T& _stride, T n) {
  ///printf("%s:%d\n", __FILE__, __LINE__);
  if (singleton) {
    _stride = 1;
    if (!minMissing) {
      _min = _max = min + offset;
    } else {
      _min = _max = 0;
    }
    return;
  }  
  ///printf("%s:%d\n", __FILE__, __LINE__);
  if (minMissing) {
    if (maxMissing) {
      if (strideMissing) {
        _min = 0;
        _max = n-1;
        _stride = 1;
      } else {
        _stride = stride;
        if (stride > 0) {
          _min = 0;
          _max = n-1;
        } else {
          _min = n-1;
          _max = 0;
        }
      }
    } else {
      // max not missing.
      if (strideMissing) {
        _min = 0;
        _max = max + offset;
        _stride = 1;        
      } else {
        _stride = stride;
        if (stride > 0) {
          _min = 0;
          _max = max + offset;
        } else {
          _min = n-1;
          _max = max + offset;
        }      
      }
    }
  } else { // min not missing.
    ///printf("%s:%d\n", __FILE__, __LINE__);
    if (min < 0) {
      _min = n+min;
    } else {
      ///printf("%s:%d\n", __FILE__, __LINE__);
      _min = min + offset;
    }
    ///printf("%s:%d\n", __FILE__, __LINE__);
    if (maxMissing) {
      _max = n-1;
    } else { // max is not missing.
      if (max < 0) {
        _max = n+max;
      } else {
      ///printf("%s:%d\n", __FILE__, __LINE__);
        _max = max + offset;      
      }
    }

    if (strideMissing) {
      _stride = 1;
    } else {
      _stride = stride;      
    }
  }
}
//##################################################################
template <typename T>
void ezSlice<T>::get(std::vector<T>& indices, T n) {
  indices.clear();

  if (singleton) {
    indices.push_back(min + offset);
    return;
  }

  T mn,mx,st;
  get(mn,mx,st,n);
  
  T i = mn;
  indices.reserve(size(n));
  
  if (st > 0) {
    for(; i <= mx; i += st)
      indices.push_back(i);    
  } else {
    for(; i >= mx; i += st)
      indices.push_back(i);
  }
}
//##################################################################
template <typename T>
void ezSlice<T>::print(T n) {
  T _min, _max, _stride;
  get(_min, _max, _stride, n);
  std::cout << _min << "," << _max << "," << _stride << "\n";
}
//##################################################################
template <typename T>
void ezSlice<T>::set(T _min, T _max, T _stride, int _offset) {
  min=_min;
  max=_max;
  stride=_stride;
  offset=_offset;
  minMissing=false;
  maxMissing=false;
  strideMissing=false;
}
//##################################################################
template <typename T>
void ezSlice<T>::set(const char* mn, const char* mx, const char* st, int _offset) {
  clear();
  offset = _offset;
  std::stringstream ss;

  if ((mn==0) || (mn[0]=='\0'))
    minMissing = true;
  else {
    ss.str(mn);
    ss >> min;
    minMissing = false;
  }
  
  if ((mx==0) || (mx[0]=='\0'))
    maxMissing = true;
  else {
    ss.clear();
    ss.str(mx);
    ss >> max;
    maxMissing = false;
  }
  
  if ((st==0) || (st[0]=='\0'))
    strideMissing = true;
  else {
    ss.clear();
    ss.str(st);
    ss >> stride;
    strideMissing = false;
  }
};
//##################################################################
template <typename T>
void ezSlice<T>::set(std::string& s, char delim, int _offset) {
  clear();
  offset = _offset;
  if (s.compare(0,3,"...")==0) {
    minMissing = 1;
    maxMissing = 1;
    strideMissing = 1;
    return;
  }
  
  int n = s.size();
  int i = 0;
  int state = 0;
  const char* c = s.c_str();
  char ch;
  int ndelim = 0;
  std::string tmp;
  tmp.reserve(8);
  std::stringstream ss;

  while (i < n) {
    ch = c[i++];
    if (ch == delim) { 
      ++ndelim; break; 
    } else
      tmp.append(1,ch);
  }
  if (tmp.size()) {
    ss.str(tmp);
    ss >> min;
    minMissing = 0;
  } else {
    minMissing = 1;
  }
  tmp.clear();

  while (i < n) {
    ch = c[i++];
    if (ch == delim) { 
      ++ndelim; break; 
    } else 
      tmp.append(1,ch);
  }

  if (tmp.size()) {
    ss.clear();
    ss.str(tmp);
    ss >> max;
    maxMissing = 0;
  } else {
    maxMissing = 1;
  }
  tmp.clear();

  while (i < n) {
    ch = c[i++];
    if (ch == delim) { 
      ++ndelim; break; 
    } else
      tmp.append(1,ch);
  }
  
  if (tmp.size()) {
    ss.clear();
    ss.str(tmp);
    ss >> stride;
    strideMissing = 0;
  } else {
    strideMissing = 1;
  }
    
  if (ndelim < 1)
    singleton = 1;
  
  //printf("%s:%d minMissing=%d, maxMissing=%d, singleton=%d, ndelim=%d\n", __FILE__, __LINE__, minMissing, maxMissing, singleton, ndelim);
}
//##################################################################
template <typename T>
T ezSlice<T>::size(T n) {
  T mn,mx,st;
  get(mn,mx,st,n);
  
  if (st > 0) {
    return ceil( (mx-mn+1.)/st );
  } else {
    if (mx > mn)
      return ceil( (mx-mn+1.)/st )*-1;
    else
      return floor( (mn-mx+1.)/st )*-1;
  }
}
//##################################################################
};
#endif // EZNCSLICE_H