/*
This file is part of ezOdometer.

ezOdometer is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ezOdometer is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with ezOdometer.  If not, see <http://www.gnu.org/licenses/>.

Copyright (C) 2011 Remik Ziemlinski
*/
/*
CHANGELOG
v0.0.0 20110514 rsz Created.
v0.1.0 20111115 rsz Refactored ++ operators, added isDone and inc* functions.
       20111123 rsz Added normalize functions.
v0.1.1 20111124 rsz Redefined isDone to be more like STL iterators.
*/
#ifndef EZODOMETER_H_
#define EZODOMETER_H_

#include <string.h>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <math.h> // ceil

namespace ez {

template<typename T>
class ezOdometer {
// Mimics custom incrementing odometer with rollover detection.
public:
  // These are the min/max number range per digit with stride (defaults are 1).
  std::vector<T> mins, maxs, strides;
  // The largest possible valid values.
  std::vector<T> limits;
  // These are the digits that will update.
  std::vector<T> counters;
  // Only digits that have non-zero mask value will be updated.
  std::vector<char> locks;
  // Index of uppermost digit. 
  unsigned char begin;
  // Index of last/lowermost digit. In left-to-right sequence, this is the rightmost digit index.
  unsigned char end;
  // Flags for truly done incrementing. Only set by ++ operator. Reset by clear and reset.
  std::vector<char> done;
  
public:
  // Allocates arrays. User must call free before calling.
  ezOdometer();
  ezOdometer(unsigned char n);
  ~ezOdometer();
  // Returns true if digit at position "i" hit its limit value. 
  bool atLimit(unsigned char i);
  // Returns true if all digits within begin/end hit their limit values. 
  bool atLimits();
  // De-allocates arrays.
  void clear();
  // Get current digits.
  void get(std::vector<T> & digits) { digits = counters; }
  // Returns number of unique digit combinations if completely iterated.
  unsigned int getNumCombinations();
  // Reserves space for arrays.
  void init(unsigned char n);
  // Returns union of done flags within begin/end.
  bool isDone();
  // If set, then digit will be fixed to a value.
  void lock(int i, T value);
  // Sets mins to zero, maxes to size minus one, and strides to one. 
  void normalize();
  // Access a digit by index.
  T& operator[] (const unsigned int index);
  // Increments.
  void operator++();
  void print(char delim=0);
  // Resets counters to all zeros so can be incremented all over again.
  void reset();
  // Set a digits range. Default is 0 to 9. Can be whacky like -10 to 123.
  void setrange(int index, T min=0, T max=9, T stride=1);
  // Returns number of digits.
  int size() { return counters.size(); }
};
//#################################################################
template<typename T>
ezOdometer<T>::ezOdometer() {
  // Make begin after end, so loops won't execute until "init" is called.
  begin = 1;
  end = 0;
}
//#################################################################
template<typename T>
ezOdometer<T>::ezOdometer(unsigned char n) {
  init(n);
};
//#################################################################
template<typename T>
ezOdometer<T>::~ezOdometer() {
  clear();
};
//#################################################################
template<typename T>
bool ezOdometer<T>::atLimit(unsigned char i) {
  if (i < counters.size())
    return counters[i] == limits[i];
  
  return true;
}
//#################################################################
template<typename T>
bool ezOdometer<T>::atLimits() {
  for(char i = begin; i <= end; ++i) {
    if(!locks[i]) {
      ///printf("counters[%d]=%d, limits[%d]=%d\n", i, counters[i], i, limits[i]);
      if (counters[i] != limits[i])
        return false;
    }
  }
  
  return true;
}
//#################################################################
template<typename T>
unsigned int ezOdometer<T>::getNumCombinations() {
  if (mins.empty()) return 0;
  
  unsigned int n = 1;
  int i = 0;
  int m = mins.size();
  for(; i < m; ++i) {
    if (locks[i] || (maxs[i] == mins[i]))
      continue;
      
    n *= ((maxs[i]-mins[i])/strides[i] + 1);
  }

  return n;
}
//#################################################################
template<typename T>
void ezOdometer<T>::init(unsigned char n) {
  counters.resize(n);
  done.resize(n);
  limits.resize(n);
  mins.resize(n);
  maxs.resize(n);
  strides.resize(n);
  locks.resize(n);
  if (n!=0) {
    begin = 0;
    end = n-1;
  } else {
    // Ensure no loops execute on zero length arrays.
    begin = 1;
    end = 0;
  }
  
  int i;
  for(i=0; i < n; ++i) {
    counters[i]=0;
    done[i]=0;
    limits[i]=0;
    mins[i]=0;
    maxs[i]=0;
    strides[i]=1;
    locks[i]=0;
  }
}
//#################################################################
template<typename T>
bool ezOdometer<T>::isDone() {
  for(char i = begin; i <= end; ++i)
    if (!done[i]) return false;
    
  return true;
}
//#################################################################
template<typename T>
void ezOdometer<T>::clear() {
  done.clear();
  limits.clear();
  mins.clear();
  maxs.clear();
  counters.clear();
  locks.clear();
}
//#################################################################
template<typename T>
void ezOdometer<T>::lock(int i, T value) {
  counters[i] = value;
  limits[i] = value;
  locks[i] = 1;
  maxs[i] = value;
  mins[i] = value;
}
//#################################################################
template<typename T>
void ezOdometer<T>::normalize() {
  int i=0;
  int n = counters.size();
  int m;
  for(; i < n; ++i) {
    m = (maxs[i]-mins[i])/strides[i];
    mins[i] = 0;
    counters[i] = 0;
    limits[i] = m;
    maxs[i] = m;
    strides[i] = 1;
  }
}
//#################################################################
template<typename T>
T& ezOdometer<T>::operator[](const unsigned int index) {
  return counters[index];
}
//#################################################################
template<typename T>
void ezOdometer<T>::operator++() {
  if (isDone()) return;
  
  if (atLimits()) {
    // All must be at limits in index range to be set as done.
    // Like STL iterator, increment past last value goes over edge.
    for(char j=begin; j <= end; ++j) done[j] = 1;
    return;
  }
  
  // Can't use unsigned for counter since when need to loop below zero.
  char i = end;
  while (i >= begin) {
    if (!locks[i]) {
      if (counters[i] == limits[i])
        counters[i] = mins[i];
      else {
        counters[i] += strides[i];
        break;
      }
    }
    --i;
  }
};
//#################################################################
template<typename T>
void ezOdometer<T>::print(char delim) {
  int i;
  if (delim) {
    for(i=0; i < counters.size()-1; ++i)
      std::cout << counters[i] << delim;
      
      std::cout << counters[i];
  } else {
    for(i=0; i < counters.size(); ++i)
      std::cout << counters[i];
  }
};
//#################################################################
template<typename T>
void ezOdometer<T>::reset() {
  for(unsigned char i=begin; i <= end; ++i) {
    if (!locks[i]) {
      counters[i] = mins[i];
    }
    done[i] = 0;
  }
};
//#################################################################
template<typename T>
void ezOdometer<T>::setrange(int i, T min, T max, T stride) {
  counters[i] = min;
  done[i] = 0;
  locks[i] = 0;
  mins[i] = min;
  maxs[i] = max;
  strides[i] = stride;
  
  int n;
  if (stride > 0) {
    n = ceil( (max-min+1.)/stride );
  } else {
    if (max > min)
      n = ceil( (max-min+1.)/stride )*-1;
    else
      n = floor( (min-max+1.)/stride )*-1;
  }
  ///printf("min=%d, max=%d, stride=%d, n=%d\n", (int)min, (int)max, (int)stride);
  limits[i] = min + (n-1)*stride;
};
//#################################################################
// Allows custom values per digit for non-linear odometer, for example
// a digit can be only even numbers, or primes, or arbitrary order.
// Also, the digit isn't limited to single digit as classic 0-9, instead
// each can be arbitrary, which in that case this odometer becomes 
// just an array combination generator.
// During increment, if vector of custom values for a digit is empty, then
// the counter itself will be the digit.
template<typename T>
class ezCustomOdometer : public ezOdometer<T> {
public:
  // List of custom values per digit.
  std::vector< std::vector<T> > values;
  // Final result you want. The custom mapped digits for the current state of the odometer. 
  std::vector<T> digits;
  
  // Ctor. 
  ezCustomOdometer() { };
  // Construct with "n" digits wide.
  ezCustomOdometer(unsigned char n);
  // Dtor.
  ~ezCustomOdometer();
  // If digit at index 'i' is at largest value.
  bool atLimit(unsigned char i);
  // If all unlocked digits within begin/end range are at largest values.
  bool atLimits();
  // Clears arrays.
  void clear();
  // Returns number of unique digit combinations if completely iterated.  
  unsigned int getNumCombinations();
  // Allocates memory. Be sure to call "clear" before using this.
  void init(unsigned char n);
  // Sets mins to zero, maxes to size minus one, and strides to one. 
  void normalize();
  // Access a digit.
  T& operator[](const unsigned int index);
  // Increments.
  void operator++();
  void print(char delim=0);
  // Sets all digits to initial values to reallow incrementing from beginning.
  void reset();
  // Sets an index to have only one locked value and sets its digit.
  void setCustomValue(int index, T v);
  // Set the vector of custom values for a digit, and reset digit to first value.
  void setCustomValues(int index, std::vector<T> & v);
  // Allow strided range.
  void setrange(int i, T min, T max, T stride);
};
//#################################################################
template<typename T>
ezCustomOdometer<T>::ezCustomOdometer(unsigned char n) {
  init(n);
};
//#################################################################
template<typename T>
ezCustomOdometer<T>::~ezCustomOdometer() {
  clear();
};
//#################################################################
template<typename T>
bool ezCustomOdometer<T>::atLimit(unsigned char i) {
  if (i < digits.size())
    return digits[i] == this->limits[i];
  
  return true;
}
//#################################################################
template<typename T>
bool ezCustomOdometer<T>::atLimits() {
///  printf("%s:%d begin=%d, end=%d\n", __FILE__, __LINE__, (int)this->begin, (int)this->end);
  for(char i = this->begin; i <= this->end; ++i) {
///    printf("%s:%d i=%d\n", __FILE__, __LINE__, (int)i);
    if(!this->locks[i]) {
///      printf("%s:%d digits[%d]=%d, limits[%d]=%d\n", __FILE__, __LINE__, i, (int)this->digits[i], i, (int)this->limits[i]);
      if (digits[i] != this->limits[i]) {
///        printf("%s:%d\n", __FILE__, __LINE__);
        return false;
      }
    }
  }
  
  return true;
}
//#################################################################
template<typename T>
void ezCustomOdometer<T>::clear() {
  digits.clear();
  values.clear();
  ezOdometer<T>::clear();
}
//#################################################################
template<typename T>
unsigned int ezCustomOdometer<T>::getNumCombinations() {
  if (digits.empty()) return 0;
  
  unsigned int n = 1;
  int i = 0;
  int m = digits.size();
  for(; i < m; ++i) {
    if (values[i].empty())
      n *= ((this->maxs[i]-this->mins[i])/this->strides[i] + 1);
    else
      n *= values[i].size();
  }
  
  return n;
}
//#################################################################
template<typename T>
void ezCustomOdometer<T>::init(unsigned char n) {
  ezOdometer<T>::init(n);
  digits.resize(n);
  values.resize(n);  
}
//#################################################################
template<typename T>
void ezCustomOdometer<T>::normalize() {
  ezOdometer<T>::normalize();
  int i=0;
  int n = this->counters.size();
  int m;
  for(; i < n; ++i) {
    digits[i] = 0;
    values[i].clear();
  }
}
//#################################################################
template<typename T>
T& ezCustomOdometer<T>::operator[](const unsigned int index) {
  return digits[index];
}
//#################################################################
template<typename T>
void ezCustomOdometer<T>::operator++() {
///printf("%s:%d\n", __FILE__, __LINE__);
  if (this->isDone()) return;
///printf("%s:%d\n", __FILE__, __LINE__);
  
  if (atLimits()) {
    // All must be at limits in index range to be set as done.
    // Like STL iterator, increment past last value goes over edge.
    for(char j=this->begin; j <= this->end; ++j) this->done[j] = 1;
    return;
  }

  // Can't use unsigned for counter since when need to loop below zero.
  char i = this->end;
  while (i >= this->begin) {
    if (!this->locks[i]) {
      if (this->digits[i] == this->limits[i]) {
        if (values[i].empty()) {
          digits[i] = this->mins[i];
        } else {
          digits[i] = values[i][0];
        }
        this->counters[i] = this->mins[i];
      } else {
        this->counters[i] += this->strides[i];
        if (values[i].empty())
          digits[i] = this->counters[i];
        else
          digits[i] = values[i][this->counters[i]];
        
        break;
      }
    }
    --i;
  }
};
//#################################################################
template<typename T>
void ezCustomOdometer<T>::print(char delim) {
  int i=0;
  int n = digits.size();
  if (delim) {
    #define PRINTARR(A) {\
      /*std::cout << #A": ";*/\
      for(; i < n-1; ++i)\
        std::cout << (int)A[i] << delim;\
      std::cout << (int)A[i];\
    }
    PRINTARR(this->digits);
    /*PRINTARR(this->counters);
    PRINTARR(this->limits);
    PRINTARR(this->locks);
    PRINTARR(this->mins);
    PRINTARR(this->maxs);
    PRINTARR(this->strides);
    */
  } else {
    for(; i < digits.size(); ++i)
      std::cout << digits[i];
  }
};
//#################################################################
template<typename T>
void ezCustomOdometer<T>::reset() {
  ezOdometer<T>::reset();

  for(unsigned char i=this->begin; i <= this->end; ++i) {
    if (!this->locks[i]) {
      if (values[i].empty())
        digits[i] = this->mins[i];
      else
        digits[i] = values[i][0];
    }
  }
};
//#################################################################
template<typename T>
void ezCustomOdometer<T>::setCustomValue(int i, T v) {
  values[i].clear();
  this->mins[i] = v;
  digits[i] = v;
  ezOdometer<T>::lock(i,v);
};
//#################################################################
template<typename T>
void ezCustomOdometer<T>::setCustomValues(int i, std::vector<T> & v) {
  this->counters[i] = 0;
  this->mins[i] = 0;
  this->maxs[i] = v.size()-1;
  values[i] = v; // copy vector.
  this->strides[i] = 1; // ensure simple incrementing of counter.
  if (v.size()) {
    digits[i] = v[0];
    this->limits[i] = v[v.size()-1];
  } else {
    digits[i] = 0;
    this->limits[i] = 0;
  }
};
//#################################################################
template<typename T>
void ezCustomOdometer<T>::setrange(int i, T min, T max, T stride) {
  ezOdometer<T>::setrange(i,min,max,stride);
  digits[i] = min;
  values[i].clear();
}
}
#endif // EZEODOMETER_H_