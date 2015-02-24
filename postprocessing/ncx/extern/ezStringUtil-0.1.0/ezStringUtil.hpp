/*
This file is part of ezStringUtil.

ezStringUtil is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ezStringUtil is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with ezStringUtil.  If not, see <http://www.gnu.org/licenses/>.

Copyright (C) 2011 Remik Ziemlinski
*/
/*
CHANGELOG

v0.0.0 20110502 rsz Created.
v0.0.0 20111119 rsz Updated, packaged, unit and memory tested.
v0.1.0 20111202 rsz Added overloaded functions.
*/

#ifndef EZ_STRINGUTIL_H
#define EZ_STRINGUTIL_H

#include <vector>
#include <string>
#include <set>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <ctype.h>
#include <stdio.h>

namespace ez {
//#####################################################################
// True if first string is less than other (lowered) only if chars are less or equal up to its shorter length than second string.
bool compare_nocase (const std::string &first, const std::string &second)
{
	unsigned int i=0;
	unsigned int n1 = first.length();
	unsigned int n2 = second.length();
	unsigned int n = ( (n1 < n2) ? n1 : n2 );
  char c1,c2;
	while ( i < n )	{
    c1 = tolower(first[i]);
    c2 = tolower(second[i]);
		if (c1 < c2) return true;
		else if (c1 > c2) return false;
		++i;
	}
	if (n1 < n2) return true;
	else return false;
};
//#####################################################################
// Make string all lower case.
void lower(std::string& string) {
  int i = 0;
  int n = string.size();
  char c;
  const char* s = string.c_str();
  for(; i < n; ++i) {
    c = s[i];
    if (c<='Z' && c>='A')
      string[i] = c+32;
  }
}
//#####################################################################
// Return index of exact first match of string in unsorted vector. -1 if no match.
int find_first(std::vector<std::string>& v, std::string& s, bool ignorecase=false) {
  if (v.empty() || s.empty()) return -1;
  int res = -1;
  int i = 0;
  int n = v.size();
  
  if (ignorecase) {
    std::string t;
    for(; i < n; ++i) {
      t = v[i];
      lower(t);
      if (t.compare(s) == 0) {
        res = i;
        break;
      }
    }
    
  } else {
    for(; i < n; ++i) {
      if (v[i].compare(s) == 0) {
        res = i;
        break;
      }
    }
  }
  
  return res;
}
//#####################################################################
int find_first(std::vector<std::string>& v, const char* s, bool ignorecase=false) {
  if (v.empty() || (s==0)) return -1;
  int res = -1;
  int i = 0;
  int n = v.size();
  
  if (ignorecase) {
    std::string t;
    for(; i < n; ++i) {
      t = v[i];
      lower(t);
      if (t.compare(s) == 0) {
        res = i;
        break;
      }
    }
    
  } else {
    for(; i < n; ++i) {
      if (v[i].compare(s) == 0) {
        res = i;
        break;
      }
    }
  }
  
  return res;
}
//#####################################################################
// Return elements that are present in one of the sets, but not in the other.
// For example, set1=1,2,3,4,5,6, set2=2,4,5,7 then result=1,3,6,7.
// Inputs are by copy because they must be sorted.
void not_in_both(std::vector<std::string> v1, std::vector<std::string> v2, std::vector<std::string> & difference) {
	std::sort(v1.begin(),v1.end());
	std::sort(v2.begin(),v2.end());

  difference.clear();
	std::set_symmetric_difference(v1.begin(), v1.end(), v2.begin(), v2.end(), std::inserter(difference, difference.end()) );
};
//#####################################################################
// Return elements that are present in the first set, but not in the second.
// For example, set1=1,2,3,4,5,6, set2=2,4,5, then result=1,3,6.
void not_in_second(std::vector<std::string> & v1, std::vector<std::string> & v2, std::vector<std::string> & difference) { 
	std::sort(v1.begin(),v1.end());
	std::sort(v2.begin(),v2.end());
  difference.clear();
	std::set_difference(v1.begin(), v1.end(), v2.begin(), v2.end(), std::inserter(difference, difference.end()));
};
//#####################################################################
// For example, set1=1,2,3,4,5,6, set2=2,4,5, then result=2,4,5.
void intersection(std::vector<std::string> & set1, std::vector<std::string> & set2, std::vector<std::string> & inter) { 
	std::sort(set1.begin(),set1.end());
	std::sort(set2.begin(),set2.end());
  inter.clear();
	std::set_intersection(set1.begin(), set1.end(), set2.begin(), set2.end(), std::inserter(inter, inter.end()));
};
//#####################################################################
void print(std::vector<std::string> & v, char delim) { 
	int n = v.size()-1;
	for(int i=0; i < n; ++i)
		std::cout << v[i] << delim;
		
	if (n>=0)
		std::cout << v[n];
};
//#####################################################################
void sort_nocase(std::vector<std::string> & v) { 
	std::sort(v.begin(), v.end(), ez::compare_nocase); 
};
//#####################################################################
// For example, in="1,2,3," then result=<"1","2","3">.
void split( const std::string& s, const char token, std::vector<std::string> & result) {
	std::string::const_iterator i = s.begin();
	std::string::const_iterator j = s.begin();
	const std::string::const_iterator e = s.end();

	while(i!=e) {
		while(i!=e && *i++!=token);
    result.push_back(std::string());
    std::string & newstr = result.back();
    newstr.assign(j, (*(i-1) == token) ? i-1 : i);
    if (newstr.empty()) result.pop_back();
		j = i;
	}
};
//#####################################################################
// For example, v1=1,2,3,4,5,6, v2=0,2,4,5,7 then result=0,1,2,3,4,5,6,7.
void vector_union(std::vector<std::string> & v1, std::vector<std::string> & v2, std::vector<std::string> & result) { 
	std::sort(v1.begin(),v1.end());
	std::sort(v2.begin(),v2.end());

  result.clear();
	std::set_union(v1.begin(), v1.end(), v2.begin(), v2.end(), std::inserter(result, result.end()));
};
}
#endif // EZ_STRINGUTIL_H
