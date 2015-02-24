/*
20111106 rsz Created.
*/
#ifndef EZNCVARIANT_H
#define EZNCVARIANT_H

#include "netcdf.h"

namespace ez {
class ezNcVariant {
// Wrap a single value and have it accessible for any cast type.
// Set its type before the value.
public:
  ezNcVariant() : type(NC_INT64), ul(0) {};
  // Get value auto-casted to desired output type.
  template <typename T>
  T get();
  // Set with type not necessarily same type as target type.
  template <typename T>
  void set(T v);
  
  // What type this really should be used as.
  nc_type type;
  // The actual value.
  union {
    char c;
    short s;
    int i;
    long long l;
    float f;
    double d;
    unsigned char uc;
    unsigned short us;
    unsigned int ui;
    unsigned long long ul;
  };
};
//#####################################################################
template <typename T>
T ezNcVariant::get() {
  switch(type) { 
    case NC_CHAR: return static_cast<T>(c); break;
    case NC_SHORT: return static_cast<T>(s); break;
    case NC_INT: return static_cast<T>(i); break;
    case NC_INT64: return static_cast<T>(l); break;
    case NC_FLOAT: return static_cast<T>(f); break;
    case NC_DOUBLE: return static_cast<T>(d); break;
    case NC_UBYTE: case NC_BYTE: return static_cast<T>(uc); break;
    case NC_USHORT: return static_cast<T>(us); break;
    case NC_UINT: return static_cast<T>(ui); break;
    case NC_UINT64: return static_cast<T>(ul); break;
    default: return 0; break;
  }
};
//#####################################################################
template <typename T>
void ezNcVariant::set(T v) {
  switch(type) { 
    case NC_CHAR: c = static_cast<char>(v); break;
    case NC_SHORT: s = static_cast<short>(v); break;
    case NC_INT: i = static_cast<int>(v); break;
    case NC_INT64: l = static_cast<long long>(v); break;
    case NC_FLOAT: f = static_cast<float>(v); break;
    case NC_DOUBLE: d = static_cast<double>(v); break;
    case NC_UBYTE: case NC_BYTE: uc = static_cast<unsigned char>(v); break;
    case NC_USHORT: us = static_cast<unsigned short>(v); break;
    case NC_UINT: ui = static_cast<unsigned int>(v); break;
    case NC_UINT64: ul = static_cast<unsigned long long>(v); break;
    default: break;    
  }
};
//#####################################################################
};
#endif // EZVARIANT_H